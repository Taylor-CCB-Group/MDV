import os
import sys
import h5py
import numpy
import pandas
import scipy
import json
import gzip
import shlex
import subprocess
import fasteners
import warnings
import shutil
import random
import string
from os.path import join, split, exists, basename
from pathlib import Path
import scipy.sparse
from werkzeug.utils import secure_filename
from shutil import copytree, ignore_patterns, copyfile
from typing import Optional, NewType, List, Union, Any, cast
from pandas.api.types import is_bool_dtype
import polars as pl
# from mdvtools.charts.view import View 
import time
import copy
import tempfile
from mdvtools.image_view_prototype import create_image_view_prototype
from mdvtools.charts.table_plot import TablePlot
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)

DataSourceName = str  # NewType("DataSourceName", str)
ColumnName = str  # NewType("ColumnName", str)
# List[ColumnName] gets tricky, `ColumnName | str` syntax needs python>=3.10
Cols = Union[List[str], NewType("Params", List[ColumnName])]

datatype_mappings = {
    "int64": "integer",
    "float64": "double",
    "float32": "double",
    "object": "text",
    "category": "text",
    "bool": "text",
    "int32": "double",
    "boolean":"text"
}

numpy_dtypes = {
    "text": numpy.ubyte,
    "text16": numpy.uint16,
    "multitext": numpy.uint16,
    "double": numpy.float32,
    "integer": numpy.float32,
    "int32": numpy.int32,
    # unique created in fly (depends on string length)
}


class MDVProject:
    def __init__(
        self,
        dir: str,
        id: Optional[str] = None,
        delete_existing=False,
        skip_column_clean=True, # todo - make this False by default, add tests etc
        backend_db = False,
        safe_file_save = True
    ):
        self.skip_column_clean = (
            skip_column_clean  # signficant speedup for large datasets
        )
        self.safe_file_save = safe_file_save
        self.dir = dir
        self.id = id if id else os.path.basename(dir)
        if delete_existing and exists(dir):
            shutil.rmtree(dir)
        self.h5file = join(dir, "datafile.h5")
        self.datasourcesfile = join(dir, "datasources.json")
        self.statefile = join(dir, "state.json")
        self.viewsfile = join(dir, "views.json")
        self.readmefile = join(dir, 'README.md')
        self.imagefolder = join(dir, "images")
        self.trackfolder = join(dir, "tracks")
        if not exists(dir):
            os.mkdir(dir)
        if not exists(self.trackfolder):
            os.mkdir(self.trackfolder)
        if not exists(self.datasourcesfile):
            with open(self.datasourcesfile, "w") as o:
                o.write(json.dumps([]))
        if not exists(self.viewsfile):
            with open(self.viewsfile, "w") as o:
                o.write(json.dumps({}))
                o.close()
        if not exists(self.statefile):
            with open(self.statefile, "w") as o:
                o.write(json.dumps({"all_views": []}))
        self._lock = fasteners.InterProcessReaderWriterLock(join(dir, "lock"))
        self.backend_db = backend_db

    @property
    def datasources(self):
        return get_json(self.datasourcesfile)

    @datasources.setter
    def datasources(self, value):
        save_json(self.datasourcesfile, value, self.safe_file_save)
    
    @property
    def readme(self):
        if not exists(self.readmefile): 
            return None
        with open(self.readmefile, 'r') as f:
            markdown_string = f.read()
        return markdown_string

    @property
    def views(self):
        return get_json(self.viewsfile)

    @views.setter
    def views(self, value):
        save_json(self.viewsfile, value, self.safe_file_save)

    @property
    def state(self):
        # check this - seeing errors in the logs?
        return get_json(self.statefile)

    @state.setter
    def state(self, value):
        save_json(self.statefile, value ,self.safe_file_save)

    def set_editable(self, edit=True):
        c = self.state
        c["permission"] = "edit" if edit else "view"
        self.state = c

    def set_chat_enabled(self, chat_enabled=True):
        # would prefer not to be adding methods in this file, maybe it could be in chat_server_extension.py
        c = self.state
        c["chat_enabled"] = chat_enabled
        self.state = c

    def lock(self, type="read"):
        return self._lock.read_lock() if type == "read" else self._lock.write_lock()

    def get_column_metadata(self, datasource, column):
        ds = self.get_datasource_metadata(datasource)
        col = [x for x in ds["columns"] if x["field"] == column]
        if len(col) == 0:
            raise AttributeError(
                f"column {column} not found in {datasource} datasource"
            )
        return col[0]

    def set_column_metadata(self, datasource, column, parameter, value):
        ds = self.get_datasource_metadata(datasource)
        col_index = [c for c, x in enumerate(ds["columns"]) if x["field"] == column]
        if len(col_index) == 0:
            raise AttributeError(
                f"column {column} not found in {datasource} datasource"
            )
        ds["columns"][col_index[0]][parameter] = value
        self.set_datasource_metadata(ds)

    def get_datasource_as_dataframe(self, datasource: str) -> pandas.DataFrame:
        ## nb, I'd quite like to have some 'DataSourceName' type alias so we know it's not just any string...
        ds = self.get_datasource_metadata(datasource)
        df = pandas.DataFrame()
        for c in ds["columns"]:
            data = self.get_column(datasource, c["field"])
            df[c["field"]] = data
        return df

    def check_columns_exist(self, datasource, columns):
        md = self.get_datasource_metadata(datasource)
        all_cols = set([x["field"] for x in md["columns"]])
        return [x for x in columns if x not in all_cols]

    def get_datasource_names(self) -> list[str]:
        """
        Get a list of all datasource names in the project.
        
        Returns:
            list[str]: A list of datasource names
        """
        return [ds["name"] for ds in self.datasources]
    
    def set_interactions(
        self,
        interaction_ds,
        parent_ds,
        pivot_column="sample_id",
        parent_column="annotation",
        is_single_region=True,
        interaction_columns=["Cell Type 1", "Cell Type 2"],
        default_parameter="Cross PCF gr20",
        node_size="cell 1 number",
        add_view=True,
    ):
        """
        See `spatialdata.md` for documentation...
        """
        # check columns exist in the appropriate data sets
        missing_cols = self.check_columns_exist(
            interaction_ds,
            [pivot_column, default_parameter, node_size] + interaction_columns,
        )
        if len(missing_cols) > 0:
            raise AttributeError(
                f'columns {",".join(missing_cols)} not found in {interaction_ds} datasource'
            )
        missing_cols = self.check_columns_exist(
            parent_ds, [pivot_column, parent_column]
        )
        if len(missing_cols) > 0:
            raise AttributeError(
                f'columns {",".join(missing_cols)} not found in {parent_ds} datasource'
            )
        # update the config
        md = self.get_datasource_metadata(interaction_ds)
        md["interactions"] = {
            "pivot_column": pivot_column,
            "is_single_region": is_single_region,
            "interaction_columns": interaction_columns,
            "spatial_connectivity_map": {
                "link_length": default_parameter,
                "link_thickness": default_parameter,
                "link_color": default_parameter,
                "node_size": node_size,
            },
            "cell_radial_chart": {"link_thickness": default_parameter},
        }
        self.set_datasource_metadata(md)
        # update the links between datasources
        self.insert_link(
            interaction_ds,
            parent_ds,
            "interactions",
            {
                "interaction_columns": interaction_columns + [parent_column],
                "pivot_column": pivot_column,
                "is_single_region": is_single_region,
            },
        )
        if add_view:
            # todo add stuff to the view
            self.set_view(
                interaction_ds, {"initialCharts": {parent_ds: [], interaction_ds: []}}
            )

    def get_datasource_metadata(self, name):
        ds = [x for x in self.datasources if x["name"] == name]
        if len(ds) == 0:
            raise AttributeError(f"{name} datasource not found")
        return ds[0]

    def set_datasource_metadata(self, ds):
        mds = self.datasources
        index = [c for c, x in enumerate(mds) if x["name"] == ds["name"]]
        if len(index) == 0:
            mds.append(ds)
        else:
            mds[index[0]] = ds
        self.datasources = mds

    def add_images_to_datasource(
        self,
        ds: str,
        image_folder: str,
        key_column: str,
        name="Images",
        image_type="png",
    ):
        """Adds images to a datasource.
        These will be copied as static assets with `convert_to_static_page()`,
        and served from their original location by `serve()`.
        Args:
            ds (str): The name of the datasource.
            image_folder (str): The folder containing the images.
            key_column (str): The name of the column in ds containing the image names.
            name (str, optional): The name of the image set. Defaults to "Images".
            image_type (str, optional): The image type. Defaults to "png".
        """
        ds_metadata = self.get_datasource_metadata(ds)
        if "images" not in ds_metadata:
            ds_metadata["images"] = {}
        image_meta = ds_metadata["images"]
        if name in image_meta:
            raise AttributeError(f"Image set {name} already exists in {ds} datasource")
        image_meta[name] = {
            "original_folder": image_folder,  # for convenience locally, maybe shouldn't be saved, may be useful
            "base_url": f"./images/{ds}/{name}/",
            "key_column": key_column,
            "type": image_type,
        }
        logger.info(f"Added image set {name} to {ds} datasource")
        self.set_datasource_metadata(ds_metadata)

    def add_or_update_image_datasource(self, tiff_metadata, datasource_name, file):
        """Add or update an image datasource in datasources.json
        returns the name of the added view so the user can navigate to it"""
        # todo consider moving this elsewhere & generally review the methods for adding spatial data
        # Load current datasources
        datasources = self.datasources
        # Check if the datasource exists
        datasource = next((ds for ds in datasources if ds["name"] == datasource_name), None)
        datasource_backup = None
        
        is_new_datasource = False

        target_folder = os.path.join(self.imagefolder, 'avivator')
        if not os.path.exists(target_folder):
            os.makedirs(target_folder)
        original_filename = file.filename
        upload_file_path = os.path.join(target_folder, original_filename)
        view_name = None
        try:
            
            # Step 1: Update or create datasource
            if datasource:
                 # Create a backup of the existing datasource before updating
                datasource_backup = datasource.copy()
                view_name = self.update_datasource_for_tiff(datasource, datasource_name, tiff_metadata, original_filename)
            else:
                is_new_datasource = True
                datasource_name = "default" if not datasource_name else datasource_name

                # Check if the default datasource exists
                datasource = next((ds for ds in datasources if ds["name"] == datasource_name), None)

                view_name = self.update_datasource_for_tiff(datasource, datasource_name, tiff_metadata, original_filename)
            
            # Step 2: Upload the image
            self.upload_image_file(file, upload_file_path)
            
            # Step 3: Add database entry (exception will propagate up if it fails)
            if self.backend_db:
                from mdvtools.dbutils.dbservice import ProjectService, FileService

                FileService.add_or_update_file_in_project(
                    file_name=file.filename,
                    file_path=upload_file_path,  # Adjust as necessary for actual file path
                    project_id=self.id
                )
                
                ProjectService.set_project_update_timestamp(self.id)
            # Print success message
            logger.info(f"Datasource '{datasource_name}' updated, TIFF file uploaded, and database entry created successfully.")

        except Exception as e:
            logger.error(f"Error in MDVProject.add_or_update_image_datasource: {e}")
            
            # Attempt rollback actions
            try:
                # Rollback the file upload
                if os.path.exists(upload_file_path):  # Check the existence of the file at the upload path
                    logger.info("Reverting file upload...")
                    self.delete_uploaded_image(upload_file_path) 
                
                # Rollback datasource creation if it was new
                if is_new_datasource and any(ds['name'] == datasource_name for ds in self.datasources):
                    logger.info("Reverting new datasource creation...")
                    self.datasources = [x for x in self.datasources if x["name"] != datasource_name]
                    #self.delete_datasource(datasource_name, False)
                elif datasource:
                    logger.info("Reverting datasource update...")
                    self.restore_datasource(datasource_backup)  # This method may need to be implemented for updates
            except Exception as rollback_error:
                logger.error(f"Error during rollback in MDVProject.add_or_update_image_datasource: {rollback_error}")
            
            # Re-raise the original exception for the caller to handle
            raise
        return view_name
    def upload_image_file(self, file, upload_file_path):
        """Upload the TIFF file to the imagefolder, saving it with the original filename."""
        try:
            # Save the file to the /images/avivator folder
            file.save(upload_file_path)
            logger.info(f"File uploaded successfully to {upload_file_path}")

        except Exception as e:
            logger.error(f"Error in MDVProject.upload_image_file: Failed to upload file to '{upload_file_path}': {e}")
            raise
    
    def delete_uploaded_image(self, file_path):
        """Delete the uploaded image file at the specified path."""
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                logger.info(f"Deleted uploaded image at: {file_path}")
            else:
                logger.info(f"File does not exist at: {file_path}")
        except Exception as e:
            logger.error(f"Error in MDVProject.delete_uploaded_image: Error deleting file at {file_path}: {e}")
            raise 
    
    def restore_datasource(self, datasource_backup):
        """Restore the datasource from the backup."""
        try:
            # Find the existing datasource by name
            existing_datasource = next(
                (ds for ds in self.datasources if ds["name"] == datasource_backup["name"]), 
                None
            )

            if existing_datasource:
                # Overwrite the existing datasource with the backup values
                existing_datasource.update(datasource_backup)  
                logger.info(f"Restored datasource '{datasource_backup['name']}' from backup. ")

                # Save the updated datasources to the JSON file
                self.datasources = self.datasources  # This will call the setter and save the data
            else:
                logger.warning(f"Warning: Could not find datasource '{datasource_backup['name']}' to restore.")
        
        except Exception as e:
            logger.error(f"Error in MDVProject.restore_datasource: {str(e)}")
            raise
    

    def update_datasource_for_tiff(self, datasource, datasource_name, tiff_metadata, region_name: str):
        """Update an existing datasource with new image metadata."""
        try:
            # Find the existing datasource by name
            existing_datasource = next((ds for ds in self.datasources if ds["name"] == datasource_name), None)

            logger.info("In update_datasource_for_tiff")
            logger.info(datasource_name)
            # print(existing_datasource)
            # print(datasource)
            #datasources empty template
            # If the datasource doesn't exist, create a new one
            if (existing_datasource is None):
                logger.info("In update_datasource_for_tiff: existing_datasource is None")
                datasource = self.create_datasource_template(datasource_name)
                self.datasources.append(datasource)  # Add the new datasource to the list
                logger.info(f"In MDVProject.update_datasource: Created new datasource template for '{datasource_name}'.")
                
                #adding default columns for new empty ds
                filename = 'mdvtools/dbutils/emptyds.csv'
                df_default = pandas.read_csv(filename)
                # self.add_datasource(project_id, datasource_name, df_default, add_to_view=None)
                self.add_datasource(datasource_name, df_default, add_to_view=None)
                datasource = self.get_datasource_metadata(datasource_name)
            
            
            #datasources-> regions section
            # Ensure datasource has a 'regions' field
            if "regions" not in datasource:
               datasource["regions"] = {}
            

            # Corrected path to access size and scale information
            pixels_data = tiff_metadata['Pixels']
            width = pixels_data['SizeX']
            height = pixels_data['SizeY']
            scale = pixels_data.get('PhysicalSizeX', 1.0)
            scale_unit = pixels_data.get('PhysicalSizeXUnit', 'µm')  # Default to µm if not present

            # Call ensure_regions_fields to ensure the required fields and values are set
            datasource['regions'] = self.ensure_regions_fields(
                datasource['regions'],  # Existing regions dictionary
                scale_unit=scale_unit,  # Pass the scale unit
                scale=scale             # Pass the scale
            )
            
            
            #datasources-> regions -> all_regions-> new entry section
            # Determine region name
            # full_name = tiff_metadata.get('Name', 'unknown')
            # region_name = full_name.split(".ome")[0] if ".ome" in full_name else full_name  # Use 'Name' from metadata or 'unknown'
            
            # Define new region with metadata
            new_region = {
                "roi": {
                    "min_x": 0,
                    "min_y": 0,
                    "max_x": width,
                    "max_y": height
                },
                "images": {},
                "viv_image": {
                    "file": region_name,
                    "linked_file": True
                }
            }

            # Update or add the region in the datasource
            datasource["regions"]["all_regions"][region_name] = new_region
            #datasource['size'] = len(datasource['regions']['all_regions'])

            # Save the updated datasource
            self.set_datasource_metadata(datasource)

            #add empty default columns
            
            # Update views and image            
            region_view_json = create_image_view_prototype(datasource_name, region_name)

            #views = self.views

            #if views and "default" in views:
            #    view_name = region_name  # If "default" exists, set view_name to region_name
            #else:
            #    view_name = "default"
            view_name = region_name
            logger.info(view_name)
            self.set_view(view_name, region_view_json)
            
            #self.add_viv_images(region_name, image_metadata, link_images=True)
            return view_name
        
        except Exception as e:
            logger.error(f"Error in MDVProject.update_datasource :  Error updating datasource '{datasource_name}': {e}")
            raise
    
    def create_datasource_template(self, datasource_name: str) -> dict:
        """Create a new datasource template with the basic structure."""
        try:
            template = {
                "name": datasource_name,
                "columns": [],          # Initialize empty columns
                "size": 0,              # Start size at 0
                "regions": {},            # Initialize regions as an empty dictionary,
                "columnGroups": []
            }
            logger.info(f"Created new datasource template '{datasource_name}'.")
            return template
        
        except Exception as e:
            logger.error(f"In MDVProject.create_datasource_template: Error creating datasource template '{datasource_name}': {e}")
            raise  # Re-raises the caught exception

    def ensure_regions_fields(self, 
                          regions, 
                          position_fields=['x', 'y'], 
                          region_field='sample_id', 
                          default_color='leiden', 
                          scale_unit='µm', 
                          scale=1.0):
        try:
            
            # Ensure that required fields exist with default values
            regions["position_fields"] = regions.get("position_fields", position_fields)
            regions["region_field"] = regions.get("region_field", region_field)
            regions["default_color"] = regions.get("default_color", default_color)
            regions["scale_unit"] = regions.get("scale_unit", scale_unit)
            regions["scale"] = regions.get("scale", scale)
            
            # Initialize regions structure if needed
            if "all_regions" not in regions:
                regions["all_regions"] = {}

            # Check if 'avivator' field is present, if not, add default settings
            if "avivator" not in regions:
                default_channels = [{'name': 'DAPI'}]
                regions["avivator"] = {
                    "default_channels": default_channels,
                    "base_url": "images/avivator/",
                }

            return regions
        
        except Exception as e:
            logger.error(f"In MDVProject.ensure_regions_fields: Error in ensure_regions_fields: {e}")
            raise  # Re-raises the caught exception

    
    def get_image(self, path: str):
        """Gets the filename of an image."""
        # assume path is of the form <ds>/<name>/<filename>
        p = path.split("/")
        ds = p[0]
        filename = p[-1]
        ds_metadata = self.get_datasource_metadata(ds)
        image_meta = ds_metadata["images"][p[1]]
        return join(image_meta["original_folder"], filename)

    def _get_h5_handle(self, read_only=False, attempt=0) -> h5py.File:
        mode = "r"
        if not exists(self.h5file):
            mode = "w"
        elif not read_only:
            mode = "a"
        try:
            return h5py.File(self.h5file, mode)
        except Exception:
            # certain environments seem to have issues with the handle not being closed instantly
            # if there is a better way to do this, please change it
            # if there are multiple processes trying to access the file, this may also help
            # (although if they're trying to write, who knows what bad things may happen to the project in general)
            time.sleep(0.1)
            attempt += 1
            logger.error(f"error opening h5 file, attempt {attempt}...")
            return self._get_h5_handle(read_only, attempt)

    def get_column(self, datasource: str, column, raw=False):
        cm = self.get_column_metadata(datasource, column)
        h5 = self._get_h5_handle()
        gr = h5[datasource]
        if not isinstance(gr, h5py.Group):
            h5.close()
            raise AttributeError(f"cannot open datasource '{datasource}' in h5 file")
        raw_data = numpy.array(gr[column])
        if raw:
            return raw_data
        dt = cm["datatype"]
        if dt == "text" or dt == "text16":
            data = [cm["values"][x] for x in raw_data]
        elif dt == "multitext":
            chunksize = raw_data.shape[0] / cm["stringLength"]
            arr = numpy.split(raw_data, chunksize)
            data = [",".join([cm["values"][x] for x in y if x != 65535]) for y in arr]
        elif dt == "unique":
            data = [x.decode() for x in raw_data]
        else:
            data = list(raw_data)
        h5.close()
        return data

    def set_column_with_raw_data(self, datasource, column, raw_data):
        """Adds or updates a column with raw data
        Args:
            datasource (str): The name of the datasource.
            column (dict): The complete metadata for the column
            raw_data (list|array): The raw binary data for the column
        """
        h5 = self._get_h5_handle()
        cid = column["field"]
        ds = h5[datasource]
        if not isinstance(ds, h5py.Group):
            raise AttributeError(f"{datasource} is not a group")
        if ds.get(cid):
            del ds[cid]
        dt = numpy_dtypes.get(column["datatype"])
        if not dt:
            dt = h5py.string_dtype("utf-8", column["stringLength"])
        ds.create_dataset(cid, len(raw_data), data=raw_data, dtype=dt)
        ds = self.get_datasource_metadata(datasource)
        cols = ds["columns"]
        ind = [c for c, x in enumerate(cols) if x["field"] == cid]
        if len(ind) == 0:
            cols.append(column)
        else:
            cols[ind[0]] = column
        self.set_datasource_metadata(ds)

    def set_column(self, datasource, column, data):
        """Adds (or replaces an existing column) with the data supplied

        Args:
            datasource (str): The name of the datasource.
            column (str|dict):  metadata for the column. Can be a string with the column's name,
                although datatype should also be included as the inferred datatype
                is not always correct
            raw_data (list|array): Anything that can be converted into a pandas Series
            The data should be in the correct order
        """
        if isinstance(column, str):
            column = {"name": column}
        if not column.get("field"):
            column["field"] = column["name"]
        ds = self.get_datasource_metadata(datasource)
        ind = [c for c, x in enumerate(ds["columns"]) if x["field"] == column["field"]]
        col_exists = len(ind) > 0
        li = pandas.Series(data)
        if not column.get("datatype"):
            column["datatype"] = datatype_mappings.get(str(li.dtype), "text")
        h5 = self._get_h5_handle()
        gr = h5[datasource]
        if not isinstance(gr, h5py.Group):
            raise AttributeError(f"{datasource} is not a group")
        if gr.get(column["field"]):
            del gr[column["field"]]
        add_column_to_group(column, li, gr, len(li), self.skip_column_clean)
        h5.close()
        if col_exists:
            ds["columns"][ind[0]] = column
        else:
            ds["columns"].append(column)
        self.set_datasource_metadata(ds)

    def remove_column(self, datasource, column):
        """Removes the specified column

        Args:
            datasource (str): The name of the data source.
            column (str): The id (field) of the column.
        """
        ds = self.get_datasource_metadata(datasource)
        cols = [x for x in ds["columns"] if x["field"] != column]
        if len(cols) == len(ds["columns"]):
            warnings.warn(f"deleting non existing column: {column} from {datasource}")
            return
        ds["columns"] = cols
        h5 = self._get_h5_handle()
        gr = h5[datasource]
        if not isinstance(gr, h5py.Group):
            raise AttributeError(f"{datasource} is not a group")
        del gr[column]
        self.set_datasource_metadata(ds)

    def add_annotations(
        self,
        datasource: str,
        data: pandas.DataFrame | str,
        separator="\t",
        missing_value="ND",
        columns=[],
        supplied_columns_only=False,
    ):
        """Adds annotations based on an existing column

        Args:
            datasource (str): The name of the data source.
            data (dataframe|str): Either a pandas dataframe or a text file. The first column
                should be the 'index' column and match a column in the datasource. The other columns should
                contain the annotations to add.
            separator (str,optional): The delimiter if a text file is supplied (tab by default)
            missing_value(str,optional): The value to put if the index value is missing in the input data.
                Default is 'ND'
        """
        if isinstance(data, str):
            data = pandas.read_csv(data, sep=separator)
        assert(isinstance(data, pandas.DataFrame))
        ds = self.get_datasource_metadata(datasource)
        index_col = data.columns[0]
        data = data.set_index(index_col)
        columns = get_column_info(columns, data, supplied_columns_only)
        col = [x for x in ds["columns"] if x["field"] == index_col]
        if len(col) == 0:
            raise AttributeError(
                f"index column {index_col} not found in {datasource} datasource"
            )
        # py-right: dictionary key must be hashable
        newdf = pandas.DataFrame({index_col: self.get_column(datasource, index_col)})  # type: ignore
        h5 = self._get_h5_handle()
        gr = h5[datasource]
        assert isinstance(gr, h5py.Group)
        for c in columns:
            d = {k: v for k, v in zip(data.index, data[c["field"]])}
            # v slow - needs improving
            # py-right: `Argument of type "Series | Unknown | DataFrame" cannot be assigned to parameter "data" of type "Series"`
            ncol = newdf.apply(lambda row: d.get(row[0], missing_value), axis=1)
            assert isinstance(ncol, pandas.Series)
            add_column_to_group(c, ncol, gr, len(ncol), self.skip_column_clean)
            ds["columns"].append(c)
        self.set_datasource_metadata(ds)
        h5.close()

    def set_column_group(self, datasource, groupname, columns):
        """Adds (or changes) a column group

        Args:
            datasource(string): The name of the datasource
            groupname(string): The name of the column group
            columns(list): The field names of columns in the group. If None, then the column
                group will be removed
        """
        ds = self.get_datasource_metadata(datasource)
        # check if columns exists
        if columns:
            colfields = set([x["field"] for x in ds["columns"]])
            missingcols = [x for x in columns if x not in colfields]
            if len(missingcols) > 0:
                raise AttributeError(f"adding non existent columns ({','.join(missingcols)}) to column group {groupname}\
                                    in datasource {datasource}")
        cg = ds.get("columnGroups")
        # create entry if absent
        if not cg:
            cg = []
            ds["columnGroups"] = cg
        # does group exist
        ind = [c for c, x in enumerate(cg) if x["name"] == groupname]
        # change (or delete) existing group
        if len(ind) == 1:
            if columns:
                cg[ind[0]]["columns"] = columns
            else:
                del cg[ind[0]]
        # add new group
        else:
            # no group to delete
            if not columns:
                raise AttributeError(f"removing non existent column group {groupname}\
                                    from datasource {datasource}")
            # add new group
            cg.append({"name": groupname, "columns": columns})
        self.set_datasource_metadata(ds)

    def delete_datasource(self, name, delete_views=True):
        print("in delete -1 ")
        h5 = self._get_h5_handle()
        del h5[name]
        h5.close()
        print("in delete -2 ")
        self.datasources = [x for x in self.datasources if x["name"] != name]
        print("in delete -3 ")
        # delete all views contining that datasource
        if delete_views:
            print("in delete -4 ")
            views = self.views
            for view in views:
                data = views[view]
                if data["initialCharts"].get(name):
                    self.set_view(view, None)

    def add_genome_browser(
        self,
        datasource,
        parameters=["chr", "start", "end"],
        name=None,
        extra_params=None,
        custom_track=None,
        overwrite=False,
    ):
        """
        args:
            datasource (string): The name of the datasource
            parameters (list, optional): The names of the columns containing the chromosome, start and end
                positions. Defaults to ["chr","start","end"]
            name (string, optional): The name of the genome browser track. Defaults to the datasource name.
        """
        try:
            check_htslib()  # will raise an error if htslib is not installed
        except Exception:
            raise Exception(
                "htslib not installed. This is not supported on Windows, other platforms will need to install e.g. via brew install htslib"
            )
        if len(parameters) != 3:
            raise AttributeError(
                "genome browser parameters should be a list of 3 column names"
            )
        if not name:
            name = datasource
        ds = self.get_datasource_metadata(datasource)
        if "genome_browser" in ds and not overwrite:
            raise AttributeError(
                f"genome browser track already exists for {datasource}"
            )
        track_name = f"{secure_filename(name)}.bed"
        if not custom_track:
            # get all the genome locations
            loc = [self.get_column(datasource, x) for x in parameters]
            # write to a bed file
            # nb - should check whether it actually improves anything adding {name} to filenames
            # bed = join(self.trackfolder, f"t_{name}.bed") # reverting to Martin's original
            bed = join(self.trackfolder, "t.bed")
            o = open(bed, "w")
            for c, (chr, start, end) in enumerate(zip(loc[0], loc[1], loc[2])):
                o.write(f"{chr}\t{start}\t{end}\t{c}\n")
            o.close()
            # indexed_bed = join(self.trackfolder, "loc_{name}.bed") # reverting to Martin's original
            indexed_bed = join(self.trackfolder, track_name)
            create_bed_gz_file(bed, indexed_bed)
            os.remove(bed)
        else:
            # copy the custom track to the tracks folder
            shutil.copy(custom_track["location"], join(self.trackfolder, f"{track_name}.gz"))
            # copy index file
            shutil.copy(
                custom_track["location"] + ".tbi",
                join(self.trackfolder, f"{track_name}.gz.tbi"),
            )

        if not name:
            name = datasource
        gb = {
            "location_fields": parameters,
            # "default_track": {"url": "tracks/loc_{name}.bed.gz", "label": name}, # reverting to Martin's original
            "default_track": {"url": f"tracks/{track_name}.gz", "label": name},
        }
        if custom_track:
            gb["default_track"]["type"] = custom_track["type"]
        if extra_params:
            gb.update(extra_params)
        ds = self.get_datasource_metadata(datasource)
        ds["genome_browser"] = gb
        self.set_datasource_metadata(ds)

    def get_genome_browser(self, datasource):
        ds = self.get_datasource_metadata(datasource)
        info = ds.get("genome_browser")
        if not info:
            raise AttributeError(f"no genome browser for {datasource}")
        default_track =  {
            "short_label": info["default_track"]["label"],
            "url": info["default_track"]["url"],
            "track_id": "_base_track",
            "decode_function": "generic",
            "height": 15,
            "displayMode": "EXPANDED"
        }
        if info.get("default_track_parameters"):
            default_track.update(info["default_track_parameters"])
        
        gb = {
            "type": "genome_browser",
            "param": info["location_fields"],
            "tracks": [default_track]
        }
        at = info.get("atac_bam_track")
        if at:
            gb["tracks"].append(
                {
                    "short_label": "Coverage",
                    "height": 400,
                    "track_id": "_atac_bam_track",
                    "url": at["url"],
                    "type": "bam_sca_track",
                }
            )
        dt = info.get("default_tracks")
        if dt:
            for t in dt:
                gb["tracks"].append(t)
        if info["default_parameters"]:
            gb.update(info["default_parameters"])
        return gb

    def add_refseq_track(self, datasource, genome):
        ds = self.get_datasource_metadata(datasource)
        gb = ds.get("genome_browser")
        if not gb:
            raise AttributeError(f"no genome browser for {datasource}")
        tdir = join(split(os.path.abspath(__file__))[0], "templates", "tracks")
        reft = join(tdir, f"{genome}.bed.gz")
        if not os.path.exists(reft):
            raise AttributeError(f"no refseq track for {genome}")
        dt = gb.get("default_tracks")
        if not dt:
            dt = gb["default_tracks"] = []
        # add to start of list
        dt.insert(
            0,
            {
                "short_label": "RefSeq",
                "height": 50,
                "displayMode": "EXPANDED",
                "decode_function": "decodeRefflat",
                "track_id": "_refseq_track",
                "url": f"tracks/{genome}.bed.gz",
            },
        )
        # copy to tracks folder
        shutil.copy(reft, join(self.trackfolder, f"{genome}.bed.gz"))
        shutil.copy(reft + ".tbi", join(self.trackfolder, f"{genome}.bed.gz.tbi"))
        self.set_datasource_metadata(ds)

    def add_tracks(self,datasource: str,tracks: list[dict]) :
        """Adds a list of tracks to the datasource's genome browser.
        
        Args:
            tracks (list[dict]): A list of track dictionaries to add.
            datasource (str): The name of the datasource to which the tracks will be added.
        """
        ds = self.get_datasource_metadata(datasource)
        gb = ds.get("genome_browser")
        if not gb:
            raise AttributeError(f"no genome browser for {datasource}")
        dt = gb.get("default_tracks", [])
        for track in tracks:
            if not isinstance(track, dict):
                raise TypeError("Each track must be a dictionary")
            if  "file" not in track:
                raise ValueError("Each track must have specify a local or remote file")
            fname = basename(track["file"])
            track_name = track.get("name", fname.split(".")[0])
            track_type= track.get("type")
            if not track_type:
                if ".bed" in fname or fname.endswith(".bb"):
                    track_type = "bed"
                elif ".bigwig" in fname or fname.endswith(".bw"):
                    track_type = "wig"
            if not track_type:
                raise AttributeError(f" the type of track {fname} cannot be deduced")
            #no need to do anything - will be served from the original location
            if track["file"].startswith("http"):
                url = track["file"]
            else:
                if not exists(track["file"]):
                    raise FileNotFoundError(f"Track file {track['file']} does not exist")
                # tracks already compressed and indexed - just copy to tracks folder
                if track_type == "wig" or fname.endswith(".gz") or fname.endswith(".bb"):
                    to_file = join(self.trackfolder, fname)
                    shutil.copy(track["file"], to_file)
                    # for .gz also need to copy index
                    if fname.endswith(".gz"):
                        i_file = track["file"] + ".tbi"
                        if not exists(i_file):
                            raise FileNotFoundError(f"Index file {i_file} not found")
                        shutil.copyfile(i_file, f"{to_file}.tbi")
                # assume its a just a bed file- compress and index it
                else:
                    t_file = join(track["file"])
                    o_file = join(self.trackfolder, fname)
                    create_bed_gz_file(t_file, o_file)        
                    fname= fname+".gz"
                url = f"./tracks/{fname}"
            #will need to adapt this for other browsers
            mtrack ={
                "short_label": track_name,
                "url": url,
                "track_id": track.get("id",track_name),
                "color": track.get("color", "black"),
            }
            if track_type== "bed":
                mtrack["featureHeight"] = track.get("featureHeight", 10)
                mtrack["height"] = mtrack["featureHeight"] + 12
                mtrack["displayMode"] = track.get("displayMode", "EXPANDED")
            elif track_type == "wig":
                mtrack["height"] = track.get("height", 50)
            for param in ["hideLabels"]:
                if track.get(param):
                    mtrack[param] = track[param]
            dt.append(mtrack)
        gb["default_tracks"] = dt
        self.set_datasource_metadata(ds)

    def add_datasource(
        self,
        # project_id: str,
        name: str,
        dataframe: pandas.DataFrame | str, # could we add xarray here - we pass things from anndata which may not be DataFrame.
        columns: Optional[list] = None,
        supplied_columns_only=False,
        replace_data=False,
        add_to_view: Optional[str] = "default",
        separator="\t"
    ) -> list[dict[str, str]]:
        """Adds a pandas dataframe to the project. Each column's datatype, will be deduced by the
        data it contains, but this is not always accurate. Hence, you can supply a list of column
        metadata, which will override the names/types deduced from the dataframe.
        
        Args:
            name (string): The name of datasource
            dataframe (dataframe|str): Either a pandas dataframe or the path of a text file
            columns (list, optional) : A list of objects containing the column name and datatype.
            supplied_columns_only (bool, optional): If True, only the the subset of columns in the columns argument
            replace_data (bool, optional): If True, the existing datasource will be overwritten, Default is False,
            add_to_view (string, optional): The datasource will be added to the specified view.
            separator (str, optional): If a path to text file is supplied, then this should be the file's delimiter.
        """
        dodgy_columns = []  # To hold any columns that can't be added
        gr = None  # Initialize the group variable
        h5 = None
        
        try:
            if isinstance(dataframe, str):
                dataframe = pandas.read_csv(dataframe, sep=separator)
            
            # Get columns to add
            columns = get_column_info(columns, dataframe, supplied_columns_only)
            
            # Check if the datasource already exists
            try:
                ds = self.get_datasource_metadata(name)
            except Exception:
                ds = None

            if ds:
                # Delete the existing datasource if replace_data is True
                if replace_data:
                    self.delete_datasource(name)
                else:
                    raise FileExistsError(
                        f"Attempt to create datasource '{name}' failed because it already exists."
                    )
            
            # Open HDF5 file and handle group creation
            try:
                h5 = self._get_h5_handle()
                
                # Check for and delete existing group with this name
                if name in h5:
                    del h5[name]
                    logger.warning(f"Deleted existing group '{name}' in HDF5 file.")
                
                
                gr = h5.create_group(name)
            except Exception as e:
                raise RuntimeError(f"Error managing HDF5 groups for datasource '{name}': {e}")
            
            # Verify columns are provided
            if not columns:
                raise AttributeError("No columns to add. Please provide valid columns metadata.")
            
            # Add columns to the HDF5 group
            dodgy_columns = []
            for col in columns:
                try:
                    add_column_to_group(col, dataframe[col["field"]], gr, len(dataframe), self.skip_column_clean) # type: ignore
                except Exception as e:
                    dodgy_columns.append(col["field"])
                    logger.warning(
                        f"Failed to add column '{col['field']}' to datasource '{name}': {repr(e)}"
                    )
            
            h5.close()  # Close HDF5 file
            columns = [x for x in columns if x["field"] not in dodgy_columns]
            #print(f" - non-dodgy columns: {columns}")
            
            # Update datasource metadata
            ds = {"name": name, "columns": columns, "size": len(dataframe)}
            #print(f'--- setting datasource metadata: {ds}')
            self.set_datasource_metadata(ds)
            
            # Add to view if specified
            if add_to_view:
                # TablePlot parameters
                title=name,
                #params = ["leiden", "ARVCF", "DOK3", "FAM210B", "GBGT1", "NFE2L2", "UBE2D4", "YPEL2"]
                params = dataframe.columns.to_list()
                size = [792, 472]
                position = [10, 10]
            
                # Create plot
                table_plot = self.create_table_plot(title, params, size, position)
                
                # Convert plot to JSON and set view
                table_plot_json = self.convert_plot_to_json(table_plot)
                
                v = self.get_view(add_to_view)
                if not v:
                    v = {"initialCharts": {}}

                # Check if the initialCharts already has entries for `name`
                if name not in v["initialCharts"] or not v["initialCharts"][name]:
                    # If empty, initialize with [table_plot_json]
                    v["initialCharts"][name] = [table_plot_json]
                else:
                    # If not empty, append table_plot_json to the existing list
                    v["initialCharts"][name].append(table_plot_json)
                    
                self.set_view(add_to_view, v)
            
            
            # Update the project's update timestamp using the dedicated method
            if self.backend_db:
                from mdvtools.dbutils.dbservice import ProjectService
                ProjectService.set_project_update_timestamp(self.id)
            
            
            logger.info(f"add_datasource: Added datasource successfully '{name}'")
            return dodgy_columns

        except Exception as e:
            logger.error(f"Error in MDVProject.add_datasource : Error adding datasource '{name}': {e}")
            raise  # Re-raise the exception to propagate it to the caller

    def add_datasource_polars(
        self,
        name: str,
        dataframe: pl.DataFrame | pl.LazyFrame | str,
        columns: Optional[list] = None,
        supplied_columns_only=False,
        replace_data=False,
        add_to_view: Optional[str] = "default",
        separator="\t"
    ) -> list[dict[str, str]]:
        """Adds a polars dataframe to the project. Each column's datatype will be deduced by the
        data it contains, but this is not always accurate. Hence, you can supply a list of column
        metadata, which will override the names/types deduced from the dataframe.
        
        Args:
            name (string): The name of datasource
            dataframe (pl.DataFrame|str): Either a polars dataframe or the path of a text file
            columns (list, optional) : A list of objects containing the column name and datatype.
            supplied_columns_only (bool, optional): If True, only the subset of columns in the columns argument
            replace_data (bool, optional): If True, the existing datasource will be overwritten, Default is False,
            add_to_view (string, optional): The datasource will be added to the specified view.
            separator (str, optional): If a path to text file is supplied, then this should be the file's delimiter.
        """
        dodgy_columns = []  # To hold any columns that can't be added
        gr = None  # Initialize the group variable
        h5 = None
        
        try:
            print("starting add_datasource_polars")
            
            # If a path is provided, use scan_csv for lazy loading to prevent
            # loading the entire file into memory at once.
            if isinstance(dataframe, str):
                dataframe = pl.scan_csv(dataframe, separator=separator, try_parse_dates=False)

            is_lazy = isinstance(dataframe, pl.LazyFrame)

            # Determine the number of rows. This is memory-efficient even for lazy frames.
            if isinstance(dataframe, pl.LazyFrame):
                num_rows = dataframe.select(pl.count()).collect().item()
            else:
                num_rows = len(dataframe)

            # Get columns to add using polars-compatible function.
            # Pass num_rows to avoid recalculating it for LazyFrames.
            columns = get_column_info_polars(columns, dataframe, supplied_columns_only, num_rows)
            
            # Check if the datasource already exists
            try:
                ds = self.get_datasource_metadata(name)
            except Exception:
                ds = None

            print(f"is ds None? {ds}")

            if ds:
                # Delete the existing datasource if replace_data is True
                if replace_data:
                    self.delete_datasource(name)
                else:
                    raise FileExistsError(
                        f"Attempt to create datasource '{name}' failed because it already exists."
                    )
            
            print("got passed the ds check")
            
            # Open HDF5 file and handle group creation
            try:
                h5 = self._get_h5_handle()
                
                # Print current groups for visibility
                for group_name in h5.keys():
                    print(group_name)
                    
                # Check for and delete existing group with this name
                if name in h5:
                    del h5[name]
                    print(f"Deleted existing group '{name}' in HDF5 file.")
                
                gr = h5.create_group(name)
            except Exception as e:
                raise RuntimeError(f"Error managing HDF5 groups for datasource '{name}': {e}") from e
            
            print("created h5 group without error")
            
            # Verify columns are provided
            if not columns:
                raise AttributeError("No columns to add. Please provide valid columns metadata.")
            
            # Add columns to the HDF5 group
            dodgy_columns = []
            for col in columns:
                try:
                    print(f"- adding column '{col['field']}' to datasource '{name}'")
                    
                    # If working with a lazy frame, select and collect one column at a time.
                    # This reads only one column from the file into memory.
                    if is_lazy:
                        polars_series = dataframe.select(col["field"]).collect().get_columns()[0]
                    else:
                        polars_series = dataframe[col["field"]]

                    add_column_to_group_from_polars(col, polars_series, gr, num_rows, self.skip_column_clean)
                except Exception as e:
                    print(f" ++++++ DODGY COLUMN: {col['field']}")
                    dodgy_columns.append(col["field"])
                    warnings.warn(
                        f"Failed to add column '{col['field']}' to datasource '{name}': {repr(e)}"
                    )
            
            h5.close()  # Close HDF5 file
            columns = [x for x in columns if x["field"] not in dodgy_columns]
            
            # Update datasource metadata
            ds = {"name": name, "columns": columns, "size": num_rows}
            self.set_datasource_metadata(ds)
            
            print("Updated datasource metadata")
            
            # Add to view if specified
            if add_to_view:
                # TablePlot parameters
                title = name
                # Use collect_schema().names() for LazyFrame, .columns for DataFrame
                if isinstance(dataframe, pl.LazyFrame):
                    params = dataframe.collect_schema().names()
                else:
                    params = dataframe.columns
                size = [792, 472]
                position = [10, 10]
            
                # Create plot
                table_plot = self.create_table_plot(title, params, size, position)
                
                # Convert plot to JSON and set view
                table_plot_json = self.convert_plot_to_json(table_plot)
                
                v = self.get_view(add_to_view)
                if not v:
                    v = {"initialCharts": {}}

                # Check if the initialCharts already has entries for `name`
                if name not in v["initialCharts"] or not v["initialCharts"][name]:
                    # If empty, initialize with [table_plot_json]
                    v["initialCharts"][name] = [table_plot_json]
                else:
                    # If not empty, append table_plot_json to the existing list
                    v["initialCharts"][name].append(table_plot_json)
                    
                self.set_view(add_to_view, v)
            
            # Update the project's update timestamp using the dedicated method
            if self.backend_db:
                from mdvtools.dbutils.dbservice import ProjectService
                ProjectService.set_project_update_timestamp(self.id)
            
            print(f"In MDVProject.add_datasource_polars: Added datasource successfully '{name}'")
            return dodgy_columns

        except Exception as e:
            print(f"Error in MDVProject.add_datasource_polars : Error adding datasource '{name}': {e}")
            raise  # Re-raise the exception to propagate it to the caller


    def create_table_plot(self, title, params, size, position):
        """Create and configure a TablePlot instance with the given parameters."""
        plot = TablePlot(
            title=title,
            params=params,
            size=size,
            position=position
        )
        
        return plot
    
    def convert_plot_to_json(self, plot):
        """Convert plot data to JSON format."""
        return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

    def insert_link(self, datasource: str, linkto: str, linktype: str, data):
        """
        Adds a link between two datasources.
        datasource is the name of the datasource to which the link will be added.
        linkto is the name of the datasource to which the link will point.
        linktype is the type of link to add.
        The data argument is a dictionary containing the data for the link,
        the format of which depends on the linktype.
        """
        ds = self.get_datasource_metadata(datasource)
        if not ds:
            raise AttributeError(f"datasource '{datasource}' does not exist")
        to_ds = self.get_datasource_metadata(linkto)
        if not to_ds:
            raise AttributeError(f"datasource for linkto '{linkto}' does not exist")
        links = ds.get("links")
        if not links:
            links = {}
            ds["links"] = links
        llink = links.get(linkto)
        if not llink:
            llink = {}
            links[linkto] = llink
        llink[linktype] = data
        self.set_datasource_metadata(ds)

    def add_rows_as_columns_link(
        self, ds_row: str, ds_col: str, column_name: str, name: str
    ):
        """
        Adds a link between two datasources, such that columns may be added to the `ds_row` datasource
        based on the values in the `ds_col` datasource dynamically at runtime. The values in the `column_name` column of the
        `ds_col` datasource will be used as the names for the columns added to the `ds_row` datasource.

        Args:
            ds_row (string): The name of the datasource into which the link will be added
            ds_col (string): The name of the datasource containing the columns
            column_name (string): The name of a column in the `ds_col` datasource, the row-value of which will be used
                as the column name for columns dynamically added to the `ds_row` datasource
            name (string): The name of the link that will appear in the interface to describe the nature of the data being added, e.g. `'Gene Expr'`
        """
        to_ds = self.get_datasource_metadata(ds_col)
        if not to_ds:
            raise AttributeError(f"datasource for ds_col '{ds_col}' does not exist")
        if column_name not in [c["name"] for c in to_ds["columns"]]:
            raise AttributeError(
                f"column '{column_name}' does not exist in datasource '{ds_col}'"
            )
        data = {"name_column": column_name, "name": name, "subgroups": {}}
        self.insert_link(ds_row, ds_col, "rows_as_columns", data)

    def add_rows_as_columns_subgroup(
        self,
        row_ds: str,
        col_ds: str,
        stub: str,
        data,
        name: Optional[str] = None,
        label: Optional[str] = None,
        sparse: Optional[bool] = None, # this should be inferred from the data
        chunk_data=False
    ):
        """Add rows as columns in a subgroup."""
        name = name if name else stub
        label = label if label else name
        h5 = self._get_h5_handle()
        ds = h5[row_ds]
        if not isinstance(ds, h5py.Group):
            raise AttributeError(f"{row_ds} is not a group")
        if name in ds:
            raise ValueError(f"Group '{name}' already exists in {row_ds}.")
        
        gr = ds.create_group(name)
        
        if sparse is None:
            # Infer if the data is sparse or dense unless specified
            sparse = scipy.sparse.issparse(data)

        if sparse:
            # Handle sparse matrix
            density = f"{len(data.data) / (data.shape[0] * data.shape[1]):.3f}"
            logger.info(f"Adding sparse matrix to {row_ds} with {len(data.data)} elements (density {density})")
            gr.create_dataset(
                "x", (len(data.data),), data=data.data, dtype=numpy.float32
            )
            gr.create_dataset(
                "i", (len(data.indices),), data=data.indices, dtype=numpy.uint32
            )
            gr.create_dataset("p", (len(data.indptr),), data=data.indptr)
        else:
            # Fallback to dense or convertible
            logger.info(f"Adding dense matrix to {row_ds}")
            try:
                #Requires a lot of RAM leads to memory errors on smaller machines
                if not chunk_data:
                    dense_data = data.toarray() if hasattr(data, 'toarray') else numpy.asarray(data)
                    total_len = dense_data.shape[0] * dense_data.shape[1]
                    gr.create_dataset(
                        "x", (total_len,),
                        data=dense_data.flatten(order="F"),
                        dtype=numpy.float32
                    )
                    gr["length"] = [dense_data.shape[0]]
                else:
                    # Create the dataset with chunking and compression - slow
                    num_rows, num_cols = data.shape
                    chunk_size = num_cols
                    total_len = num_rows * num_cols
                
                    dset = gr.create_dataset(
                        "x", (total_len,),
                        dtype=numpy.float32,
                        chunks=(chunk_size,),
                        compression="gzip"
                    )
                    gr.create_dataset("length", data=[num_rows])

                    # Process the data in chunks to avoid high memory usage
                    # A flat single array is not the ideal storage but in theory 
                    # only small arrays will be dense. However large arrays are common due to
                    # normalization creating (different) non zero values in each gene for 0 reads
                    # Transposing the array without loading into memmory would probably take the 
                    # same amount of time
                    # Adapt chunk_size so large row length won't consume too much memory.
                    chunk_size = self._get_optimal_chunk_size(num_rows, num_cols)
                    for i in range(0, num_cols,chunk_size):
                        # protect against out-of-bounds end_col here 
                        # - although in testing, no difference observed, because
                        # total_len of dset is the right size & both slices will truncate correctly
                        end_col = min(i + chunk_size, num_cols)
                        dset[i*num_rows:(end_col)*num_rows] = data[:,i:end_col].flatten("F")
            except Exception as e:
                raise TypeError(f"Unsupported data type for dense processing: {type(data)}. Original error: {e}")

        # Update metadata
        ds = self.get_datasource_metadata(row_ds)
        ds["links"][col_ds]["rows_as_columns"]["subgroups"][stub] = {
            "name": name,
            "label": label,
            "type": "sparse" if sparse else "dense",
        }
        self.set_datasource_metadata(ds)
        h5.close()

    def _get_optimal_chunk_size(self, rows: int, cols: int) -> int:
        """Determine optimal chunk size based on matrix dimensions and available memory."""
        try:
            import psutil
            available_memory_mb = psutil.virtual_memory().available / 1024 / 1024
        except ImportError:
            # Fallback if psutil is not available
            available_memory_mb = 1000  # Assume 1GB available
        
        # Estimate memory per element (float32 = 4 bytes)
        bytes_per_element = 4
        
        # Calculate how many elements we can process with available memory
        # Use 25% of available memory to be very safe for large datasets
        safe_memory_bytes = available_memory_mb * 1024 * 1024 * 0.25
        max_elements = safe_memory_bytes / bytes_per_element
        
        if rows * cols <= max_elements:
            return cols
        else:
            return max(1, int(max_elements / rows))

    def get_links(self, datasource, filter=None):
        ds = self.get_datasource_metadata(datasource)
        links = []
        lnks = ds.get("links")
        if lnks:
            for lnkto in lnks:
                lnk = lnks[lnkto]
                if filter is None or lnk.get(filter):
                    links.append({"datasource": lnkto, "link": lnk})
        return links

    def serve(self, **kwargs):
        from mdvtools.server import create_app
        from mdvtools.server_extension import MDVServerOptions

        options = kwargs.pop("options", None)

        if options is not None:
            if kwargs:
                warnings.warn(
                    "Redundant keyword arguments passed to serve() along with 'options' object. "
                    "These arguments will be ignored.",
                    UserWarning
                )
        else:
            # `options` was not provided, create it from kwargs.
            options = MDVServerOptions(**kwargs)
        
        create_app(self, options=options)
        

    def delete(self):
        # todo - remove from project routes, set a flag indicating it's been deleted
        shutil.rmtree(self.dir)

    def get_configs(self):
        config = {
            "datasources": self.datasources,
            "state": self.state,
        }
        
        # legacy
        hyperion_conf = join(self.dir, "hyperion_config.json")
        if os.path.exists(hyperion_conf):
            config["hyperion_config"] = get_json(hyperion_conf)
        # end
        return config

    def convert_to_static_page(self, outdir, include_sab_headers=True):
        fdir = split(os.path.abspath(__file__))[0]
        # consider adding an option to do a js build here
        tdir = join(fdir, "templates")
        # copy everything except the data (todo options for images etc)
        copytree(self.dir, outdir, ignore=ignore_patterns(*("*.h5", "*.ome.tiff")))
        # copy the js and images
        self.copy_images(outdir)
        copytree(join(fdir, "static"), join(outdir, "static"))
        # create the static binary files
        self.convert_data_to_binary(outdir)
        # write out the index file
        page = "static_page.html"
        template = join(tdir, page)
        page = open(template).read()
        # make sure the static files are referenced correctly
        page = page.replace("/static", "static")
        # call init with no route, will be interpreted as static page (at /)
        page = page.replace(
            "_mdvInit('{{route}}')", "_mdvInit()"
        )  # not relevant for build-flask-vite build
        # correct config
        conf = self.state
        # can't edit static page
        conf["permission"] = "view"
        # consider using this flag for determining front-end behaviour
        conf["static"] = True
        # throttle the dataloading so don't get network errors
        conf["dataloading"] = {"split": 5, "threads": 2}
        save_json(join(outdir, "state.json"), conf,self.safe_file_save)
        # add service worker for cross origin headers
        if include_sab_headers:
            page = page.replace("<!--sw-->", '<script src="serviceworker.js"></script>')
            copyfile(join(tdir, "serviceworker.js"), join(outdir, "serviceworker.js"))
        with open(join(outdir, "index.html"), "w") as o:
            o.write(page)

    def copy_images(self, outdir):
        for ds in self.datasources:
            name = ds["name"]
            if "images" not in ds:
                # print("no images for", name)
                continue
            image_metadata = ds["images"]
            dsdir = join(outdir, "images", name)
            if not os.path.exists(dsdir):
                os.makedirs(dsdir)
            for image_set_name in image_metadata:
                image_set = image_metadata[image_set_name]
                print("copying", image_set_name)
                original_folder = image_set["original_folder"]
                copytree(original_folder, join(outdir, image_set["base_url"]))
            # TODO also copy any linked avivator images

    def save_state(self, state):
        # update/add or view
        # view will be deleted if view is null
        if state.get("currentView"):
            self.set_view(state["currentView"], state["view"])
        ud = state.get("updatedColumns")
        # update/add/delete any columns
        if ud:
            for ds in ud:
                item = ud[ds]
                for data in item["colors_changed"]:
                    self.set_column_metadata(
                        ds, data["column"], "colors", data["colors"]
                    )
                for data in item["columns"]:
                    self.set_column_with_raw_data(ds, data["metadata"], data["data"])
                for col in item["removed"]:
                    self.remove_column(ds, col)
        # update any datasource metadata
        md = state.get("metadata")
        if md:
            for ds in md:
                datasource = self.get_datasource_metadata(ds)
                for param in md[ds]:
                    datasource[param] = md[ds][param]
                self.set_datasource_metadata(datasource)

    def add_image_set(self, datasource, setname, column, folder, type="png",large=False):
        """Adds a set of images to a datasource. The images should be in a folder, with the same name as the column
        Args:
            datasource (str): The name of the datasource.
            column (str): The name of the column describing the images. The folder
                should contain images with the same name as the values in this column
                minus the file extension (typically .png or .jpg).
            type (str, optional): The type of the images. Default is 
                'png'. Other options are 'jpg', 'jpeg', etc.
            large (bool, optional): If True, the images will be avialable to the Row Summmary Box
                where only a single image is shown and can be panned and zoomed
                if false (default) the images will be available in the Image Table and should
                be thumbnails of the same size
            setname (str): The name of the image set. More than one set of images can be associated
                with a datasource. The name should be unique within the datasource and is used tp 
                create a folder (with a sanitized name) in the images directory     
            folder (str): The path to the folder containing the images.
        """
        ds = self.get_datasource_metadata(datasource)
        # col = self.get_column_metadata(datasource, column)

        images = [x for x in os.listdir(folder) if x.endswith(type)]
        # create the image folder
        fname = secure_filename(setname)
        imdir = join(self.imagefolder, fname)
        if not exists(imdir):
            os.makedirs(imdir)
        # copy the images
        for im in images:
            copyfile(join(folder, im), join(imdir, im))
        # update the metadata
        imt = "large_images" if large else "images"
        if not ds.get(imt):
            ds[imt] = {}
        ds[imt][setname] = {
            "key_column": column,
            "type": type,
            "base_url": f"./images/{fname}/",
        }

        self.set_datasource_metadata(ds)

    def get_view(self, view):
        views = self.views
        return views.get(view)

    def set_view(self, name: str, view: Optional[Any], make_default=False):
        """Sets the view with the given name to the supplied view data.
        If the view is None, then the view will be deleted.
        """
        views = self.views
        # update or add the view
        if view:
            # clone, in case a calling script might get confused by the view being changed
            view2 = copy.deepcopy(view)
            # generate ids for the views if not present
            # may consider more robust id generation & other checks here
            initialCharts = view2["initialCharts"]
            for ds in initialCharts:
                # todo check that there is a ds with that name
                for c in initialCharts[ds]:
                    if "id" not in c:
                        c["id"] = str(
                            "".join(random.choices(string.ascii_letters, k=6))
                        )
            views[name] = view2
        # remove the view
        else:
            if views.get(name):
                del views[name]
        self.views = views

        state = self.state
        # add to list and make default
        if view:
            if name not in state["all_views"]:
                state["all_views"].append(name)
            if make_default:
                state["initial_view"] = name
        # delete from list
        else:
            state["all_views"].remove(name)
            iv = state.get("initial_view")
            # if the deleted view is the default view then
            # change the default view to the first view in the list
            if iv:
                state["initial_view"] = state["all_views"][0]
        self.state = state

    def convert_data_to_binary(self, outdir=None):
        if not outdir:
            outdir = self.dir
        h5 = h5py.File(self.h5file)
        dss = self.datasources
        for ds in dss:
            n = ds["name"]
            gr = h5[n]
            if not isinstance(gr, h5py.Group):
                raise TypeError("Expected 'gr' to be of type h5py.Group.")
            dfile = join(outdir, "{}.gz".format(n))
            o = open(dfile, "wb")
            index = {}
            current_pos = 0
            for c in ds["columns"]:
                dt = gr.get(c["field"])
                if not dt:
                    continue
                arr = numpy.array(dt)
                comp = gzip.compress(arr.tobytes())
                o.write(comp)
                new_pos = current_pos + len(comp)
                index[c["field"]] = [current_pos, new_pos - 1]
                current_pos = new_pos

            # add rows to columns gene score / tensors etc
            lnks = self.get_links(n, "rows_as_columns")
            for ln in lnks:
                rc = ln["link"]["rows_as_columns"]
                for sg in rc["subgroups"]:
                    info = rc["subgroups"][sg]
                    sgrp = gr[info["name"]]
                    sparse = info.get("type") == "sparse"
                    # get number of rows in linked datasource
                    plen = [x["size"] for x in dss if x["name"] == ln["datasource"]][0]
                    for i in range(0, plen):
                        comp = gzip.compress(get_subgroup_bytes(sgrp, i, sparse))
                        o.write(comp)
                        new_pos = current_pos + len(comp)
                        index[f"{sg}{i}"] = [current_pos, new_pos - 1]
                        current_pos = new_pos

            o.close()
            ifile = dfile[: dfile.rindex(".")] + ".json"
            i = open(ifile, "w")
            i.write(json.dumps(index, allow_nan=False))
            i.close()

    def get_byte_data(self, columns, group):
        h5 = h5py.File(self.h5file, "r")
        byte_list = []
        gr = h5[group]
        if not isinstance(gr, h5py.Group):
            raise TypeError("Expected 'gr' to be of type h5py.Group.")
        for column in columns:
            sg = column.get("subgroup")
            if sg:
                sgindex = int(column["sgindex"])
                byte_list.append(
                    get_subgroup_bytes(
                        gr[sg], sgindex, column.get("sgtype") == "sparse"
                    )
                )
            else:
                data = gr[column["field"]]
                byte_list.append(numpy.array(data).tobytes())
        h5.close()
        return b"".join(byte_list)

    def set_region_data(
        self,
        datasource,
        data,
        region_field="sample_id",
        default_color="annotations",
        position_fields=["x", "y"],
        scale_unit="µm",
        scale=1.0,
    ):
        md = self.get_datasource_metadata(datasource)
        cols = set([x["field"] for x in md["columns"]])
        missing = [
            x
            for x in [region_field] + [default_color] + position_fields
            if x not in cols
        ]
        if len(missing) > 0:
            raise AttributeError(
                f"setting region data on {datasource} but the specified columns({','.join(missing)}) are missing"
            )
        md["regions"] = {
            "position_fields": position_fields,
            "region_field": region_field,
            "default_color": default_color,
            "scale_unit": scale_unit,
            "scale": scale,
        }
        # convert to dict
        if not isinstance(data, dict):
            df = pandas.read_csv(data, sep="\t")
            df.set_index(df.columns[0], inplace=True)
            data = df.to_dict("index")
        all_regions = {}
        for k, v in data.items():
            x = v.get("x_offset", 0)
            y = v.get("y_offset", 0)
            all_regions[k] = {
                "roi": {
                    "min_x": x,
                    "min_y": y,
                    "max_y": v["height"] + y,
                    "max_x": v["width"] + x,
                },
                # "images": v.get("images", {}),
                "images": {},
            }
            if "json" in v:
                f = v["json"]
                assert isinstance(f, str)  # in future we may allow dict or list
                if exists(f):
                    # copy the json file to the project
                    name = os.path.basename(f)
                    rel = join(self.dir, "json", name)
                    Path(rel).parent.mkdir(parents=True, exist_ok=True)
                    try:
                        shutil.copyfile(f, rel)
                        all_regions[k]["json"] = join("json", name)
                    except Exception as e:
                        logger.warning(
                            f"Skipping json for region {k} because of error copying {f} to {rel}\n{repr(e)}"
                        )
                else:
                    raise FileNotFoundError(f"json file '{f}' not found")
        # maybe warn if replacing existing regions
        # or add to existing regions
        md["regions"]["all_regions"] = all_regions
        self.set_datasource_metadata(md)

    def add_region_images(self, datasource: DataSourceName, data):
        """Adds images to regions in a datasource.

        Args:
            datasource (str): The name of the datasource.
            data (dict|str): A dictionary containing data about which images should be associated with
        """
        imdir = join(self.dir, "images", "regions")
        if not exists(imdir):
            os.makedirs(imdir)
        md = self.get_datasource_metadata(datasource)

        md["regions"]["base_url"] = "images/regions/"
        # convert flat file to dict
        if not isinstance(data, dict):
            df = pandas.read_csv(data, sep="\t")
            df.set_index(df.columns[0], inplace=True)
            data = df.to_dict("index")
        all_regions = md["regions"]["all_regions"]
        for k, v in data.items():
            region = all_regions.get(k)
            if not region:
                raise AttributeError(
                    f"adding image to non existant region ({k}) in {datasource}"
                )
            roi = region["roi"]
            name = v.get("name")
            x = v.get("offset_x", roi["min_x"])
            y = v.get("offset_y", roi["min_y"])
            region["default_image"] = name
            reg = {
                "position": [x, y],
                "height": v.get("height", roi["max_y"] - roi["min_y"]),
                "width": v.get("width", roi["max_x"] - roi["min_x"]),
                "name": name,
            }
            # simple url
            if v["path"].startswith("http"):
                reg["url"] = v["path"]
            # local file - need to copy to images directory
            else:
                im = split(v["path"])[1]
                im_details = im.split(".")
                newname = get_random_string() + "." + im_details[1]
                shutil.copyfile(v["path"], join(imdir, newname))
                reg["file"] = newname
            all_regions[k]["images"][name] = reg
        self.set_datasource_metadata(md)

    def add_viv_viewer(self, datasource_name, default_channels):
        """Add a Viv viewer to the specified datasource with default channels."""
        try:
            md = self.get_datasource_metadata(datasource_name)
            reg = md.get("regions")
            if not reg:
                raise AttributeError(
                    f"Adding viv viewer to {datasource_name}, which does not contain regions"
                )

            # Create the directory for storing images if it doesn't exist
            imdir = join(self.dir, "images", "avivator")
            if not exists(imdir):
                os.makedirs(imdir)

            # Add avivator configuration to the regions
            reg["avivator"] = {
                "default_channels": default_channels,
                "base_url": "images/avivator/",
            }
            
            # Save the updated metadata
            self.set_datasource_metadata(md)

        except Exception as e:
            logger.error(f"Error in MDVProject.add_viv_viewer: {e}")
            raise  # Re-raise the exception after logging


    def add_viv_images(self, datasource, data, link_images=True):
        md = self.get_datasource_metadata(datasource)
        try:
            md["regions"]["avivator"]
        except Exception:
            raise AttributeError(
                "Adding viv images when viv viewer has not been specified"
            )
        all_regions = md["regions"]["all_regions"]
        imdir = join(self.dir, "images", "avivator")
        if not isinstance(data, dict):
            df = pandas.read_csv(data, sep="\t")
            df.set_index(df.columns[0], inplace=True)
            data = df.to_dict("index")
        for k, v in data.items():
            region = all_regions.get(k)
            if not region:
                raise AttributeError(
                    f"adding image to non existant region ({k}) in {datasource}"
                )
            if v["path"].startswith("http"):
                region["viv_image"] = {"url": v["path"]}
            # local file - default is to copy to images directory
            # or link if link_images is True
            # todo consider situations where the image is in a different volume etc and may not be linkable.
            else:
                if not link_images:
                    newname = get_random_string() + ".ome.tiff"
                    shutil.copyfile(v["path"], join(imdir, newname))
                    region["viv_image"] = {"file": newname}
                else:
                    try:
                        dest = join(imdir, os.path.basename(v["path"]))
                        # we might want to use the same image in multiple regions / datasources
                        # in which case, it will already be linked
                        # this assumes that having the same name means it is actually the same image...
                        if not os.path.exists(dest):
                            os.symlink(
                                os.path.abspath(v["path"]),
                                join(imdir, os.path.basename(v["path"])),
                            )
                    except Exception as e:
                        logger.error(
                            f"Cannot link viv image '{v['path']}' to '{datasource}'\n{repr(e)}"
                        )
                    region["viv_image"] = {
                        "file": os.path.basename(v["path"]),
                        "linked_file": True,
                        # "original_folder": os.path.dirname(v["path"])
                    }
        self.set_datasource_metadata(md)

    def get_interaction_matrix(
        self, datasource, group, interaction_metric, square_size=20
    ):
        """
        Args:
            datasource (str): The name of the datasource.
            group (str): The name of the group.
            interaction_metric (str): The name of the interaction metric.
        """
        md = self.get_datasource_metadata(datasource)
        im_info = self.get_column_metadata(datasource, interaction_metric)
        i = md.get("interactions")
        if not i:
            raise AttributeError(f"no interactions in {datasource}")
        icd = self.get_column_metadata(datasource, i["interaction_columns"][0])[
            "values"
        ]
        side = len(icd) * square_size
        chart = {
            "type": "single_heat_map",
            "param": i["interaction_columns"]
            + [interaction_metric]
            + [i["pivot_column"]],
            "category": group,
            "title": f"{group} - {im_info['name']}",
            "size": [side, side],
            "axis": {
                "x": {
                    "textSize": 13,
                    "label": "",
                    "size": 101,
                    "tickfont": 11,
                    "rotate_labels": True,
                },
                "y": {"textSize": 13, "label": "", "size": 94, "tickfont": 10},
            },
        }
        return chart

    def get_selection_dialog(self, datasource, selections):
        filters = {}
        for s in selections:
            sel = s.get("filter")
            if sel:
                col = self.get_column_metadata(datasource, s["column"])
                if col["datatype"] not in ["text", "text16", "multitext"]:
                    if sel[0] is None:
                        sel[0] = col["minMax"][0]
                    if sel[1] is None:
                        sel[1] = col["minMax"][1]
                else:
                    if isinstance(sel, list):
                        sel = {"category": sel}
                filters[s["column"]] = sel
        return {
            "type": "selection_dialog",
            "param": [x["column"] for x in selections],
            "filters": filters,
        }

    def get_image_plot(self, datsource, image_set):
        md = self.get_datasource_metadata(datsource)
        ims = md.get("images")
        if not ims:
            raise AttributeError(f"no images in {datsource}")
        img = ims.get(image_set)
        if not img:
            raise AttributeError(f"no image set {image_set} in {datsource}")

        return {
            "type": "image_table_chart",
            "title": image_set,
            "param": [img["key_column"]],
            "image_set": image_set,
        }

    def get_centroid_plot(
        self, datasource, region, background_image="_default", scale=0.5
    ):
        """
        Args:
            datasource (str): The name of the datasource.
            region (str): The name of the region.
            background_image (str, optional): The name of the background image. Default is '_default'
            scale (float, optional): The scale of the image. Default is 0.5

        Returns:
            dict: The chart specification.
        """
        md = self.get_datasource_metadata(datasource)
        regions = md.get("regions")
        if not regions:
            raise AttributeError("no regions in specifeid")
        r_info = regions["all_regions"].get(region)
        if not r_info:
            raise AttributeError(f"no region {region} in regions")
        chart = {
            "type": "image_scatter_plot",
            "param": regions["position_fields"] + [regions["default_color"]],
            "background_filter": {
                "column": regions["region_field"],
                "category": region,
            },
            "title": region,
            "radius": 3.5,
            "color_by": regions["default_color"],
            "color_legend": {"dsiplay": False},
            "region": region,
            "roi": r_info.get("roi"),
        }
        dims = r_info["roi"]
        mx = dims.get("min_x", 0)
        my = dims.get("min_y", 0)
        size = [dims["max_x"] - mx, dims["max_y"] - my]
        size = [x * scale for x in size]
        chart["size"] = size
        if background_image:
            if background_image == "_default":
                background_image = r_info["default_image"]
            chart["background_image"] = r_info["images"][background_image]

        return chart

def get_json(file):
    with open(file) as f:
        return json.load(f)


def save_json(file, data, safe= True):
    try:
        if safe:
            save_json_atomic(file, data)
        else:
            with open(file,"w") as f:
                json.dump(data,f,indent=2,allow_nan=False)
    except Exception as e:
        logger.error(
            f"Error saving json to '{file}': some data cleaning may be necessary... project likely to be in a bad state."
        )
        raise (e)

def save_json_atomic(path, data):
    """
    Save JSON data to a file atomically.
    Hopefully this will be safer - in certain situations we were ending up with truncated output files.
    This method was suggested by ChatGPT: https://chatgpt.com/share/6813337b-9acc-800b-a6cd-6d058f339cd5
    """
    dir_name = os.path.dirname(path)
    with tempfile.NamedTemporaryFile("w", dir=dir_name, delete=False) as tmp:
        json.dump(data, tmp, indent=2, allow_nan=False)
        tmp.flush()
        os.fsync(tmp.fileno())
        temp_name = tmp.name
    # potential issues particularly in Docker where the files are on a different volume
    # safest option is to sync like this before and after, but may be overkill
    # 'sync' will fail silently on windows
    os.system("sync")
    os.replace(temp_name, path)  # Atomic move on most OSes
    os.system("sync")
    #file is saved with restictive permissions (this may be intentional)
    #
    # this method is lower overhead than os.system("sync") 
    # but stress testing indicates it is less robust
    # dir_fd = os.open(dir_name, os.O_DIRECTORY)
    # os.fsync(dir_fd)
    # os.close(dir_fd)

def get_subgroup_bytes(grp, index, sparse=False):
    if sparse:
        offset = grp["p"][index : index + 2]
        _len = offset[1] - offset[0]
        _indexes = numpy.array(grp["i"][offset[0] : offset[1]])
        _values = numpy.array(grp["x"][offset[0] : offset[1]], numpy.float32)
        return (
            numpy.array([_len], numpy.uint32).tobytes()
            + numpy.array(_indexes).tobytes()
            + numpy.array(_values).tobytes()
        )
    else:
        _len = grp["length"][0]
        offset = index * _len
        return numpy.array(grp["x"][offset : offset + _len], numpy.float32).tobytes()


def add_column_to_group(
    col: dict,
    data: pandas.Series,
    group: h5py.Group,
    length: int,
    skip_column_clean: bool,
):
    """
    col (dict): The column metadata (may be modified e.g. to add values)
    data (pandas.Series): The data to add
    group (h5py.Group): The group to add the data to
    length (int): The length of the data
    """

    if (
        col["datatype"] == "text"
        or col["datatype"] == "unique"
        or col["datatype"] == "text16"
    ):
        #in pandas missing values are represented by NaN
        #which cause problems when co-ercing into text, therefore replace with ND
        if isinstance(data.dtype, pandas.CategoricalDtype):
            # Handle pandas Categorical data
            if "ND" not in data.cat.categories:
                # see test_categorical_missing_values_edge_cases()
                data = data.cat.add_categories("ND")
            data = data.fillna("ND")
        #no boolean datatype in MDV at the moment, have to co-erce to text
        elif is_bool_dtype(data):
            cat = data.apply(lambda x: "True" if x is True else "False" if x is False else "ND")
            assert isinstance(cat, pandas.Series)
            data = cat
        else:
            # may need to double-check this...
            data = data.fillna("ND")
        values = data.value_counts()

        if len(values) < 65537 and col["datatype"] != "unique":
            t8 = len(values) < 257
            col["datatype"] = "text" if t8 else "text16"
            dtype = numpy.ubyte if t8 else numpy.uint16
            if not col.get("values"):
                col["values"] = [x for x in values.index if values[x] != 0]
            vdict = {k: v for v, k in enumerate(col["values"])}
            group.create_dataset(
                col["field"],
                length,
                dtype=dtype,
                data=data.map(vdict),  # type: ignore
            )
            # convert to string
            col["values"] = [str(x) for x in col["values"]]

        else:
            max_len = max(data.str.len())
            utf8_type = h5py.string_dtype("utf-8", int(max_len))
            col["datatype"] = "unique"
            col["stringLength"] = max_len
            group.create_dataset(col["field"], length, data=data, dtype=utf8_type)

    elif col["datatype"] == "multitext":
        delim = col.get("delimiter", ",")
        value_set: set[str] = set()
        maxv = 0
        # first parse - get all possible values and max number
        # of values in a single field
        for v in data:
            if not isinstance(v, str):
                continue
            vs = v.split(delim)
            value_set.update([x.strip() for x in vs])
            maxv = max(maxv, len(vs))
        if "" in value_set:
            value_set.remove("")
        ndata = numpy.empty(shape=(length * maxv,), dtype=numpy.uint16)
        ndata.fill(65535)
        values = list(value_set)
        # dict more efficient than index list
        vmap = {k: v for v, k in enumerate(values)}
        for i in range(0, length):
            b = i * maxv
            try:
                v = data[i]  # may raise KeyError if data is None at this index
                if not isinstance(v, str) or v == "":
                    continue
                vs = v.split(delim)
                vs = [x.strip() for x in vs]
            except Exception:
                continue
            vs.sort()
            for n in range(0, len(vs)):
                ndata[b + n] = vmap[vs[n]]
        col["values"] = values
        col["stringLength"] = maxv
        group.create_dataset(
            col["field"], length * maxv, data=ndata, dtype=numpy.uint16
        )
    else:
        dt = numpy.int32 if col["datatype"] == "int32" else numpy.float32
        clean = (
            data
            if skip_column_clean
            #this is pretty fast now
            else
                pandas.to_numeric(
                    data.iloc[:, 0] if isinstance(data, pandas.DataFrame) else data,
                    errors="coerce"
                )
        )  # this is slooooow?
        # faster but non=numeric values have to be certain values
        # clean=data.replace("?",numpy.NaN).replace("ND",numpy.NaN).replace("None",numpy.NaN)
        ds = group.create_dataset(col["field"], length, data=clean, dtype=dt)
        # remove NaNs for min/max and quantiles - this needs to be tested with 'inf' as well.
        na = numpy.array(ds)
        na = na[numpy.isfinite(na)]
        col["minMax"] = [float(str(numpy.amin(na))), float(str(numpy.amax(na)))]
        quantiles = [0.001, 0.01, 0.05]
        col["quantiles"] = {}
        for q in quantiles:
            #quantiles must be of type float else won't serialise to json
            col["quantiles"][str(q)] = [
                float(numpy.percentile(na, 100 * q)),
                float(numpy.percentile(na, 100 * (1 - q))),
            ]

def get_column_info(columns, dataframe, supplied_columns_only):
    if columns:
        for col in columns:
            if not col.get("field"):
                col["field"] = col["name"]

    if not supplied_columns_only:
        cols = [
            {"datatype": datatype_mappings[d.name], "name": c, "field": c}
            for d, c in zip(dataframe.dtypes, dataframe.columns)
        ]
        # replace with user given column metadata
        if columns:
            col_map = {x["field"]: x for x in columns}
            cols = [col_map.get(x["field"], x) for x in cols]
        columns = cols
    return columns


def check_htslib():
    try:
        subprocess.run(["tabix", "--version"])
    except Exception:
        raise AttributeError(
            "htslib not found, needed for preparing genome browser data"
        )


##!! will not work in windows and requires htslib installed
def create_bed_gz_file(infile, outfile):
    # need to sort
    command = "sort -k1,1V -k2,2n -k3,3n {} > {}".format(
        shlex.quote(infile), shlex.quote(outfile)
    )
    os.system(command)
    subprocess.run(["bgzip", outfile])
    subprocess.run(["tabix", outfile + ".gz"])


def get_random_string(length=6):
    return "".join(
        random.choices(
            string.ascii_uppercase + string.ascii_lowercase + string.digits, k=length
        )
    )


def map_polars_to_mdv_type(polars_dtype):
    """Map Polars data types to MDV project data types"""
    import polars as pl
    
    # Map Polars types to our internal types
    if polars_dtype in [pl.Int8, pl.Int16, pl.Int32, pl.Int64, pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64]:
        return "integer"
    elif polars_dtype in [pl.Float32, pl.Float64]:
        return "double"
    elif polars_dtype == pl.Boolean:
        return "text"  # Map boolean to text as per original mapping
    elif polars_dtype == pl.Categorical:
        return "text"
    elif polars_dtype in [pl.Utf8, pl.String]:
        return "text"  # Default to text, might change to unique based on cardinality
    else:
        # Default to text for any other types
        return "text"

def add_column_to_group_from_polars(col_info, polars_series, h5_group, num_rows, skip_column_clean):
    """Add a column to HDF5 group using Polars Series data"""
    import numpy as np
    
    field = col_info["field"]
    
    # Handle different column types
    if col_info["datatype"] in ["text", "text16", "multitext"]:
        # For categorical-like columns
        if col_info["datatype"] == "multitext":
            # This requires special handling for delimited strings
            # First, extract unique values across all rows
            if isinstance(polars_series, pl.Series):
                # Convert to pandas for compatible processing with existing code
                data = polars_series.to_pandas()
                # Use existing add_column_to_group since it handles complex multitext logic
                add_column_to_group(col_info, data, h5_group, num_rows, skip_column_clean)
                return
        else:
            # For regular categorical columns
            unique_values = polars_series.unique().sort()
            # Only use text16 if we have too many categories for text
            if len(unique_values) >= 256:
                col_info["datatype"] = "text16"
                dtype = np.uint16
            else:
                col_info["datatype"] = "text"
                dtype = np.ubyte
            
            # Store the unique values in the column metadata
            col_info["values"] = [str(v) for v in unique_values.to_list()]
            
            # Create mapping from values to indices
            value_to_index = {str(v): i for i, v in enumerate(col_info["values"])}
            
            # Create an empty dataset and write to it in chunks
            dset = h5_group.create_dataset(field, (num_rows,), dtype=dtype)
            chunk_size = 100_000  # Process 100,000 rows at a time

            for i in range(0, num_rows, chunk_size):
                chunk = polars_series[i : i + chunk_size]
                encoded_chunk = np.array([value_to_index.get(str(v), 0) for v in chunk], dtype=dtype)
                dset[i : i + len(encoded_chunk)] = encoded_chunk
    
    elif col_info["datatype"] == "unique":
        # For string columns with many unique values
        if isinstance(polars_series, pl.Series):
            # Calculate max length more efficiently
            max_len = polars_series.str.len_bytes().max()
            if max_len is None:
                max_len = 1
            
            # Create a fixed-length string dataset (HDF5 requirement)
            dt = h5py.string_dtype("utf-8", max_len)
            dset = h5_group.create_dataset(field, (num_rows,), dtype=dt)
            
            # Update metadata
            col_info["stringLength"] = max_len
            
            # Process and write in chunks
            chunk_size = 100_000  # Process 100,000 rows at a time
            for i in range(0, num_rows, chunk_size):
                chunk = polars_series[i : i + chunk_size]
                processed_chunk = chunk.fill_null("").to_list()
                dset[i : i + len(processed_chunk)] = processed_chunk
    
    else:
        # For numeric columns (integer, double)
        if isinstance(polars_series, pl.Series):
            # Use to_numpy() for more direct and memory-efficient conversion
            if col_info["datatype"] == "int32":
                data = polars_series.to_numpy(zero_copy_only=False).astype(np.int32)
            else:
                # fillna is necessary because NaN cannot be cast to int
                data = polars_series.fill_null(np.nan).to_numpy(zero_copy_only=False).astype(np.float32)

            # Create dataset
            h5_group.create_dataset(field, data=data)

            # Calculate stats for metadata
            finite_mask = np.isfinite(data)
            if np.any(finite_mask):
                finite_data = data[finite_mask]
                col_info["minMax"] = [float(np.min(finite_data)), float(np.max(finite_data))]
                
                col_info["quantiles"] = {}
                for q in [0.001, 0.01, 0.05]:
                    col_info["quantiles"][str(q)] = [
                        float(np.percentile(finite_data, 100 * q)),
                        float(np.percentile(finite_data, 100 * (1 - q)))
                    ]

def get_column_info_polars(columns, dataframe: "pl.DataFrame | pl.LazyFrame", supplied_columns_only, num_rows: int):
    """Polars version of get_column_info function, supports LazyFrames."""
    
    if columns:
        for col in columns:
            if not col.get("field"):
                col["field"] = col["name"]

    if not supplied_columns_only:
        # Get column info from polars dataframe or lazyframe
        cols = []
        # .schema works on both DataFrame and LazyFrame without loading data
        for col_name, polars_dtype in dataframe.collect_schema().items():
            mdv_dtype = map_polars_to_mdv_type(polars_dtype)
            
            # Check if this should be unique based on cardinality
            if mdv_dtype == "text":
                # Handle both DataFrame and LazyFrame cases
                if isinstance(dataframe, pl.LazyFrame):
                    unique_count = dataframe.select(pl.col(col_name).n_unique()).collect().item()
                else:
                    # For DataFrame, don't call collect()
                    unique_count = dataframe.select(pl.col(col_name).n_unique()).item()
                total_count = num_rows
                # If most values are unique, treat as unique type
                if unique_count > min(65536, total_count * 0.8):
                    mdv_dtype = "unique"
            
            cols.append({
                "datatype": mdv_dtype,
                "name": col_name,
                "field": col_name
            })
        
        # Replace with user given column metadata
        if columns:
            col_map = {x["field"]: x for x in columns}
            cols = [col_map.get(x["field"], x) for x in cols]
        columns = cols
    
    return columns

if __name__ == "__main__":
    path = os.getcwd()
    if len(sys.argv) > 1:
        path = sys.argv[1]
    if os.path.exists(os.path.join(path, "datasources.json")):
        MDVProject(path).serve()
