import os
import sys
import h5py
import numpy
import pandas
import json
import gzip
import shlex
import subprocess
import fasteners
import warnings
import shutil
import random
import string
from os.path import join, split, exists
from werkzeug.utils import secure_filename
from shutil import copytree, ignore_patterns, copyfile
from typing import Optional, NewType, List, Union
import time

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
}

numpy_dtypes = {
    "text": numpy.ubyte,
    "multitext": numpy.uint16,
    "double": numpy.float32,
    "integer": numpy.float32,
    "int32": numpy.int32,
    # unique created in fly (depends on string length)
}


class MDVProject:
    def __init__(self, dir: str, delete_existing=False):
        self.dir = dir
        self.name = dir.split("/")[-1]
        if delete_existing and exists(dir):
            shutil.rmtree(dir)
        self.h5file = join(dir, "datafile.h5")
        self.datasourcesfile = join(dir, "datasources.json")
        self.statefile = join(dir, "state.json")
        self.viewsfile = join(dir, "views.json")
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
                o.write(json.dumps({"all_views": [], "popouturl": "popout.html"}))
        self._lock = fasteners.InterProcessReaderWriterLock(join(dir, "lock"))

    @property
    def datasources(self):
        return get_json(self.datasourcesfile)

    @datasources.setter
    def datasources(self, value):
        save_json(self.datasourcesfile, value)

    @property
    def views(self):
        return get_json(self.viewsfile)

    @views.setter
    def views(self, value):
        save_json(self.viewsfile, value)

    @property
    def state(self):
        return get_json(self.statefile)

    @state.setter
    def state(self, value):
        save_json(self.statefile, value)

    def set_editable(self, edit):
        c = self.state
        c["permission"] = "edit" if edit else "view"
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
            df[c["name"]] = data
        return df

    def check_columns_exist(self, datasource, columns):
        md = self.get_datasource_metadata(datasource)
        all_cols = set([x["field"] for x in md["columns"]])
        return [x for x in columns if x not in all_cols]

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
        print(f"Added image set {name} to {ds} datasource")
        self.set_datasource_metadata(ds_metadata)

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
            print(f"error opening h5 file, attempt {attempt}...")
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
        if dt == "text":
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
        add_column_to_group(column, li, gr, len(li))
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
            add_column_to_group(c, ncol, gr, len(ncol))
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
        h5 = self._get_h5_handle()
        del h5[name]
        h5.close()
        self.datasources = [x for x in self.datasources if x["name"] != name]
        # delete all views contining that datasource
        if delete_views:
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
        # get all the genome locations
        loc = [self.get_column(datasource, x) for x in parameters]
        # write to a bed file
        # nb - should check whether it actually improves anything adding {name} to filenames
        bed = join(self.trackfolder, f"t_{name}.bed")
        o = open(bed, "w")
        for c, (chr, start, end) in enumerate(zip(loc[0], loc[1], loc[2])):
            o.write(f"{chr}\t{start}\t{end}\t{c}\n")
        o.close()
        indexed_bed = join(self.trackfolder, f"loc_{name}.bed")
        create_bed_gz_file(bed, indexed_bed)
        os.remove(bed)
        gb = {
            "location_fields": parameters,
            "default_track": {"url": f"tracks/loc_{name}.bed.gz", "label": name},
        }
        if extra_params:
            gb.update(extra_params)
        ds = self.get_datasource_metadata(datasource)
        ds["genome_browser"] = gb
        self.set_datasource_metadata(ds)

    def get_genome_browser(self, datasource):
        ds = self.get_datasource_metadata(datasource)
        info = ds.get("genome_browser")
        gb = {
            "type": "genome_browser",
            "param": info["location_fields"],
            "tracks": [
                {
                    "short_label": info["default_track"]["label"],
                    "url": info["default_track"]["url"],
                    "track_id": "_base_track",
                    "decode_function": "generic",
                    "height": 15,
                    "displayMode": "EXPANDED",
                }
            ],
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

    def add_datasource(
        self,
        name: str,
        dataframe: pandas.DataFrame | str,
        columns: Optional[list] = None,
        supplied_columns_only=False,
        replace_data=False,
        add_to_view: Optional[str] = "default",
        separator="\t",
    ) -> list[dict[str, str]]:
        """Adds a pandas dataframe to the project. Each column's datatype, will be deduced by the
        data it contains, but this is not always accurate. Hence, you can supply a list of column
        metadata, which will override the names/types deduced from the dataframe.

        Args:
            name (string): The name of datasource
            dataframe (dataframe|str): Either a pandas dataframe or the path of a text file
            columns (list, optional) : A list of objects containing the column name and datatype.
                e.g. [{"name":"column_1","datatype":"double"},]. If you want the column to have a
                different label, the object requires a field (the column name in the dataframe) and
                a name (the label seen by the user) e.g. {"field":"column_1","datatype":"double","name":"My Column 1"}
                In the case of "multitext" columns, you can also supply a "separator" field, otherwise it will
                default to a comma. <<< check this >>>
            supplied_columns_only(bool, optional): If True, only the the subset of columns in the columns argument
                will be added to the datasource. Default is False
            replace_data(bool, optional): If True, the existing datasource will be overwritten, Default is False,
                in which case, trying to add a datasource which already exists, will throw an error.
            add_to_view (string, optional): The datasource will be added to the specified view. The view will
                be created if it does not exist. The default is 'default'. If None, then it will not be added to
                a view.
            separator (str, optional): If a path to text file is supplied, then this should be the file's delimiter.
                Defaults to a tab.
        """
        if isinstance(dataframe, str):
            dataframe = pandas.read_csv(dataframe, sep=separator)
        # get the columns to add
        columns = get_column_info(columns, dataframe, supplied_columns_only)
        # does the datasource exist
        try:
            ds = self.get_datasource_metadata(name)
        except Exception:
            ds = None
        if ds:
            # delete the datasource
            if replace_data:
                self.delete_datasource(name)
            else:
                raise FileExistsError(
                    f"Trying to create {name} datasource, which already exits"
                )
        # create the h5 group
        h5 = self._get_h5_handle()
        gr = h5.create_group(name)
        size = len(dataframe)
        dodgy_columns = []
        if not columns:
            # we could set columns to an empty list, but it's probably better to throw an error
            # seems unlikely a user would want to add a datasource with no columns
            raise AttributeError("no columns to add")
        for col in columns:
            try:
                add_column_to_group(col, dataframe[col["field"]], gr, size)
            except Exception as e:
                dodgy_columns.append(col["field"])
                warnings.warn(
                    f"cannot add column '{col['field']}' to datasource '{name}':\n{repr(e)}"
                )

        h5.close()
        columns = [x for x in columns if x["field"] not in dodgy_columns]
        # add the metadata
        ds = None
        ds = {"name": name, "columns": columns, "size": size}
        self.set_datasource_metadata(ds)
        # add it to the view
        if add_to_view:
            v = self.get_view(add_to_view)
            if not v:
                v = {"initialCharts": {}}
            v["initialCharts"][name] = []
            self.set_view(add_to_view, v)
        return dodgy_columns

    def insert_link(self, datasource, linkto, linktype, data):
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

    def add_rows_as_columns_link(self, ds_row, ds_col, column_name, name):
        """
        Adds a link between two datasources, such that columns may be added to the `ds_row` datasource
        based on the values in the `ds_col` datasource dynamically at runtime. The values in the `column_name` column of the
        `ds_col` datasource will be used as the names for the columns added to the `ds_row` datasource.

        Args:
            ds_row (string): The name of the datasource into which the link will be added
            ds_col (string): The name of the datasource containing the columns
            column_name (string): The name of a column in the `ds_col` datasource, the row-value of which will be used
                as the column name for columns dynamically added to the `ds_row` datasource
            name (string): The name of the link
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
        self, row_ds, col_ds, stub, data, name=None, label=None, sparse=False
    ):
        name = name if name else stub
        label = label if label else name
        h5 = self._get_h5_handle()
        ds = h5[row_ds]
        if not isinstance(ds, h5py.Group):
            raise AttributeError(f"{row_ds} is not a group")
        gr = ds.create_group(name)
        if sparse:
            gr.create_dataset(
                "x", (len(data.data),), data=data.data, dtype=numpy.float32
            )
            gr.create_dataset(
                "i", (len(data.indices),), data=data.indices, dtype=numpy.uint32
            )
            gr.create_dataset("p", (len(data.indptr),), data=data.indptr)
        else:
            length = data.shape[0]
            total_len = data.shape[0] * data.shape[1]
            gr.create_dataset(
                "x", (total_len,), data=data.flatten("F"), dtype=numpy.float32
            )
            gr["length"] = [length]
        ds = self.get_datasource_metadata(row_ds)
        ds["links"][col_ds]["rows_as_columns"]["subgroups"][stub] = {
            "name": name,
            "label": label,
            "type": "sparse" if sparse else "dense",
        }
        self.set_datasource_metadata(ds)
        h5.close()

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
        from .server import create_app

        create_app(self, **kwargs)

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
        page = "page.html"
        template = join(tdir, page)
        page = open(template).read()
        # dummy popout page
        copyfile(join(tdir, "popout.html"), join(outdir, "popout.html"))
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
        # throttle the dataloading so don't get network errors
        conf["dataloading"] = {"split": 5, "threads": 2}
        save_json(join(outdir, "state.json"), conf)
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

    def add_image_set(self, datasource, setname, column, folder, type="png"):
        """Adds a set of images to a datasource. The images should be in a folder, with the same name as the column
        Args:
            datasource (str): The name of the datasource.
            column (str): The name of the column.
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
        if not ds.get("images"):
            ds["images"] = {}
        ds["images"][setname] = {
            "key_column": column,
            "type": type,
            "base_url": f"./images/{fname}/",
        }

        self.set_datasource_metadata(ds)

    def get_view(self, view):
        views = self.views
        return views.get(view)

    def set_view(self, name, view, make_default=False):
        views = self.views
        # update or add the view
        if view:
            views[name] = view
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
                "images": {},
            }
        # maybe warn if replacing existing regions
        # or add to existing regions
        md["regions"]["all_regions"] = all_regions
        self.set_datasource_metadata(md)

    def add_region_images(self, datasource, data):
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

    def add_viv_viewer(self, datasource, default_channels):
        md = self.get_datasource_metadata(datasource)
        reg = md.get("regions")
        if not reg:
            raise AttributeError(
                f"Adding viv viewer to {datasource}, which does not contain regions"
            )
        imdir = join(self.dir, "images", "avivator")
        if not exists(imdir):
            os.makedirs(imdir)

        reg["avivator"] = {
            "default_channels": default_channels,
            "base_url": "images/avivator/",
        }
        self.set_datasource_metadata(md)

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
                        print(
                            f"Cannot link viv image '{v['path']}' to '{datasource}'\n{repr(e)}"
                        )
                    region["viv_image"] = {
                        "file": os.path.basename(v["path"]),
                        "linked_file": True,
                        # "original_folder": os.path.dirname(v["path"])
                    }
        self.set_datasource_metadata(md)


def get_json(file):
    return json.loads(open(file).read())


def save_json(file, data):
    o = open(file, "w")
    try:
        o.write(json.dumps(data, indent=2, allow_nan=False))
    except Exception as e:
        print(
            f"Error saving json to '{file}': some data cleaning may be necessary... project likely to be in a bad state."
        )
        raise (e)
    o.close()


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
    col: dict, data: pandas.Series | pandas.DataFrame, group: h5py.Group, length: int
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
        if data.dtype == "category":
            data = data.cat.add_categories("ND")
            data = data.fillna("ND")
        else:
            # may need to double-check this...
            data = data.fillna("NaN")
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
                if v == "" or not isinstance(v, str):
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
        clean = data.apply(pandas.to_numeric, errors="coerce")
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
            col["quantiles"][str(q)] = [
                numpy.percentile(na, 100 * q),
                numpy.percentile(na, 100 * (1 - q)),
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


if __name__ == "__main__":
    path = os.getcwd()
    if len(sys.argv) > 1:
        path = sys.argv[1]
    if os.path.exists(os.path.join(path, "datasources.json")):
        MDVProject(path).serve()
