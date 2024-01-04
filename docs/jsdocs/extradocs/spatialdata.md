
## Spatial Data
Your data should contain an x and y column , which is the centroid for each cell/item and a column specifying which region each cell/item is from e.g. sample_id, roi, slide_number. In addition it needs to contain a column that will be used to color centroid plots e.g. annotations, cell types, cluster designations.

### Adding Regions
use the `set_region_data` method of the MDVProject. For example to set regions for the cell's datasource
 ```python
p.set_region_data("cells", "/path/to/file",
                region_field="sample_id",
                default_color="annotations",
                position_fields=["x", "y"],
                scale_unit="mm",
                scale=0.001)
```
The first argument is the name of the datasource and the second is a path to a tab delimited table or a dictionary describing the regions e.g.

|sample_id|width|height|
|---------|-----|------|
|sample_1_roi_1|2000|2000|
|sample_1_roi_2|1250|1500|

The first column should contain the name of the region (can be any header). The other columns are:-
* width - the width of the region 
* height -the height of the region
* x_offset, y_offset - usually these are 0 and in such cases, these columns can be missing. In some circumstances however the region may be part of a larger region and the hence the x, y co-ordinates will not be 0 and need to be supplied.

Alternatively instead of a path to a file, a dictionary of region name to data can be used e.g.
```json
{
    "sample_1":{
        "width":2000,
        "height":2000
    }
}
```
The other arguments are:-

* *region_field* - the column in your data specifying regions e.g. ROI, sample_id
* *default_color* - the default color of the  centroid plots, usually annotations or cluster designations
* *position_fields* - the fields describing the x, y co-ordinates of each cell/object
* *scale unit* - the unit you want the data described in
* *scale* - the actual value of x and y in terms of the scale unit


### Adding Images

Use `add_region_images' to add background images:-
```python
p.add_region_images("cells", "examples/data/spatial/images.tsv")
```
The second argument should be a path to a tab delimited table or a dictionary describing the images

|sample_id|path|name|
|---------|-----|------|
|sample_1_roi_1|/path/to/images/sample_1_roi_1/cellmask.png|cellmask|
|sample_1_roi_1|/path/to/images/sample_1_roi_1/he_stain.png|HE|
|sample_2_roi_1|/path/to/images/sample_2_roi_1/cellmask.png|cellmask|

The table should contain the region name (first column) then a path column which should contain either a url to an image or a local path to an image and name for the image. Multiple images for the same region are allowed but they must have different names.

x_offset, y_offset, width and height columns can also be included if the background image is larger or offset from the region. If these values are not present, they will be the same as those of the region. The user is also able to try and align the image with centroids.

If a local file path is given in 'path', the image will be copied to the project and the correct url for displaying it will be constructed.

## Adding a Viv view

Use `add_viv_viewer` of the `MDVProject'

```python
p.add_viv_viewer("cells",
    [
        {"name": "DNA1", "color": "#23133"},
        {"name": "collagen", "color": [223, 21, 123]}
    ]
)
```

The first argument is the datasource name, and the second the default channels. Contrast limits can be added , but usually these are calculated based on the image. Colors can be either hex or an RGB list. They can also be omitted and default colors will be assigned. These are only the initial values and channels can be added/removed and edited by the user. The user is also able to set new defaults


To add ome-tiffs use `add_viv_images`. The format is the same as for the normal regions images i.e. a path to a file or a dictionary and again can either be a url or a path to a local images

```python
p.add_viv_images("cells", {
     "COVID_OP_SAMPLE_ID_6_ROI_1": {
        "path": "/path/to/COVID_OP_SAMPLE_ID_6_ROI_1.ome-tiff"
    },
      "COVID_OP_SAMPLE_ID_6_ROI_2": {
        "path": "https://mysite/images/C_ID_6_2.ome.tiff"
      }
})
```


## Adding Interactions

First add a datasource, that describes each cell/cell interaction for each ROI or condition (group of ROIs). The datasource should contain the following:-
 * two columns describing the two cell/object types in the interaction these two columns should contain values that are
 present in a column in the parent dataset
 * two columns with the number of cells/objects of each type
 * a column designating the ROI or condition (also known as the pivot column)
 * any number of columns with  stats about the interaction

|Cell Type 1|Cell Type 2|cell 1 number|cell 2 number|sample_id|gr20|contacts|
|-----------|-----------|-------------|-------------|---------|----|--------|
|B cell|T cell |143|36|sample_6_roi_1|2.3|100|
|B cell|CD24   |143|24|sample_6_roi_1|0.8|27|


Then use the set_interactions methods to specify that the dataset contains interaction data and to link it to parent dataset e.g. cells
```python
p.set_interactions("my_interactions","cells",
                    pivot_column="sample_id",
                    parent_column="annotation",
                    is_single_region=True,
                    interaction_columns=["Cell Type 1","Cell Type 2"],
                    default_parameter="Cross PCF gr20",
                    node_size="cell 1 number",
                    add_view=True)
```

 * interaction_columns - a list of the two columns which specify the two cell types/items
 * node_size - the name of column in the interaction dataset that describes the number of cells/items in the interaction e.g. cell 1 number 
 * parent_column - the name of the column in the parent dataset that defines the  interacting items e.g. annotation, cell type
 * pivot_column - the column in both the interaction and parent dataset that specifies the extent of the interaction e.g.sample_id, condition, roi, disease
 * is_single_region - if True, the pivot column specifies a single roi i.e. it is not composed of several regions e.g.condition,disease state



