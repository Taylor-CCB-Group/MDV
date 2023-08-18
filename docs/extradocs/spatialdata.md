
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

## Adding an Avivator view

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



