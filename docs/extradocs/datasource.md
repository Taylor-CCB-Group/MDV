
The DataStore is created with a size (number of rows) and a config 

```
const ds =  new DataStore(1000,{....})
```
In theory The config is optional, columns and their data can be added later, but practically a config is required and  can contain the following  

## Config

### name

The name of the DataStore - should be human readable

### columns

An array of column objects 

The DataStore can contain a number of columns. Each column can have the following

* **name** The column's label that is exposed to the user
* **field** The column's id - should be unique

* **datatype**

  * text - string can't have more than 256 values 
  * integer - an integer (represented by float32)
  * double - a floating point number (represented by float32)
  * unique - any text
  * multitext - A field that can have more than one value eg 'red', 'red,blue'
  * int32 - for larger integers e.g. genomic co-ordinates (represented by an int32)

* **values** - For text and multitext colummns, this will be an array of the possible values, the raw data will contain the index(es) of the value(s) in the array. Text columns can't have more than 256 values and multitext no more than 65536

* **colors** - An array of colors (in hex format). Only used for text/multitext and
integer/double columns. In the former the colors array is mapped to the values array. In the latter, the colors are interpolated between max and min values (or quantiles if specified). Default colors will be used if not supplied

* **stringLength** - for unique columns, this is the length of the longest value and for multitext the maximum number of values associated with a data item.

* **colorLogScale** - true of false - for number columns, whether the color should be a log scale 

* **editable** - specifies a column can be edited.

* **is_url** - the column contains url links (unique and multitext columns)

* **minMax** - for integer and double columns, the minimum/maximum values in an array of length 2

* **quantiles** - an object with 0.05,0.01 and 0.001 as keys, with each entry and array of the x and 1-x quantiles
```
    {
        "quantiles":{
            "0.05":[x,y]
            "0.01":[x,y]
            "0.001":[x,y]
        }   
    }
```


* **data** -  the actual data of the column either an array or SharedArrayBuffer. Usually the column data is loaded in later on demand see column data



If the column's data is supplied as an array, then values and stringLength need not be supplied and can be calculated. However if the column's data is supplied as a SharedArray Buffer these values are required. 

minMax and quantiles will be calculated if not supplied but this may be slow for large datasets and hence should given in the config

### columnGroups

A list of logically grouped columns, in the following format. Columns should be referenced by their field values. A column can belong to more than one group
```
    "columnGroups":[
        {
            "columns":[
                "area",
                "eccentricity",
                "perimeter"
            ],
            "name":"cell Stats"
        },
         {
            "columns":[
                "CD34",
                "CD38",
                "GranzymeB"
            ],
            "name":"markers"
        }
    ]
```


### links
A dictionary where the keys are the dataStores that this dataStore should link with and the entries are objects describing the links. See the section for linking to other DataStores

```
    "links":{
        "ds1":{.....},
        "ds2":{.....}
    }
```

### images

Represent thumbnails for each data item. It should be an object containing 'image sets', with the key being the name of the set. Each set should have:-

* base_url - the base url of the image -  the images key and type will be appended to this url
* key_column - the column in the dataStore that contains the keys for each image
* type - the type of image

```
    {
        "images":{
            "set1":{
                "base_url":"/my_thumbnails/set1/"
                "key_column":"set1_image_key",
                "type":"png"
            }
        }
    }
```
Charts that use the images in set1 would then associate the row whose set1 _image_key was 26373 with the following url:-  
    
    /my_thumbnails/set1/26373.png


### large_images

These have exactly the same format as above, but represent much larger images, which are intended only to be displayed one or two at a time.

### background_images

Images which represent a category and have x, y co-ordinates

```json

"rois":{
    "coordinates":["x","y"],
    "column":"sample_id",
    "color_by":"cell_type",
    "images":{
        "sample_1":[
           {
                "roi":{
                    "min_x":0,
                    "mix_y":0,
                    "max_x":200,
                    "max_y":200
                },
                "images":[
                    {
                        "url":"/...",
                        "name":"my_image",
                        "position":[0,0],
                        "height":200,
                        "width":200
                    }
                ]
           }
        ]
    }

}

```


### offsets

Sometimes, the  x and y co-ordinates of certain groups within the data need to be changed to align with each other. 
This is enabled by adding an offsets parameter to the config which is comprised of :-

* param -  the columns that can be offset (e.g. "x" and "y" )
* groups - which column specifies the groups that can be offset
* background_filter - if the data consists of multiple spatial data, a background_filter referring to the column which represents each spatial entity e.g."ROI"
* values -  a nested dictionary for each filter and group describing its offset and rotation. This is optional, and can be omitted if there are currently no offsets or rotations. Each value should contain the following:-
    * rotation - rotation +/- in degrees
    * offset - array with x and y offsets
    * rotation_center - the relative rotation point - the actual point is +/- the offsets. if omitted, a default value based on the center if of all the points is used



An example for x and y columns, where each panel can be offset in each ROI and the offset/rotation for panel1 in ROI_1 is specified is shown below:-

```
    "offsets":{
        "param":["x","y"]
        "groups":"panels",
        "background_filter":"ROI",
        "values":{
            "ROI_1":{
                "panel1":{
                    "offset":[10,20],
                    "rotation":12
                }
            }
        }
    }
```
If no background filter is specified then a single entry "all" should be used in values:-
```
    "offsets":{
        "param":["x","y"]
        "groups":"panels",
        "values":{
            "all":{
                "panel1":{
                    "offset":[10,20],
                    "rotation":12
                }
            }
        }
    }
```


To set the values use the `setColumnOffset` method 

```
//rotate by 45 deg
let  data = {group:"panel1",rotation:45,filter:"ROI_1"};
ds.setColumnOffset(data);
//translate by 30,30
data = {group:"panel1",offsets:[30,30],filter:"ROI_1"};
//passing true as a second parameter will update all listeners
ds.setColumnOffset(data,true);
```
If no background filter has been specified, you should omit the filter value or set it to "all"

To reset a group's offsets to the original values use `resetColumnOffsets`
```
ds.resetColumnOffsets("ROI_1","panel1",true)
```
The first parameter should be null or "all" if there is no background filter


### Genome Browser

Allows a genome browser, to be added as a chart. The default feature track will show each item in the DataStore.Has the following parameters:-

* **default_parameters** (optional) an object containing parameters for the Genome Browser

* **default_track** The track associated with the data. Needs to be  gzipped bed file, indexed using tabix and have 4 columns:-  chromosome, start, end and the id of the row in the DataStore
    * **label** The name of this track
    * **url** The url of track

* **default_track_parameters** (optional) - An object containing parameters for the feature track

* **location_fields** An array containing the id of the corresponding columns in the DataStore for chr start and end

* **default_tracks** (optional) A list of track configs for any tracks which will be displayed when the genome browser is added (A ruler track is always added and need not be specified here)

* **atac_bam_track** (optional) This track will cluster reads based on the barcode tag in the bam file. The DataSource needs linked to another DataStore that has these barcodes for the atac_bam_track to display e.g.
     ```
     "links":{
        "cells":{
            "index":"cell_id",
            "access_data":true
        }
     }
     ```
     * **url** the url of the bam track
     * **cluster_read** the field in the linked dataset with which to cluster the reads


Example:
```js
    


```

## Column data

Internally every column is a typed array, which makes manipulation much quicker and allows webworkers to work on them without holding up the main thread

Column data can be added by the following

* Adding a data parameter to the column objects in the columns parameter of the config passed to the DataStores's constructor

* Using the DataStore's `addColumn` method passing the data parameter as the second parameter
    ```
        ds.addColumn({datatype:"text","field":"myvals",name:"My Vals"},data)

    ```

* If a column has already been added,  use  `setColumn` method with the column's field and data
    ```
        ds.setColumnData("myvals",data)
    ```

In each case, data can either be an array of values or a SharedArrayBuffer containing the raw format of the column. Using an array, although more convenient, will be slow, as it has to be converted to the native format


### datatype - text
JavaScript Array -  strings - there should be no more than 256 unique values.

SharedArrayBuffer - A uint8, each position in the underlying array represents the index in the column's value array 

```
    ["blue","green","green","yellow","blue","green"]
    //would be converted to 
    values:["green","blue",yellow]
    data:[1,0,0,2,1,0] //(Uint8Array)
```


### datatype- mulitext

uint16 each data item is stringLength portion of the array with each position being the index of the value's array. 65535 represents no value

```
    [ "A,B,C", "B,A", "A,B", "D,E", "E,C,D" ]
    //would be converted to
    values:["A","B","C","D","E"]
    stringLength:3
    data:[0,1,2, 1,0,65535, 0,1,65535, 3,4,65535, 0,2,3] //(Uint16Array)
```


### datatype - unique
A unique column represents text that can contain more than 256 values and is represented by a Uint8Array. Each value is encoded by stringLength digits, padded by 0 values .  
```
    ["ZX1212","X21","F232","D1"]
    //would be converted to
    data:[90,88,49,50,49,50, 88,50,49,0,0,0, 70,50,51,50,0,0, 68,49,0,0,0,0] //(Uint8Array)
    stringLength:6
```


### datatype integer/double
These are treated the same and are represented by Float32Array
```
[1.2,3,4.5.7]
//would be converted to
data:[1.2,3,4.5.7] //(Float32Array)
```

### datatype integer/double
This is represented by a Int32Array and is better far larger integers. e.g.  genomic locations. However, it can't hold an NaN which represents missing data
```
[1.2,3,4.5.7]
//would be converted to
data:[1.2,3,4.5.7] //(Float32Array)
```


## Linking to other DataStores

The links parameter of a DataStore's config  allows different types of interaction between DataStores



### rows_as_columns
This specifies that the DataStore can contain data linked to the rows another DataStore. 
For example, in single cell data, the cell DataStore would be linked to the gene DataStore. Data such as gene expression per cell could then be added to the cell DataStore as columns (on demand - not all at once). It should have the following parameters:-

* **name_column** - the column in the linked dataset which identifies the row e.g. gene name (must be unique)
* **name** - the human readable name of the dataset e.g Gene Scores
* **subgroups** - this specifies the  different groups of data available. It consists of a dictionary, with the identifying stub of the subgroup pointing to the name, type and label (what the user will see) of the subgroup.

        "subgroups": {
            "gs": {
                "name": "gene_scores",
                "type": "sparse",
                "label": "Gene Score"
            },
            "igs": {
                "name": "imp_gene_scores",
                "label": "Imputed Gene Score"
            }
        }

The number of columns could potentially be very large, depending on the number of rows in the linked dataset, thus unlike other columns these are not explicitly specified in the DataStore's config. The column's field will have the following format:-

    <stub>|<name>|<index>

* **stub** is the identifier of the subgroup
* **name** is the name shown to the user
* **index** the index of the associated row in the linked dataset and should also be used retrieve the data for this column. For example, a config for a histogram of imputed gene score for gene EMB

        {   
            "type":"bar_chart",
            "param":"igs|EMB(igs)|342"
        }


The columns's data will be loaded in the normal way by the specified dataloader, but the column object passed to the dataloader will also have the following fields:-

* **sgindex** the unique identifier 
* **sgtype** the type of data e.g. sparse
* **subgroup** the stub of subgroup


example
```

   [
        {
            "name":"cells",
            "columns":[],
            "links":{
                "genes":{
                    "rows_as_columns":{
                        "name_column":"gene_name",
                        "name":"Gene Scores",
                        "subgroups": {
                            "gs": {
                                "name": "gene_scores",
                                "type": "sparse",
                                "label": "Gene Score"
                            },
                            "igs": {
                                "name": "imp_gene_scores",
                                "label": "Imputed Gene Score"
                            }
                        }
                    }
                }
            }
        },
        {
            "name":"genes",
            "columns":[
                {
                    "name":"Gene name",
                    "field":"name",
                    "datatype":"unique"
                }
            ]
        }
   ]
``` 


### Linking to Columns

You can link columns from two DataStores in order to minimize the duplication of data. For example if a cell dataset contained cells from different samples and a sample dataset contained information about these samples, the cell DataStore should have access to this metadata, without the duplication of data. To achieve this, include a link to the other DatStore, specifying which columns are needed from that DataStore and the index, which is the column that the two DataStores share and must be of datatype text. You must also describe the  linked columns in the 'columns' parameter with the fields matching the linked DataStore.  The colors/values can be omitted as these will be taken from the linked DataStore.  However column names(what the user sees) can differ

```
 {
        "name":"cells",
        "size":101737,
        "columns":[
            {
                "name":"sample_id",
                "datatype2:"text"
            },
            //the following columns will get their data from the samples DataStore
            {
                "name":"sample condition",
                "datatype":"text",
                "field":"condition"
            },
            {
                "name":"sample sex",
                "datatype":"text",
                "field":"sex"
            }
        ],
        "links":{
            "samples":{
                "columns":["condition","sex"],
                "index":"sample_id"
            }
        }
    },
    {
        "name":"samples",
        "size":11,
        "columns": [
            {
                "name":"sample_id",
                "datatype":"text",
                "values":["1","2","3","4"]
            },
            {
                "name":"condition",
                "datatype":"text",
                "values":["sick","healthy"],
                "colors":["red","green"]
            },
            {
                "name":"sex",
                "datatype":"text",
                "values":["male","female"],
                "colors":["blue","pink"]
            }
        ]
    }

```


### Accessing columns


For a dataStore to access columns of another it needs to specify this in the links and index column

```
links:{
      "cells":{
        "index":"cell_id",
        "access_data":true
      }

}
```

### synchronizing column colors

To link column colors , specify the column from the linked DataStore you wish the colors to match 'link_to" and the column you want to match this column 'col'. You can link one column to many

```
{
    name:"cell interactions",
    links:{
        "cells":{
            "sync_column_colors":[
                {
                    link_to:"Cell Type",
                    col:"Cell Type 1"

                },
                {
                    link_to:"Cell Type",
                    col:"Cell Type 2"
                }
            ]
        }

    }
}
```

## Listeners

Listeners can be added/removed from the DataStore with the `addListener` and `removeListener` methods. The callback passed to the addListener will be invoked with the type of event and any data relevant to the event.

* data_changed
    *  columns - a list of column fields where the data had changed
    *  hasFiltered - true if the changing the data has caused the refiltering of the DataStore 

* filtered - The DataStore has been filtered, no data is passed to the listener unless all filters have been removed, in which case 'all_removed' is passed.

* column_removed - called when a column is removed passing the column's field (id)

* data_highlighted 
    * indexes - a list of indexes of the highlighted items(rows)
    * source - the object doing the highlighting




## Filtering

This is achieved by getting a dimension from the datastore, then calling the filter() on the dimension, specifying the filter method, columns and arguments

### category dimension

This has two filtering methods filterCategories works on a single column whereas filterMultipleCategories works on multiple columns

    

```
    const dim = datastore.getDimension("category_dimension");
    //filter items with red in color column
    dim.filter("filterCategories",["color"],"blue");
    dim.removeFilter();
    //filter items that are red or blue
    dim.filter("filterCategories",["color"],["red","blue"]);
    dim.removeFilter();
    //filter items that are red and have a shape value of square
    dim.filter("filterMultipleCategories",["color","shape"],["red","square"])
``` 


### range dimension

```

    const dim = datastore.getDimension("range_dimension");
    //filter on between 5 and 10
    dim.filter("filterRange",["x"],{min:5,min:10});
    dim.removeFilter();
    //filter on x and y (a square)
    dim.filter("filterSquare",["x","y"],{
        range1:[5,10],
        range2:[5,10]
    });
    dim.removeFilter();
    //filter based on a polygon- just supply the points
    dim.filter("filterPoly",["x","y"],[[0,0],[2,0],[1,2]])
```

All dimensions can filter on an arbitrary set of items using the item's index

```

    const indexes = new Set([4,5,10,12]);
    const dim = datastore.getDimension("range_dimension");
    dim.filter("filterOnIndex",null,indexes)
``` 

When a filter is called on a dimension - the datastore informs all listeners and hence any updates are performed, which is likely to be expensive. Therefore, if you want to perform multiple filters, pass false as the fourth argument and then call triggerFilter on the datastore

```
    rangeDim.filter("filterRange",["x"],{min:5,min:10},false);
    catDim.filter("filterCategories",["color"],"blue",false);
    datastore.triggerFilter();
```

    