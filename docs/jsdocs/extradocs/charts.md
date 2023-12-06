# Chart Types

## GenomeBrowser

* **type** - genome_browser

* **param** - the columns the describe a genomic location e.g. chr, start and finish

* **tracks** An array of track objects. The very minimum it should contain a track with id `'_base_track'` that corresponds to the DataSource's data i.e. a gzipped bed file with a tabix index containing the columns chr, start. end and the index(id) of the row in the DataStore

* **feature_label** - the column to use in the DataStore to label the features (optional)

* **color_wig_tracks_by** - the column to color any tracks by. The track must have an id of `<column_field>`|`<value>`

* **view_margins**
    * type - either percentage (of the feature being viewed) or absolute (in bp)
    * value - either a percentage or bp value


# creating images from the chart

Charts should have a `getImage()` method that takes a function and a type, either 'svg' or 'png'. The callback should pass a serialized svg string in the case of svg or a canvas element in the case of a png.

```js
getImage(callback, type){
    if (type === "svg"){
        this.getSerializedSVG().then((svg)=>{
            callback(svg);
        });
    }
    else if type === "png"{
        callback(this.canvas);
    }
}

```
The chart will then download the image at the user's request


For charts which extend `SVGChart`, this function has already been implemented and the svg object will be saved as an svg or png image. However, if there are other elements that are not part of the svg but need to be added to the image, you can implement the methods addToImage and removeFromImage.

```js

addToImage(){
    this.svg.append("text").attr("text", "Added Text");
}

removeFromImage(){
    this.svg.select("text").remove();
}

```

# Downloading data
The chart should have a `getChartData` method which returns a text blob e.g. tab delimited table