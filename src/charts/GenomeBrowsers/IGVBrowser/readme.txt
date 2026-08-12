### IGV Browser
This is a replacement for the built in genome browser which was based on an early version
of ivg.js


Because only the BaseTrack class is exported in igv, a custom track
will subclass BaseTrack and then dynamically create a proxy track of the desired
type and and use this proxy to modify how features are retrieved and rendered.

The default track, based on the associated datasource, no longer requires a bigbed or
tabix file, but the features are loaded in upon chart creation from the datastore.
This works with hundreds of thousands of rows but has not been tested on millions.

In full genome view (all chromosomes), filtering is applied to the downsampled features
passed to the getFeatures method of the track, hence the features shown are probably
not a true representation of the filter.

Any track type that IGV supports can also be added 