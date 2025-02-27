.. |expand| raw:: html

   <i class='fas fa-plus-circle'></i>

.. |zegami| raw:: html

   <img src='http://martindev1.molbiol.ox.ac.uk/static/img/icons/zegami.png'/>

.. |clone| raw:: html

   <i class='fas fa-clone'></i>

.. |download| raw:: html

   <i class='fas fa-download'></i>

.. |history_icon| raw:: html

   <i class='fas fa-history'></i>

.. |question| raw:: html

   <i class='fas fa-question'></i>

.. |tss_icon| raw:: html

   <i class='fas fa-exchange-alt'></i>

.. |stream| raw:: html

   <i class='fas fa-stream'></i>

.. |sort_up| raw:: html

   <i class='fas fa-sort-alpha-up'></i>

.. |color_palette| raw:: html

   <i class='fas fa-palette'></i>

.. |save| raw:: html

   <i class='fas fa-save'></i>

.. |camera| raw:: html

   <i class='fas fa-camera-retro'></i>

.. |cog| raw:: html

   <i class='fas fa-cog'></i>

.. |filter| raw:: html

   <i class='fas fa-filter'></i>

.. |table| raw:: html

   <i class='fas fa-table'></i>

.. |tags| raw:: html

   <i class='fas fa-tags'></i>

.. |images| raw:: html

   <i class='far fa-images'></i>

.. |call_peaks| raw:: html

   <i class='fas fa-signature'></i>

.. |share_icon| raw:: html

   <i class='fas fa-share'></i>

.. |cluster| raw:: html

   <i class='fab fa-cloudsmith'></i>



.. _multi-locus-view:

Multi Dimensional Viewer
#########################

Summary
===================

Multi Dimensional Viewer (MDV) enables the user to visually inspect and query the underlying data from omics experiments e.g. scRNA-seq, scATAC-seq, Spatial trancriptomics / Proteomics
Initially, the user uploads a table of data (.csv/.tsv), anndata hd5 files (.h5 or .h5ad) or images (.ome.tiff) and associated data (`Uploading a File`_), which can then intuitively filtered, sorted and viewed to drill down to regions of interest. 

Each MDV project comprises a series of Views which correspond to table(s) of data. There can be multiple tables in one View and these can be joined. For example, in a single cell data set there is a 'cells' table and 'genes' table.  
This allows you to query your favourite gens and looks at which cells have ahigh expression value (for example).

Creating a Project
===================

Register and Login to MDV and you will see the MDV Project Viewer. This shows all your projects as folders and can 

To create a project click on 'My Projects' (1) on the top navigation bar and then the  Multi Dimensional Viewer panel (2). This will take you to a page where you have to fill in the name and description of the project (3).


Uploading Data
----------------

Input File(s) 
++++++++++++++++++++++++
MDV supports (.csv/.tsv), anndata hd5 files (.h5 or .h5ad) or images (.ome.tiff). 


Uploading a File
+++++++++++++++++++++

Go to your newly created project and drag and drop the file in to the Upload dialog box. You should see a preview of the file contents to verify the data
and then you can press "OK" to upload it. We support files up to 2GB via the user interface. Larger than that we recommend you upload via our API.

Once you are satisfied that the file contents appear to be correct in the preview screen, press the upload button. There will be a delay as the data is uploaded and processed and if there are no problems 
you will be presented with an initial view in the project showing the contents of raw data in the file.

For tabular data such as csv or csv this will be a single table.
For anndata this will normally we shown as a split screen with a "cell" table and a "gene" table.
For an image you should see the image viewer. 



Saving A Project
----------------
When tracks, graphs or columns are added and you have edit permissions (see `Permissions`_), they are permanantly added to the project. 
However, when you alter a graph or tracks's settings or reize/move them, these are not saved. Similary any changes you make to the table (column resizing/ordering, sorting etc.)
will not automatically be saved. In order to save these changes you need to save the layout  |save| > Save Layout.

If you wish to make changes to public project then you will have to clone the project |save| > Save As. 
The project will be copied into your name and you can them make any changes you wish.


.. _MDV-view-data:


Adding Graphs/Charts
=====================

Charts help you get a picture of the data as a whole and also help you filter the data. 
By selecting regions (dragging with the mouse) on scatter plots and histograms or clicking on sections in pie charts, row charts and box plots, the data can be intutitively filtered. 
With each filtering step, all charts will update (as well as the table and browser) to reflect the filtered data. 
Filters on individual charts can be removed by clicking the reset button which appears on the chart's title bar when a filter is applied or filters on all charts can be removed with the 'Reset All' button. 

Charts can be moved by dragging on the title bar, and resized by dragging on the resize icon, which appears in the bottom right hand corner of chart when you hover over a chart. 

Initially the only chart visible will be the raw data view that was uploaded that was uploaded so you need to add other charts to get a better insight into your data (see below)

.. _MDV-adding-a-chart:

Adding a Chart
----------------

Clicking on the 'Add Chart' button will show a dialog where you have to select the type of chart, the fields to use in the chart and its name. Once created you can change the chart's settings (|cog| icon), which differ according to the chart's type and with some charts color it (|color_palette| icon). Charts can moved by dragging them via the title bar and resized by the resize icon which appears in bottom left hand corner when the mouse is over the chart. The chart can be removed by clicking the trash icon, which appears when you hover over the graph's title. Once charts have been added and the appropriate settings/colors added, they can be saved using the |save| icon above the table. The following chart types are available

.. _MDV-scatter-plot:

Scatter Plot
+++++++++++++++++

A standard scatter plot requiring numerical fields for the x and y axis. Once created, the points can be coloured (|color_palette| on the title bar). Also by opening up the settings dialog (|cog| icon) you can alter the point size (3). By default the graph will show all the points, but you can zoom in and out using the mouse wheel and pan by pressing shift and dragging with the mouse. However if you want the default view to be a particular region, you can set this using the inputs in the Main Region section (4)and pressing show. The x and/or y axis can also be set to a log scale (5). After zooming and panning, the Centre Plot button (6) will restore the plot to show all the points or the region specified in (4). Normal mouse dragging (without shift) will cause pruduce a brush that filters the data, once created the brus can be dragged to different regions of the plot.  

Histogram
+++++++++++++++

Shows the distribution of a numerical field in the data. The x range is automatically set to include the largest and smallest values. However, this will often lead to the chart looking skewed due to low numbers of extreme outliers. Therefore, you can use the cog icon |cog| (1) to open up the settings dialog, where an upper and/or lower limit can be set (3). Values higher or lower than these limits will be lumped together in the outermost bin (4). The y axis can also be capped (5) in order to get a better handle on bins conatining fewer counts. The number of bins can be adjusted using the appropriate spinner(6). Each bar can be coloured by categorical data use the |color_palette| icon (2). 


Pie Chart
+++++++++++++
Shows categorical data. By default the maximium number of categories shown are the 8 largest ones, any reamining categories are lumped into 'Others'. This can be changed by opening up the settings dialog (|cog|). Clicking on a segment (category) will select that category and clicking on further segements will add these to the filter. To filter again with a single category, use the reset button.

Row Chart
+++++++++++++
A chart showing categories on the y axis and usually the number of records belonging to this category on the y axis. You can also choose a numerical field for the x axis, in which case the values of this field will be summed for each category. However a boxplot is usually more informative for this kind of information as the average and quartile ranges of the values are shown instead of the sum. As with the pie chart, the maximium number of categories shown are the 8 the largest ones, but this can be changed by opening up the settings dialog (|cog|)

BoxPlot
++++++++++
A chart showing categeories and average/quartile ranges of the values of another field for that category. Box plots work best for fields that contain only a small number of categories. They are scaled to include all the datapoints, so if there are extreme outliers, the boxes will appear squashed.

Bar Chart
++++++++++
A Bar Chart showing the column average of any number of supplied fields. Because fields may differ in scale,to ensure differing values can fit on the same scale, the average is scaled between the median +/- 1.5*IQR (the same as the whiskers on a boxplot). 
The graph changes to only include those datapints in the current filter. No Selection is possible with this chart, as it would make no sense to filter on a column 




The Genome Browser
====================

The browser shows the genomic location of the currently selected table row (or image). The distance either side of the region to also show can be controlled using the margin spinner (1)  above the browser.


.. _MDV-adding-tracks:

Adding Tracks
-------------

Initially only two tracks will be displayed, the genomic locations you uploaded and if you didn't select 'other' for the genome, a track is displaying the genes. Other tracks can be added with the  'Add Tracks' button (2), which shows a dialog where you need to enter the url of a publically accessible track. The hosting server of the track should allow Cross Origin Resource Sharing (CORS). The type of track will try and be ascertained based on the url, although you can manually overide this by clicking on one of the radio buttons

Tracks that can be added are:-

* bed(tabix) - A bed file that has been gzipped and indexed with tabix
* BigBed 
* BigWig
* Bam - A bam index is also required
* UCSC session - either cut and paste the url from the UCSC browser or use a session url. The latter will be more stable as the former uses a temporary id, which is only valid for a short period.


Altering Track Appearance
----------------------------

Clicking on the track label in the legend (3) will open a dialog for that track. The contents of the dialog will vary
according to the type of track. The track height can be altered from this dialog

Zooming/Panning
------------------------

There are five ways you can navigate using the browser:-

   * You can zoom in and out using the mouse wheel and scroll left and right by deagging the mouse
   * Use shift and drag to highlight and zoom into a region on the browser
   * Use the magnifying glass icons (4), the zoom amount can be controlled by the adjacent spinner (5)
   * Type the location co-ordinates (chr:<start>-<stop>) in the location text box (6)
   * Click on a row or image in the right hand table to go to that feature. The margin spinner (1) shows how many bp  either side of the feature will be displayed. 

Feature track
-------------------
This shows the uploaded regions(features) displayed in the right hand table. Clicking on the settings |cog| icon (7) will bring up a dialog where the following can be adjsuted:-

    * Set the field you wish to a label the features with
    * Set the field to color the feature by
    * Set the field with which to position the feature on the y axis. By default the feature layout depends on the layout (Collapsed, expanded or squished)
      but can be a numeric field 
    * Choose the margins (distance either side of the feature) that will be displayed when you click on an image or a row in thr table


Saving the Browser Layout
-----------------------------------
Use the disk icon |save| above the the table to save all settings including the current layout of the browser (tracks and track settings)

Capturing An Image
---------------------------
Use the |camera| icon (9) to download an image of the current browser view. The image format (png, pdf or svg) can altered using the adjacent dropdown (10)



The Table/Images
===================
The default table behaves as a typical spread sheet, you can alter the column width by dragging the header's left and right borders and move columns by dragging the column's header. Clicking on the header will sort by that column. Clicking on that row will select it and update the browser.

Table Mode
----------------


If your project contains images (see `Adding Images`_) then then you change how the table is displayed using the table icon (|table|). Three choices are table (1), images (2) and table with thumbnials (3).In image mode, the genomic location can be selected by clicking on the image and using the arrow keys to select the next/previous image. In this mode, the data can be sorted and filtered using the icons (|filter| |sort_up|) in the menu above the table.
Also in image mode you can alter image size using the slider in the table menu and also color the border around the image by a field (|color_palette|). This opens up a dialog where you can choose the field and the color scale to use


Filtering Data
----------------------
It is often more intuitive to filter using graphs (see `Adding Graphs/Charts`_ ), however data can be filtered by clicking on filter icon |filter| in each column header. To filter on multiple columns or when the table is only showing images,press the filter icon on the top table menu. This will bring up a dialog showing filtering options for all fields in the data. Whenever any filters are added or changed, any charts will update accordingly,but the filters are not added to the charts or existing filters on the charts updated as they are completely independant.

Sorting Data
-----------------
The data  can be sorted on columns by clicking the column header (shift click to sort on multiple columns). The data can be also be filtered by clcking the sort icon |sort_up| in the table menu. In the sort dialog,the columns to be sorted on are added usng the plus icon and then either Ascending (ASC) or descending (DESC) can be chosen . The sort order can be changed by dragging the labels or columns removed from the sort by clicking on the trash icon

.. _MDV-tagging-locations:

Tagging Locations
====================

Sometimes it may be useful to catgeorise or tag the genomic location based on a trend theat that you have discovered. This can be done by opening up the tagging dialog with tag icon |tags| (1) in the menu above the table. Initially only the none category is present. To add other ones type a name in the text box (2) and press the add button (3) . The category will then be added to the list at the top of the dialog. By selecting the radio button next to it, then clicking on an image or a cell in the tagging column in the table will tag that genomic location. Multiple locations can be tagged by clicking and image/cell and the shift clicking another one and all the images/rows in between will be tagged. The 'Tag All' will tag all the currently filtered locations with the currently selected catogry. Another way to tag is to use the arrow keys to go to the next previous image/row and then press the shortcut key shown in brackets next to the category to tag the currently selected items with that category. The category color can be changed by clicking on the appropriate color chooser (7). The category can be removed (which will remove all tags of this category from the data using the trash icon next to the category (8)

*N.B.* To permanantly save the tags press the Save button (5) which will commit the changes to the database.

.. _MDV-adding-images:

Adding Images
=====================



Creating Subsets
=================

Clicking on |clone| brings up a dialog, which allows you to create a subset of the currently selected locations

You can choose to create the subset either from the currently filtered locations (1) or from a random subset (2) with a specified number (3). 
After filling in a name and description (4) and (5), press 'Create' and the subset will be created. 
Once this has happened, you will get a link to the subset. You can create a susbset of any project you have viewing rights too, including public projects, and you will 
be the owner of the subset. All graphs/tracks/columns are copied, although the graphs may look different as there will be fewer locations in the new project.
 

Exporting Data
===============

Click on the download icon |download| to download the currently filtered locatons. The data is just downloaded as a text file, 
although you get a choice to download in either tsv or csv format. 
Only the currently displayed columns will be downloaded, so expand any column groups by clicking the plus icon |expand| if you want these in your file. 



Permissions
============
There are two types of permission for a project, view and edit. 

If you have a view permission for a project (anyone has view permission for a public project), you can open the project and add tracks and charts, as well as edit existing charts and tracks. However, you cannot save any updates or run any jobs such as finding TSS's or creating images.
If you want to do this, you will have to copy the project (you need to be logged in) - click the  disk icon (|save|)  and then select 'save as'. 
This will clone the project in your name and then you can make any changes you wish.

If you have editing rights to a project you can make any changes you want , run any jobs and save the layout |save| . You automatically have editing rights to a project you own it or if you are an adminisrator.  You can also be assigened editing rights to a project (see below)

The |share_icon| icon on the menu allows you to share the project with another user and assign them view or editing rights 

Sharing a Project
-----------------

Click on share icon (|share_icon|) above the table and select 'share project' 

Start typing the name of user you want to share the project with in the text box (1) and then select the name from the drop down.
When you press 'Add' (2), the project will be shared and the name of the user will be added to the dialog (3). You can change the editing rights of the user to view or edit using the dropdown (4) or stop sharing by pressing the trash icon (5)

Making a Project Public
------------------------

Click on share icon (|share_icon|) above the table and select 'share project', you will be prompted to see whether you really want the 
project to become public. If you click OK then then anyone (including non users) will be able to view the project. You can share the project by sharing the link in the browser's address bar.

Submitting an Issue
====================

An Issue or question, can be asked within MDV (if you are logged in) using the help link  in the top navigation bar |question|  > 'Send Question'.
You can also submit an issue to he `GitHub page <https://github.com/Hughes-Genome-Group/MDV/issues>`_


Frequently Asked Questions
============================

Can MDV be viewed on a mobile/small screen device?
---------------------------------------------------

No. The whole idea is to see how each component ie. the graphs,tracks and images change as you filter the data, which would not be possible,
if only one component was displayed at once. 
All panels and individual tracks/graphs/table columns/images can be resized, to get the exact layout that the user requires, rather than relying on adapative screen size techniques which limit viewing to a single compnent/panel on small screen sizes

