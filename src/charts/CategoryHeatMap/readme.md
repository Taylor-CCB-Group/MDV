### CategoryHeatMap

This chart displays the categorical values from a column (text/text16 or multitext) on each axis
Each  cells show the number of rows which contain both values. When there are more than 2000 cells , a warning message is displayed  instead of the heatmap
Apart from filtering, you can restrict the categorical values in the chart using local settings. Although background filters could have been used, this was
confusing for multitext , because categories not selected would also be shown if 
they co-occurred with the selected categories
Clicking on a cell will filter for the two categories represented by that cell