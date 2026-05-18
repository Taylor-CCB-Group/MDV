Background Filters apply at the chart level and are permanent filters, i.e. they are independent of a chart's local filter and global filters.

They were used initially only in scatter plots to show specific regions in spatial data, but there was no way for a user to add/modify them.

Before, a dimension had a separate array to store the background filter, but now it is integrated into the local filter. This way no extra data has to be passed to the worker. Local filters can now have a fourth state (local filter + background filter).
0 - no filter
1 - local filter
2 - background filter
3 - local filter + background filter


They are applied with setBackgroundFilter and clearBackgroundFilter of Dimension. At the moment, only one column can be used and it has to be text,text16 or multitext. 

The filterPredicate function of the base dimension now takes background filters into consideration.

Individual charts/dimensions have to decide how they are going to handle the background filter when aggregating/displaying data. To date only dotplot and the classic scatterplot implement background filters.