def add_centroid_plots(project, datasource, split_by, ignore_values=[], scale=0.5):
    col = project.get_column_metadata(datasource, split_by)
    md = project.get_datasource_metadata(datasource)
    val_to_region = {x: set() for x in col["values"] if x not in ignore_values}
    for sp_val, r_val in zip(
        project.get_column(datasource, split_by),
        project.get_column(datasource, md["regions"]["region_field"]),
    ):
        regions = val_to_region.get(sp_val)
        if regions is None:
            continue
        regions.add(r_val)

    for val, regions in val_to_region.items():
        regions = sorted(list(regions))
        if len(regions) > 16:
            hlf = int(len(regions) / 2)
            all_regions = [regions[:hlf], regions[hlf:]]
            region_names = [f"{val} centroid plots {x}" for x in ["1", "2"]]
        else:
            all_regions = [regions]
            region_names = [f"{val} centroid plots"]
        for c, regions in enumerate(all_regions):
            charts = []
            max_y = 0
            x = 5
            y = 5
            for region in regions:
                chart = project.get_centroid_plot(datasource, region, scale=scale)
                chart["position"] = [x, y]
                x += chart["size"][0] + 5
                max_y = max(max_y, chart["size"][1])
                if x > 1000:
                    x = 5
                    y += max_y + 5
                    max_y = 0
                charts.append(chart)
            project.set_view(
                region_names[c],
                {
                    "initialCharts": {
                        datasource: charts,
                    }
                },
            )


def create_interaction_matrix_view(
    project,
    interaction_datasource,
    exclude_groups=[],
    name="interaction view",
    cells_datasource="cells",
    interaction_metrics=["gr20"],
    square_size=20,
    images="pcf charts",
    selections=[
        {"column": "MH_FDR", "filter": [None, 0.05]},
        {"column": "MH_PC", "filter": [0, None]},
    ],
):
    md = project.get_datasource_metadata(interaction_datasource)
    i = md.get("interactions")
    if not i:
        raise ValueError(
            f"Datasource {interaction_datasource} does not have interaction data"
        )
    col = project.get_column_metadata(interaction_datasource, i["pivot_column"])
    groups = col["values"]
    charts = []
    y = 5
    x = 5
    for metric in interaction_metrics:
        x = 5
        chart = {}  # for pyright - not ideal.
        for group in groups:
            if group not in exclude_groups:
                chart = project.get_interaction_matrix(
                    interaction_datasource, group, metric, square_size=square_size
                )
                chart["position"] = [x, y]
                charts.append(chart)
                x += chart["size"][1] + 5
        y += chart["size"][0] + 5

    im_chart = project.get_image_plot(interaction_datasource, images)
    im_chart["position"] = [x, 5]
    im_chart["size"] = [300, 300]
    charts.append(im_chart)
    sel_dialog = project.get_selection_dialog(interaction_datasource, selections)
    sel_dialog["position"] = [x, 310]
    sel_dialog["size"] = [250, 250]
    charts.append(sel_dialog)

    project.set_view(
        name,
        {
            "initialCharts": {cells_datasource: [], interaction_datasource: charts},
            "dataSources": {
                cells_datasource: {"panelWidth": 20},
                interaction_datasource: {"panelWidth": 80},
            },
        },
    )


def create_roi_interaction_view(
    project,
    datasource,
    lchm=None,
    pcfs="pcf charts",
    cells_datasource="cells",
    extra_filters=[],
    table_metrics=["gr20"],
    name="roi interaction view",
):
    md = project.get_datasource_metadata(datasource)
    i = md.get("interactions")
    if not i:
        raise ValueError(f"Datasource {datasource} does not have interaction data")
    charts = []
    x = 5
    y = 5
    if lchm:
        chart = project.get_image_plot(datasource, lchm)
        chart["position"] = [x, y]
        chart["size"] = [300, 300]
        x += chart["size"][0] + 5
        charts.append(chart)
    if pcfs:
        chart = project.get_image_plot(datasource, pcfs)
        chart["position"] = [x, y]
        chart["size"] = [300, 300]
        charts.append(chart)
        x += chart["size"][0] + 5

    cols = i["interaction_columns"] + [i["pivot_column"]]
    filters = [{"column": x} for x in cols]
    filters += extra_filters
    sel = project.get_selection_dialog(datasource, filters)
    sel["position"] = [x, 5]
    sel["size"] = [250, 250]
    charts.append(sel)
    table = {
        "type": "table_chart",
        "param": cols + table_metrics,
        "column_widths": {
            i["interaction_columns"][0]: 200,
            i["interaction_columns"][1]: 200,
        },
        "size": [800, 200],
        "position": [5, 310],
    }
    charts.append(table)

    project.set_view(
        name,
        {
            "initialCharts": {cells_datasource: [], datasource: charts},
            "dataSources": {
                cells_datasource: {"panelWidth": 20},
                datasource: {"panelWidth": 80},
            },
        },
    )
