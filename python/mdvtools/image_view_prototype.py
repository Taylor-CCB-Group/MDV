import json

# Define the default JSON string (as a multi-line string).
default_image_name = "Unvaccinated-Lesioned_SAMPLE_45_ROI_1"
default_ds = "cells"

json_string = '''
{
    "initialCharts": {
        "cells": [
            {
                "course_radius": 1,
                "radius": 10,
                "opacity": 1,
                "color_by": "leiden",
                "color_legend": {
                    "display": false
                },
                "tooltip": {
                    "show": false
                },
                "category_filters": [],
                "zoom_on_filter": false,
                "point_shape": "circle",
                "contour_fill": false,
                "contour_bandwidth": 0.1,
                "contour_intensity": 1,
                "contour_opacity": 0.5,
                "title": "Unvaccinated-Lesioned_SAMPLE_45_ROI_1",
                "legend": "",
                "type": "VivMdvRegionReact",
                "param": [
                    "x",
                    "y",
                    "leiden"
                ],
                "region": "Unvaccinated-Lesioned_SAMPLE_45_ROI_1",
                "background_filter": {
                    "column": "sample_id",
                    "category": "Unvaccinated-Lesioned_SAMPLE_45_ROI_1"
                },
                "roi": {
                    "min_x": 0,
                    "min_y": 0,
                    "max_y": 9323,
                    "max_x": 5661
                },
                "json": "json/Unvaccinated-Lesioned_SAMPLE_45_ROI_1_whole-cell_DAPI_AVGMARKER_200_75.tif.s1.json",
                "viv": {
                    "channels": [
                        {
                            "name": "DAPI"
                        }
                    ],
                    "channelsStore": {
                        "channelsVisible": [
                            true,
                            true,
                            true,
                            true
                        ],
                        "colors": [
                            [
                                0,
                                0,
                                255
                            ],
                            [
                                0,
                                255,
                                0
                            ],
                            [
                                255,
                                0,
                                255
                            ],
                            [
                                255,
                                255,
                                0
                            ]
                        ],
                        "contrastLimits": [
                            [
                                1417,
                                45716
                            ],
                            [
                                1424,
                                65534
                            ],
                            [
                                1,
                                6029
                            ],
                            [
                                1,
                                30251
                            ]
                        ],
                        "brightness": [
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.5
                        ],
                        "contrast": [
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.5
                        ],
                        "domains": [
                            [
                                1367,
                                51707
                            ],
                            [
                                1273,
                                65535
                            ],
                            [
                                0,
                                6029
                            ],
                            [
                                0,
                                30251
                            ]
                        ],
                        "selections": [
                            {
                                "c": 0,
                                "t": 0,
                                "z": 0
                            },
                            {
                                "c": 1,
                                "t": 0,
                                "z": 0
                            },
                            {
                                "c": 2,
                                "t": 0,
                                "z": 0
                            },
                            {
                                "c": 3,
                                "t": 0,
                                "z": 0
                            }
                        ]
                    },
                    "viewerStore": {
                        "viewState": {
                            "target": [
                                2752.5072108998293,
                                4267.595596456018
                            ],
                            "zoom": -3.372198311959078
                        }
                    },
                    "imageSettingsStore": {}
                },
                "channel": 0,
                "contourParameter": "leiden",
                "id": "Uzvm8X",
                "size": [
                    1510,
                    767
                ],
                "gsposition": [
                    0,
                    0
                ],
                "gssize": [
                    12,
                    5
                ]
            }
        ]
    },
    "dataSources": {
        "cells": {
            "layout": "gridstack",
            "panelWidth": 100
        }
    }
}
'''

# Function to replace image_name and ds_name
def create_image_view_prototype(new_ds_name, new_image_name):
    updated_json_str = json_string.replace(default_image_name, new_image_name).replace(default_ds, new_ds_name)
    
    # Parse the updated string back into a JSON object
    updated_json_obj = json.loads(updated_json_str)
    
    return updated_json_obj

# Example usage
new_ds_name = "new_cells"
new_image_name = "Vaccinated-Sample_32_ROI_2"

updated_json = create_image_view_prototype(new_ds_name, new_image_name)

# Print the updated JSON object
print(json.dumps(updated_json, indent=4))
