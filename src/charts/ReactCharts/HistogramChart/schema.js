
import {getColumnEntry,baseChartProperties } from "../SchemaComponents";


export default  {
  "title": "Histogram",
  "type": "object",
  "properties": {
    ...baseChartProperties,
    "x":getColumnEntry("number","Frequency Data"),
    },
    "required":["x","type"]
}

