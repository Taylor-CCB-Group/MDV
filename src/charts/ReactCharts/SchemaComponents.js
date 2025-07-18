


const baseChartProperties = {
    "legend": { "type": "string"},
    "title": { "type": "string"},
    "type": { "type": "string"},
    "size": {
        "type": "array",
        "items": [
            { "type": "integer", "default": 500 },
            { "type": "integer", "default": 500 }
        ],
        "minItems": 2,
        "maxItems": 2
    },
    "position": {
        "type": "array",
        "items": [
            { "type": "integer", "default": 5 },
            { "type": "integer", "default": 5 }
        ],
        "minItems": 2,
        "maxItems": 2
    },
    "experimental": { "type": "boolean", "default": true },
}

function getColumnEntry(type,number){
    return {
      "type":"object",
        "title":"__column__",
        "properties":{
            "type": type,
            "name": number,
        }
    }
}


export {getColumnEntry, baseChartProperties};