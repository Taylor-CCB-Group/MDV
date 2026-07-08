### UCSC Genome Browser

This chart displays an image from a UCSC session at the specified location 
and requires a session URL ,which can be copied from the browser's navigation bar


Many approaches were tried e.g. iframes and rendering an image directly from the UCSC
genome server, but due to security restriction e.g. CORS, none were successful. Hence a
proxy is required that will pass on the genome co-ordinates to the
UCSC genome browser instance and return the genome image. By default the chart uses
the app/API root `/ucsc_proxy` endpoint (resolved via `mdv_api_root` in the frontend).
Otherwise one can be specified in the datasource in genome.ucsc_proxy_url
