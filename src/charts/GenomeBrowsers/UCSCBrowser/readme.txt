### UCSC Genome Browser

This chart displays an image from a UCSC session at the specified location 
and requires a session URL ,which can be copied from the browser's navigation bar


Many approaches were tried e.g. iframes and rendering an image directly from the UCSC
genome server, but due to security restriction e.g. CORS, none were successful. Hence a
proxy is required that will pass on the genome co-ordinates to the
UCSC genome browser instanceand return the genome image. The serverlite component
contains such a proxy at  /ucsc_proxy, which is the default value used by this chart.
Otherwise one can be specified in the datasource in genome.ucsc_proxy_url
