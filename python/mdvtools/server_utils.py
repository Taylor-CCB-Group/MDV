"""
Utility functions for Flask-based servers in the MDV project.

This module provides common server utilities including security headers,
file serving helpers, and HTTP range request handling.
"""

import re
import os
import mimetypes
from flask import send_file as sf, Response

# consider using flask_cors...
def add_safe_headers(resp):
    """
    Add security and CORS headers to a Flask response.

    Args:
        resp: Flask response object to modify

    Returns:
        Modified response object with added headers
    """
    # headers required for web workers
    resp.headers["Cross-Origin-Opener-Policy"] = "same-origin"
    resp.headers["Cross-Origin-Embedder-Policy"] = "require-corp"
    # headers required if serving endpoints for another server e,g dev server
    resp.headers["Access-Control-Allow-Origin"] = "*"
    resp.headers["Access-Control-Allow-Headers"] = "Content-Type"
    #required for vite dev
    resp.headers["Cross-Origin-Resource-Policy"] ="cross-origin"
    return resp


# flask send_file can't always cope with relative paths
# sets the cwd to the python path for some reason
def send_file(f):
    """
    Send a file using Flask's send_file with absolute path handling.

    Args:
        f (str): File path (relative or absolute)

    Returns:
        Flask response object for the file
    """
    if not os.path.isabs(f):
        f = os.path.join(os.getcwd(), f)
    return sf(f)


def get_range(file_name, range_header):
    """
    Handle HTTP range requests for partial file content.

    Args:
        file_name (str): Path to the file to serve
        range_header (str): HTTP Range header value

    Returns:
        Flask Response object with partial content (status 206)

    Raises:
        ValueError: If the range header is invalid
    """
  
    size = os.path.getsize(file_name)
    byte1, byte2 = 0, None

    m = re.search(r"(\d+)-(\d*)", range_header)
    if not m:
        raise ValueError("Invalid Range Header")
    g = m.groups()

    if g[0]:
        byte1 = int(g[0])
    if g[1]:
        byte2 = int(g[1])

    length = size - byte1
    if byte2 is not None:
        length = byte2 - byte1 + 1

    with open(file_name, "rb") as file:
        file.seek(byte1)
        data = file.read(length)
        
    rv = Response(
        data, 206, mimetype=mimetypes.guess_type(file_name)[0], direct_passthrough=True
    )
    rv.headers.add(
        "Content-Range", f"bytes {byte1}-{byte1 + length - 1}/{size}"
    )
    rv.headers.add("Accept-Ranges", "bytes")
    return rv