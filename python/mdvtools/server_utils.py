"""
Utility functions for Flask-based servers in the MDV project.

This module provides common server utilities including security headers,
file serving helpers, and HTTP range request handling.
"""

import re
import os
import mimetypes
from zlib import adler32
from typing import BinaryIO
from flask import request, send_file as sf, Response
from werkzeug.exceptions import RequestedRangeNotSatisfiable
from werkzeug.wrappers import Response as WerkzeugResponse


# Entry files of builds that don't hash them (see entry_assets.py). Matched by suffix
# because the static URL prefix varies: /static/ by default, /flask/ in the Docker app.
UNCACHED_STATIC_ENTRY_PATHS = {
    "/js/mdv.js",
    "/js/catalog.js",
    "/js/login.js",
    "/assets/mdv.css",
    "/assets/catalog.css",
}


def _is_uncached_static_entry(path: str) -> bool:
    return any(path.endswith(entry_path) for entry_path in UNCACHED_STATIC_ENTRY_PATHS)


def _set_shell_cache_headers(resp):
    if resp.mimetype == "text/html":
        resp.headers["Cache-Control"] = "no-store"
        return

    if _is_uncached_static_entry(request.path):
        resp.headers["Cache-Control"] = "no-cache, max-age=0, must-revalidate"

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
    # range requests for parquet (& OME-TIFF/other things, but noticing issues with points parquet)
    resp.headers["Access-Control-Allow-Headers"] = "Content-Type, Range"
    resp.headers["Access-Control-Expose-Headers"] = "Content-Range, Content-Length, Accept-Ranges"
    resp.headers["Access-Control-Max-Age"] = "86400"
    #required for vite dev
    resp.headers["Cross-Origin-Resource-Policy"] ="cross-origin"
    _set_shell_cache_headers(resp)
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


class _SeekableBlocks:
    """A file that yields fixed-size blocks and can be seeked.

    Werkzeug only seeks to the start of a range when the iterable reports `seekable()`.
    The `wsgi.file_wrapper` that gunicorn supplies does not, so werkzeug reads and
    discards everything before the range instead, which makes a request near the end
    of a multi-GB zip take hundreds of milliseconds.
    """

    def __init__(self, file: BinaryIO, block_size: int = 64 * 1024) -> None:
        self._file = file
        self._block_size = block_size

    def seekable(self) -> bool:
        return True

    def seek(self, offset: int) -> int:
        return self._file.seek(offset)

    def tell(self) -> int:
        return self._file.tell()

    def __iter__(self) -> "_SeekableBlocks":
        return self

    def __next__(self) -> bytes:
        block = self._file.read(self._block_size)
        if not block:
            raise StopIteration
        return block

    def close(self) -> None:
        self._file.close()


def send_file_ranged(f: str) -> WerkzeugResponse:
    """
    Like send_file, but a Range request seeks straight to the requested bytes
    under any WSGI server. Conditional (ETag / If-Range) handling is unchanged.
    """
    if "Range" not in request.headers:
        return send_file(f)
    path = os.path.abspath(f)
    # One handle supplies the validators and the body, so a file replaced mid-request
    # can't pair its old ETag and length with new bytes.
    file = open(path, "rb")
    try:
        stat = os.fstat(file.fileno())
        response = sf(
            file,
            download_name=os.path.basename(path),
            # same value werkzeug derives from a path, so ETags agree with non-range responses
            etag=f"{stat.st_mtime}-{stat.st_size}-{adler32(path.encode()) & 0xFFFFFFFF}",
            last_modified=stat.st_mtime,
            conditional=False,
        )
        response.content_length = stat.st_size
        response.response = _SeekableBlocks(file)
        response.make_conditional(
            request.environ, accept_ranges=True, complete_length=stat.st_size
        )
    except RequestedRangeNotSatisfiable as e:
        # returned rather than raised so the 416 survives the project router's catch-all
        file.close()
        return Response(status=416, headers={"Content-Range": f"bytes */{e.length}"})
    except Exception:
        file.close()
        raise
    return response


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
