/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


import {jszlib_inflate_buffer,arrayCopy} from "./vendor/inflate.js";
import {Zlib} from "./vendor/zlib_and_gzip.js";

let is_node=false;
try{
    navigator;
}catch(e){
    is_node=true;
}



const NONE = 0;
const GZIP = 1;
const BGZF = 2;
class igvxhr {
   
    // Compression types
   

    static load(url, options) {

       

        return new Promise(function (fulfill, reject) {
            var xhr = new XMLHttpRequest(),
                sendData = options.sendData,
                method = options.method || (sendData ? "POST" : "GET"),
                range = options.range,
                responseType = options.responseType,
                contentType = options.contentType,
                mimeType = options.mimeType,
                headers = options.headers,
                isSafari = is_node?false:navigator.vendor.indexOf("Apple") == 0 && /\sSafari\//.test(navigator.userAgent),
                withCredentials = options.withCredentials,
                header_keys, key, value, i;

            // Support for GCS paths.
           //url = url.startsWith("gs://") ? igv.Google.translateGoogleCloudURL(url) : url;
        

           /* if (igv.Google.isGoogleURL(url)) {

                url = igv.Google.addApiKey(url);

                // Add google headers (e.g. oAuth)
                headers = headers || {};
                igv.Google.addGoogleHeaders(headers);

                // Hack to prevent caching for google storage files.  Get weird net:err-cache errors otherwise
                if (range) {
                    url += url.includes("?") ? "&" : "?";
                    url += "someRandomSeed=" + Math.random().toString(36);
                }
            }
            */
    

            xhr.open(method, url);

            if (range) {
                var rangeEnd = range.size ? range.start + range.size - 1 : "";
                xhr.setRequestHeader("Range", "bytes=" + range.start + "-" + rangeEnd);
            }
            if (contentType) {
                xhr.setRequestHeader("Content-Type", contentType);
            }
            if (mimeType) {
                xhr.overrideMimeType(mimeType);
            }
            if (responseType) {
                xhr.responseType = responseType;
            }
            if (headers) {
                header_keys = Object.keys(headers);
                for (i = 0; i < header_keys.length; i++) {
                    key = header_keys[i];
                    value = headers[key];
                    // console.log("Adding to header: " + key + "=" + value);
                    xhr.setRequestHeader(key, value);
                }
            }

            // NOTE: using withCredentials with servers that return "*" for access-allowed-origin will fail
            if (withCredentials === true) {
                xhr.withCredentials = true;
            }
            xhr.timeout=30000;

            xhr.onload = function (event) {
                // when the url points to a local file, the status is 0 but that is no error
                if (xhr.status == 0 || (xhr.status >= 200 && xhr.status <= 300)) {

                    if (range && xhr.status != 206) {
                        handleError("ERROR: range-byte header was ignored for url: " + url);
                    }
                    else {
                      
                        fulfill(xhr.response,xhr);
                     
                    }
                }
                else {

                    //
                    if (xhr.status === 416) {
                        //  Tried to read off the end of the file.   This shouldn't happen, but if it does return an
                        handleError("Unsatisfiable range");
                    }
                    else {// TODO -- better error handling
                        handleError("Error accessing resource: " + xhr.status);
                    }

                }

            };

            xhr.onerror = function (event) {
               handleError("Error accessing resource: " + url + " Status: " + xhr.status);    
            }


            xhr.ontimeout = function (event) {
                handleError("Timed out");
            };

            xhr.onabort = function (event) {
                console.log("Aborted");
                reject(new igv.AbortLoad());
            };

            try {
               
                xhr.send(sendData);
                
            } catch (e) {
                reject(e);
            }


            function handleError(message) {
                if (reject) {
                    reject(message);
                }
                else {
                    throw Error(message);
                }
            }
        });
    }

    static loadArrayBuffer (url, options) {

        if (options === undefined) options = {};
        options.responseType = "arraybuffer";
        return igvxhr.load(url, options);
    };

    static loadJson (url, options) {

        var method = options.method || (options.sendData ? "POST" : "GET");

        if (method == "POST") options.contentType = "application/json";

        return new Promise(function (fulfill, reject) {

            igvxhr.load(url, options).then(
                function (result) {
                    if (result) {
                        fulfill(JSON.parse(result));
                    }
                    else {
                        fulfill(result);
                    }
                }).catch(reject);
        })
    }

    /**
     * Load a "raw" string.
     */
    static loadString(url, options) {

        var compression, fn, idx;

        if (options === undefined) options = {};

        // Strip parameters from url
        // TODO -- handle local files with ?
        idx = url.indexOf("?");
        fn = idx > 0 ? url.substring(0, idx) : url;

        if (options.bgz) {
            compression = BGZF;
        }
        else if (fn.endsWith(".gz")) {
            compression = GZIP;
        }
        else {
            compression = NONE;
        }

        if (compression === NONE) {
            options.mimeType = 'text/plain; charset=x-user-defined';
            return igvxhr.load(url, options);
        }
        else {
            options.responseType = "arraybuffer";

            return new Promise(function (fulfill, reject) {

                igvxhr.load(url, options).then(
                    function (data) {
                        var result = igvxhr.arrayBufferToString(data, compression);
                        fulfill(result);
                    }).catch(reject)
            })
        }

    };

    static loadStringFromFile(localfile, options) {

        return new Promise(function (fulfill, reject) {

            var fileReader = new FileReader(),
                range = options.range;


            fileReader.onload = function (e) {

                var compression, result;

                if (options.bgz) {
                    compression = BGZF;
                }
                else if (localfile.name.endsWith(".gz")) {

                    compression = GZIP;
                }
                else {
                    compression = NONE;
                }

                result = igvxhr.arrayBufferToString(fileReader.result, compression);

                fulfill(result, localfile);

            };

            fileReader.onerror = function (e) {
                console.log("reject uploading local file " + localfile.name);
                reject(null, fileReader);
            };

            fileReader.readAsArrayBuffer(localfile);

        });
    }

    static isCrossDomain(url) {

        var origin = window.location.origin;

        return !url.startsWith(origin);

    }

    static arrayBufferToString (arraybuffer, compression) {

        var plain, inflate;

        if (compression === GZIP) {
            inflate = new Zlib.Gunzip(new Uint8Array(arraybuffer));
            plain = inflate.decompress();
        }
        else if (compression === BGZF) {
            plain = new Uint8Array(unbgzf(arraybuffer));
        }
        else {
            plain = new Uint8Array(arraybuffer);
        }

        var result = "";
        for (var i = 0, len = plain.length; i < len; i++) {
            result = result + String.fromCharCode(plain[i]);
        }
        return result;
    }

}


//**********js/bam/bgzf.js***************************************


const BLOCK_HEADER_LENGTH = 18;
const BLOCK_LENGTH_OFFSET = 16;  // Location in the gzip block of the total block size (actually total block size - 1)
const BLOCK_FOOTER_LENGTH = 8; // Number of bytes that follow the deflated data
const MAX_COMPRESSED_BLOCK_SIZE = 64 * 1024; // We require that a compressed block (including header and footer, be <= this)
const GZIP_OVERHEAD = BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH + 2; // Gzip overhead is the header, the footer, and the block size (encoded as a short).
const GZIP_ID1 = 31;   // Magic number
const GZIP_ID2 = 139;  // Magic number
const GZIP_FLG = 4; // FEXTRA flag means there are optional fields


    // Uncompress data,  assumed to be series of bgzipped blocks
    // Code is based heavily on bam.js, part of the Dalliance Genome Explorer,  (c) Thomas Down 2006-2001.
let unbgzf = function (data, lim) {

        var oBlockList = [],
            ptr = [0],
            totalSize = 0;

        lim = lim || data.byteLength - 18;
        const blockLengths=[];

        while (ptr[0] < lim) {

            var ba = new Uint8Array(data, ptr[0], 18);

            var xlen = (ba[11] << 8) | (ba[10]);
            var si1 = ba[12];
            var si2 = ba[13];
            var slen = (ba[15] << 8) | (ba[14]);
            var bsize = (ba[17] << 8) | (ba[16]) + 1;
            blockLengths.push(bsize)

            var start = 12 + xlen + ptr[0];    // Start of CDATA
            var length = data.byteLength - start;

            if (length < (bsize + 8)) break;

            var unc = jszlib_inflate_buffer(data, start, length, ptr);

            ptr[0] += 8;    // Skipping CRC-32 and size of uncompressed data

            totalSize += unc.byteLength;
            oBlockList.push(unc);
        }

        // Concatenate decompressed blocks
        if (oBlockList.length == 1) {
            oBlockList.blockLengths=blockLengths;
            return oBlockList[0];
        } else {
            var out = new Uint8Array(totalSize);
            var cursor = 0;
            for (var i = 0; i < oBlockList.length; ++i) {
                var b = new Uint8Array(oBlockList[i]);
                arrayCopy(b, 0, out, cursor, b.length);
                cursor += b.length;
            }
            out.buffer.blockLengths=blockLengths;
            return out.buffer;
        }
}




export {igvxhr,unbgzf};

