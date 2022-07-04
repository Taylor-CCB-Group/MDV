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

import {igvxhr} from "./igvxhr.js";
//import {pako} from "./vendor/pako_inflate.js";
import {Zlib} from "./vendor/zlib_and_gzip.js";


//***********js/bigwig/bufferedReader.js*****************
class BufferedReader{
    constructor(config, contentLength, bufferSize) {
        this.path = config.url;
        this.contentLength = contentLength;
        this.bufferSize = bufferSize ? bufferSize : 512000;
        this.range = {start: -1, size: -1};
        this.config = config;
    }

    /**
     *
     * @param requestedRange - byte rangeas {start, size}
     * @param fulfill - function to receive result
     * @param asUint8 - optional flag to return result as an UInt8Array
     */
    dataViewForRange(requestedRange, asUint8) {

        var self = this;

        return new Promise(function (fulfill, reject) {
            var hasData = (self.data && (self.range.start <= requestedRange.start) &&
                ((self.range.start + self.range.size) >= (requestedRange.start + requestedRange.size))),
                bufferSize,
                loadRange;

            if (hasData) {
                subbuffer(self, requestedRange, asUint8);
            }
            else {
                // Expand buffer size if needed, but not beyond content length
                bufferSize = Math.max(self.bufferSize, requestedRange.size);

                if (self.contentLength > 0 && requestedRange.start + bufferSize > self.contentLength) {
                    loadRange = {start: requestedRange.start};
                }
                else {
                    loadRange = {start: requestedRange.start, size: bufferSize};
                }

                igvxhr.loadArrayBuffer(self.path, Object.assign(self.config, {range: loadRange}))
                    .then(function (arrayBuffer) {
                    self.data = arrayBuffer;
                    self.range = loadRange;
                    subbuffer(self, requestedRange, asUint8);
                }).catch(reject);

            }


            function subbuffer(bufferedReader, requestedRange, asUint8) {

                var len = bufferedReader.data.byteLength,
                    bufferStart = requestedRange.start - bufferedReader.range.start,
                    result = asUint8 ?
                        new Uint8Array(bufferedReader.data, bufferStart, len - bufferStart) :
                        new DataView(bufferedReader.data, bufferStart, len - bufferStart);
                fulfill(result);
            }
        });

    }

}


//**********js/bigwig/bwSource.js***************
class BWSource{

    constructor(config,create_feature_function) {
        this.reader = new BWReader(config);
        this.bufferedReader = new BufferedReader(config);
        if (!create_feature_function){
            this.create_feature=BWSource.createFeature;
        }
        else{
            this.create_feature=create_feature_function;
        }
    }
    
    /**
	* Creates a panel
	* @param {string} chr - The chromosome
	* @param {int} bpStart The starting postition 
	* @param {int} bpEnd - The end of the region to show
	* @param {boolean} use_existing - If true then the cached feature will be used- only used
	* if the co-oridinates have not changed. Although, the BWreader has a cache, it is sometimes
	* ignored and features are re-fetched for the same region
	* @param {object} data - Should contain pixelWidth- the width of the entire canvas and 
	* bpPerPixel.
	*/
    getFeatures(chr, bpStart, bpEnd,use_existing,data) {
        this.st = new Date().getTime();
        var self = this;
        return new Promise(function (fulfill, reject) {
            if (self.features && use_existing){
                fulfill(self.features);
                return;
            }
            self.reader.getZoomHeaders().then(function (zoomLevelHeaders) {

                // Select a biwig "zoom level" appropriate for the current resolution
                var bwReader = self.reader,
                    bufferedReader = self.bufferedReader,
                    bpp =data.bpPerPixel,
                    zoomLevelHeader=BWSource.zoomLevelForScale(bpp, zoomLevelHeaders),
                    treeOffset
                
               
                if (zoomLevelHeader && bwReader.type==="BigWig") {
                    treeOffset = zoomLevelHeader.indexOffset;
                    self.decodeFunction = BWSource.decodeZoomData;
                } else {
                    treeOffset = bwReader.header.fullIndexOffset;
                    if (bwReader.type === "BigWig") {
                        self.decodeFunction =BWSource.decodeWigData;
                    }
                    else {
                        self.decodeFunction =self.decodeBedData;
                    }
                }

                bwReader.loadRPTree(treeOffset).then(function (rpTree) {

                    var chrIdx = self.reader.chromTree.dictionary[chr];
                    if (chrIdx === undefined) {
                        fulfill(null);
                    }
                    else {

                        rpTree.findLeafItemsOverlapping(chrIdx, bpStart, bpEnd).then(function (leafItems) {

                            var promises = [];

                            if (!leafItems || leafItems.length == 0) fulfill([]);

                            leafItems.forEach(function (item) {

                                promises.push(new Promise(function (fulfill, reject) {
                                    var features = [];

                                    bufferedReader.dataViewForRange({
                                        start: item.dataOffset,
                                        size: item.dataSize
                                    }, true).then(function (uint8Array) {
                                  
                                       var inflate = new Zlib.Inflate(uint8Array);
                                       var plain= inflate.decompress();
                               
                                      
                                       //var plain = pako.inflate(uint8Array);
                                        self.decodeFunction(new DataView(plain.buffer), chr, chrIdx, bpStart, bpEnd, features);

                                        fulfill(features);

                                    }).catch(reject);
                                }));
                            });


                            Promise.all(promises).then(function (featureArrays) {
                                var en = new Date().getTime();
                                var e = en-self.st;
                                var a = bpp;
                                //console.log(e);
                                var i, allFeatures = featureArrays[0];
                                if(featureArrays.length > 1) {
                                   for(i=1; i<featureArrays.length; i++) {
                                       allFeatures = allFeatures.concat(featureArrays[i]);
                                   }
                                }  
                                allFeatures.sort(function (a, b) {
                                    return a.start - b.start;
                                })
                                self.features=allFeatures;
                                fulfill(allFeatures)
                            }).catch(reject);

                        }).catch(function(error){
                            reject(error);
                        });
                    }
                }).catch(function(error){
                    reject(error)
                });
            }).catch(function(error){
                reject(error);
            }
            );


        });
    }
    
    
    getDefaultRange() {
        
        if(this.reader.totalSummary != undefined) {
            return this.reader.totalSummary.defaultRange;
        }
        else {
            return undefined;
        }
        
    }


    static zoomLevelForScale(bpPerPixel, zoomLevelHeaders) {

        var level = null, i, zl;

        for (i = 0; i < zoomLevelHeaders.length; i++) {

            zl = zoomLevelHeaders[i];

            if (zl.reductionLevel > bpPerPixel) {
                level = zl;
                break;
            }
        }

        if (null == level) {
            level = zoomLevelHeaders[zoomLevelHeaders.length - 1];
        }

        return (level && level.reductionLevel < 4 * bpPerPixel) ? level : null;
    }


    static decodeWigData(data, chr, chrIdx, bpStart, bpEnd, featureArray) {

        var binaryParser = new BinaryParser(data),
            chromId = binaryParser.getInt(),
            chromStart = binaryParser.getInt(),
            chromEnd = binaryParser.getInt(),
            itemStep = binaryParser.getInt(),
            itemSpan = binaryParser.getInt(),
            type = binaryParser.getByte(),
            reserved = binaryParser.getByte(),
            itemCount = binaryParser.getUShort(),
            value;

        if (chromId === chrIdx) {

            while (itemCount-- > 0) {

                switch (type) {
                    case 1:
                        chromStart = binaryParser.getInt();
                        chromEnd = binaryParser.getInt();
                        value = binaryParser.getFloat();
                        break;
                    case 2:
                        chromStart = binaryParser.getInt();
                        value = binaryParser.getFloat();
                        chromEnd = chromStart + itemSpan;
                        break;
                    case 3:  // Fixed step
                        value = binaryParser.getFloat();
                        chromEnd = chromStart + itemSpan;
                        chromStart += itemStep;
                        break;

                }

                if (chromStart >= bpEnd) {
                    break; // Out of interval
                } else if (chromEnd > bpStart && Number.isFinite(value)) {
                    featureArray.push({chr: chr, start: chromStart, end: chromEnd, value: value});
                }


            }
        }

    }

    static decodeZoomData(data, chr, chrIdx, bpStart, bpEnd, featureArray) {

        var binaryParser = new BinaryParser(data),
            minSize = 8 * 4,   // Minimum # of bytes required for a zoom record
            chromId,
            chromStart,
            chromEnd,
            validCount,
            minVal,
            maxVal,
            sumData,
            sumSquares,
            value;

        while (binaryParser.remLength() >= minSize) {
            chromId = binaryParser.getInt();
            if (chromId === chrIdx) {

                chromStart = binaryParser.getInt();
                chromEnd = binaryParser.getInt();
                validCount = binaryParser.getInt();
                minVal = binaryParser.getFloat();
                maxVal = binaryParser.getFloat();
                sumData = binaryParser.getFloat();
                sumSquares = binaryParser.getFloat();
                value = validCount == 0 ? 0 : sumData / validCount;

                if (chromStart >= bpEnd && chromStart<1000000000) {
                     console.log("should have broken")

                    break; // Out of interval
                   
                } else if (chromEnd > bpStart && Number.isFinite(value)) {
                    featureArray.push({chr: chr, start: chromStart, end: chromEnd, value: value});
                }

            }
        }

    }




    decodeBedData(data, chr, chrIdx, bpStart, bpEnd, featureArray) {

        var binaryParser = new BinaryParser(data),
            minSize = 3 * 4 + 1,   // Minimum # of bytes required for a bed record
            chromId,
            chromStart,
            chromEnd,
            rest,
            tokens,
            feature,
            exonCount, exonSizes, exonStarts, exons, eStart, eEnd;


        while (binaryParser.remLength() >= minSize) {

            chromId = binaryParser.getInt();
            if (chromId != chrIdx) continue;

            chromStart = binaryParser.getInt();
            chromEnd = binaryParser.getInt();
            rest = binaryParser.getString();

            feature = {chr: chr, start: chromStart, end: chromEnd};

            if (chromStart < bpEnd && chromEnd >= bpStart) {
                featureArray.push(feature);

                tokens = rest.split("\t");
                this.create_feature(tokens,feature);

                
            }
        }

    }
    
    static createFeature(tokens,feature){
        if (tokens.length > 0) {
                    feature.name = tokens[0];
                }

                if (tokens.length > 1) {
                    feature.score = parseFloat(tokens[1]);
                }
                if (tokens.length > 2) {
                    feature.strand = tokens[2];
                }
                if (tokens.length > 3) {
                    feature.cdStart = parseInt(tokens[3]);
                }
                if (tokens.length > 4) {
                    feature.cdEnd = parseInt(tokens[4]);
                }
                if (tokens.length > 5) {
                    //if (tokens[5] !== "." && tokens[5] !== "0")
                        //feature.color = igv.createColorString(tokens[5]);
                }
                if (tokens.length > 8) {
                    exonCount = parseInt(tokens[6]);
                    exonSizes = tokens[7].split(',');
                    exonStarts = tokens[8].split(',');
                    exons = [];

                    for (var i = 0; i < exonCount; i++) {
                        eStart = start + parseInt(exonStarts[i]);
                        eEnd = eStart + parseInt(exonSizes[i]);
                        exons.push({start: eStart, end: eEnd});
                    }

                    feature.exons = exons;
                }

    }
}



//************js/bigwig/bwReader.js*****************


const BIGWIG_MAGIC_LTH = 0x888FFC26; // BigWig Magic Low to High
const BIGWIG_MAGIC_HTL = 0x26FC8F66; // BigWig Magic High to Low
const BIGBED_MAGIC_LTH = 0x8789F2EB; // BigBed Magic Low to High
const BIGBED_MAGIC_HTL = 0xEBF28987; // BigBed Magic High to Low
const BBFILE_HEADER_SIZE = 64;


class BWReader{
    constructor(config) {
        this.path = config.url;
        this.headPath = config.headURL || this.path;
        this.rpTreeCache = {};
        this.config = JSON.parse(JSON.stringify(config))
    };

    getZoomHeaders() {

        var self = this;

        return new Promise(function (fulfill, reject) {
            if (self.zoomLevelHeaders) {
                fulfill(self.zoomLevelHeaders);
            }
            else {
                self.loadHeader().then(function () {
                    fulfill(self.zoomLevelHeaders);
                }).catch(function (error) {
                    reject(error);
                });
            }
        });
    }

    loadHeader() {

        var self = this;

        return new Promise(function (fulfill, reject) {
            igvxhr.loadArrayBuffer(self.path, Object.assign(self.config, {range: {start: 0, size: BBFILE_HEADER_SIZE}}))
                .then(function (data) {

                if (!data) return;

                // Assume low-to-high unless proven otherwise
                self.littleEndian = true;

                var binaryParser = new BinaryParser(new DataView(data));

                var magic = binaryParser.getUInt();

                if (magic === BIGWIG_MAGIC_LTH) {
                    self.type = "BigWig";
                }
                else if (magic == BIGBED_MAGIC_LTH) {
                    self.type = "BigBed";
                }
                else {
                    //Try big endian order
                    self.littleEndian = false;

                    binaryParser.littleEndian = false;
                    binaryParser.position = 0;
                    var magic = binaryParser.getUInt();

                    if (magic === BIGWIG_MAGIC_HTL) {
                        self.type = "BigWig";
                    }
                    else if (magic == BIGBED_MAGIC_HTL) {
                        self.type = "BigBed";
                    }
                    else {
                        // TODO -- error, unknown file type  or BE
                    }

                }
                // Table 5  "Common header for BigWig and BigBed files"
                self.header = {};
                self.header.bwVersion = binaryParser.getUShort();
                self.header.nZoomLevels = binaryParser.getUShort();
                self.header.chromTreeOffset = binaryParser.getLong();
                self.header.fullDataOffset = binaryParser.getLong();
                self.header.fullIndexOffset = binaryParser.getLong();
                self.header.fieldCount = binaryParser.getUShort();
                self.header.definedFieldCount = binaryParser.getUShort();
                self.header.autoSqlOffset = binaryParser.getLong();
                self.header.totalSummaryOffset = binaryParser.getLong();
                self.header.uncompressBuffSize = binaryParser.getInt();
                self.header.reserved = binaryParser.getLong();

                self.loadZoomHeadersAndChrTree().then(fulfill).catch(reject);
            }).catch(function (error) {
                    reject(error);
                });

        });
    }


   loadZoomHeadersAndChrTree() {


        var startOffset = BBFILE_HEADER_SIZE,
            self = this;

        return new Promise(function (fulfill, reject) {
            
            var range = {start: startOffset, size: (self.header.fullDataOffset - startOffset + 5)};

            igvxhr.loadArrayBuffer(self.path, Object.assign(self.config, {range: range}))
                .then(function (data) {

                var nZooms = self.header.nZoomLevels,
                    binaryParser = new BinaryParser(new DataView(data)),
                    i,
                    len,
                    zoomNumber,
                    zlh;

                self.zoomLevelHeaders = [];

                self.firstZoomDataOffset = Number.MAX_VALUE;
                for (i = 0; i < nZooms; i++) {
                    zoomNumber = nZooms - i;
                    zlh = new ZoomLevelHeader(zoomNumber, binaryParser);
                    self.firstZoomDataOffset = Math.min(zlh.dataOffset, self.firstZoomDataOffset);
                    self.zoomLevelHeaders.push(zlh);
                }

                // Autosql
                if (self.header.autoSqlOffset > 0) {
                    binaryParser.position = self.header.autoSqlOffset - startOffset;
                    self.autoSql = binaryParser.getString();
                }

                // Total summary
                if (self.header.totalSummaryOffset > 0) {
                    binaryParser.position = self.header.totalSummaryOffset - startOffset;
                    self.totalSummary = new BWTotalSummary(binaryParser);
                }

                // Chrom data index
                if (self.header.chromTreeOffset > 0) {
                    binaryParser.position = self.header.chromTreeOffset - startOffset;
                    self.chromTree = new BPTree(binaryParser, startOffset);
                }
                else {
                    // TODO -- this is an error, not expected
                }

                //Finally total data count
                binaryParser.position = self.header.fullDataOffset - startOffset;
                self.dataCount = binaryParser.getInt();

                fulfill();

            }).catch(function(error){
                reject(error);
            });
        });
    }

    loadRPTree(offset) {

        var self = this;

        return new Promise(function (fulfill, reject) {
            var rpTree = self.rpTreeCache[offset];
            if (rpTree) {
                fulfill(rpTree);
            }
            else {
                rpTree = new RPTree(offset, self.contentLength, self.config, self.littleEndian);
                self.rpTreeCache[offset] = rpTree;
                rpTree.load().then(function () {
                    fulfill(rpTree);
                }).catch(reject);
            }
        });
    }
}

 class ZoomLevelHeader{
     constructor(index, byteBuffer) {
        this.index = index;
        this.reductionLevel = byteBuffer.getInt();
        this.reserved = byteBuffer.getInt();
        this.dataOffset = byteBuffer.getLong();
        this.indexOffset = byteBuffer.getLong();
    }
 }


const RPTREE_MAGIC_LTH = 0x2468ACE0;
const RPTREE_MAGIC_HTL = 0xE0AC6824;
const RPTREE_HEADER_SIZE = 48;
const RPTREE_NODE_LEAF_ITEM_SIZE = 32;   // leaf item size
const RPTREE_NODE_CHILD_ITEM_SIZE = 24;  // child item size
const BUFFER_SIZE = 512000; 

//***********js/bigwig/RPTree*******************


     //  buffer

class RPTree{

    constructor (fileOffset, contentLength, config, littleEndian) {

        this.config = config;
        this.filesize = contentLength;
        this.fileOffset = fileOffset; // File offset to beginning of tree
        this.path = config.url;
        this.littleEndian = littleEndian;
    }


    load() {

        var self = this;

        return new Promise(function (fulfill, reject) {
            var rootNodeOffset = self.fileOffset + RPTREE_HEADER_SIZE,
                bufferedReader = new BufferedReader(self.config, self.filesize, BUFFER_SIZE);

            self.readNode(rootNodeOffset, bufferedReader).then(function (node) {
                self.rootNode = node;
                fulfill(self);
            }).catch(reject);
        });
    }


    readNode(filePosition, bufferedReader) {

        var self = this;

        return new Promise(function (fulfill, reject) {

            bufferedReader.dataViewForRange({start: filePosition, size: 4}, false).then(function (dataView) {
                var binaryParser = new BinaryParser(dataView, self.littleEndian);

                var type = binaryParser.getByte();
                var isLeaf = (type === 1) ? true : false;
                var reserved = binaryParser.getByte();
                var count = binaryParser.getUShort();

                filePosition += 4;

                var bytesRequired = count * (isLeaf ? RPTREE_NODE_LEAF_ITEM_SIZE : RPTREE_NODE_CHILD_ITEM_SIZE);
                var range2 = {start: filePosition, size: bytesRequired};

                bufferedReader.dataViewForRange(range2, false).then(function (dataView) {

                    var i,
                        items = new Array(count),
                        binaryParser = new BinaryParser(dataView);

                    if (isLeaf) {
                        for (i = 0; i < count; i++) {
                            var item = {
                                isLeaf: true,
                                startChrom: binaryParser.getInt(),
                                startBase: binaryParser.getInt(),
                                endChrom: binaryParser.getInt(),
                                endBase: binaryParser.getInt(),
                                dataOffset: binaryParser.getLong(),
                                dataSize: binaryParser.getLong()
                            };
                            items[i] = item;

                        }
                        fulfill(new RPTreeNode(items));
                    }
                    else { // non-leaf
                        for (i = 0; i < count; i++) {

                            var item = {
                                isLeaf: false,
                                startChrom: binaryParser.getInt(),
                                startBase: binaryParser.getInt(),
                                endChrom: binaryParser.getInt(),
                                endBase: binaryParser.getInt(),
                                childOffset: binaryParser.getLong()
                            };
                            items[i] = item;

                        }

                        fulfill(new RPTreeNode(items));
                    }
                }).catch(reject);
            }).catch(reject);
        });
    }


    findLeafItemsOverlapping(chrIdx, startBase, endBase) {

        var self = this;

        return new Promise(function (fulfill, reject) {

            var leafItems = [],
                processing = new Set(),
                bufferedReader = new BufferedReader(self.config, self.filesize, BUFFER_SIZE);

            processing.add(0);  // Zero represents the root node
            findLeafItems(self.rootNode, 0);

            function findLeafItems(node, nodeId) {

                if (RPTree.overlaps(node, chrIdx, startBase, endBase)) {

                    var items = node.items;

                    items.forEach(function (item) {

                        if (RPTree.overlaps(item, chrIdx, startBase, endBase)) {

                            if (item.isLeaf) {
                                leafItems.push(item);
                            }

                            else {
                                if (item.childNode) {
                                    findLeafItems(item.childNode);
                                }
                                else {
                                    processing.add(item.childOffset);  // Represent node to-be-loaded by its file position
                                    self.readNode(item.childOffset, bufferedReader).then(function (node) {
                                        item.childNode = node;
                                        findLeafItems(node, item.childOffset);
                                    }).catch(reject);
                                }
                            }
                        }
                    });

                }

                if (nodeId != undefined) processing.delete(nodeId);

                // Wait until all nodes are processed
                if (processing.size===0) {
                    fulfill(leafItems);
                }
            }
        });
    }




    /**
     * Return true if {chrIdx:startBase-endBase} overlaps item's interval
     * @returns {boolean}
     */
    static overlaps(item, chrIdx, startBase, endBase) {

        //  if (chrIdx > item.endChrom || chrIdx < item.startChrom) return false;

        if (!item) {
            console.log("null item");
            return false;
        }

        return ((chrIdx > item.startChrom) || (chrIdx == item.startChrom && endBase >= item.startBase)) &&
            ((chrIdx < item.endChrom) || (chrIdx == item.endChrom && startBase < item.endBase));
    }
}


class RPTreeNode{
    constructor(items) {
        this.items = items;

        var minChromId = Number.MAX_VALUE,
            maxChromId = 0,
            minStartBase = Number.MAX_VALUE,
            maxEndBase = 0,
            i,
            item;

        for (i = 0; i < items.length; i++) {
            item = items[i];
            minChromId = Math.min(minChromId, item.startChrom);
            maxChromId = Math.max(maxChromId, item.endChrom);
            minStartBase = Math.min(minStartBase, item.startBase);
            maxEndBase = Math.max(maxEndBase, item.endBase);
        }

        this.startChrom = minChromId;
        this.endChrom = maxChromId;
        this.startBase = minStartBase;
        this.endBase = maxEndBase;

    }
}

//*******************js/binary.js******************
class BinaryParser{
    constructor(dataView, littleEndian) {

        this.littleEndian = (littleEndian ? littleEndian : true);
        this.position = 0;
        this.view = dataView;
        this.length = dataView.byteLength;
    }

    available() {
        return this.length - this.position;
    }

    remLength() {
        return this.length - this.position;
    }

    hasNext() {
        return this.position < this.length - 1;
    }

    getByte() {
        var retValue = this.view.getUint8(this.position, this.littleEndian);
        this.position++;
        return retValue;
    }

    getShort() {

        var retValue = this.view.getInt16(this.position, this.littleEndian);
        this.position += 2
        return retValue;
    }

    getUShort(){

        // var byte1 = this.getByte(),
        //     byte2 = this.getByte(),
        //     retValue = ((byte2 << 24 >>> 16) + (byte1 << 24 >>> 24));
        //     return retValue;

       //
        var retValue = this.view.getUint16 (this.position, this.littleEndian);
        this.position += 2
        return retValue;
    }


    getInt() {

        var retValue = this.view.getInt32(this.position, this.littleEndian);
        this.position += 4;
        return retValue;
    }


    getUInt() {
        var retValue = this.view.getUint32(this.position, this.littleEndian);
        this.position += 4;
        return retValue;
    }

    getLong() {

        // DataView doesn't support long. So we'll try manually

        var b = [];
        b[0] = this.view.getUint8(this.position);
        b[1] = this.view.getUint8(this.position + 1);
        b[2] = this.view.getUint8(this.position + 2);
        b[3] = this.view.getUint8(this.position + 3);
        b[4] = this.view.getUint8(this.position + 4);
        b[5] = this.view.getUint8(this.position + 5);
        b[6] = this.view.getUint8(this.position + 6);
        b[7] = this.view.getUint8(this.position + 7);

        var value = 0;
        if (this.littleEndian) {
            for (var i = b.length - 1; i >= 0; i--) {
                value = (value * 256) + b[i];
            }
        } else {
            for (var i = 0; i < b.length; i++) {
                value = (value * 256) + b[i];
            }
        }


        this.position += 8;
        return value;
    }

    getString(len) {

        var s = "";
        var c;
        while ((c = this.view.getUint8(this.position++)) != 0) {
            s += String.fromCharCode(c);
            if (len && s.length == len) break;
        }
        return s;
    }

    getFixedLengthString(len) {

        var s = "";
        var i;
        var c;
        for (i = 0; i < len; i++) {
            c = this.view.getUint8(this.position++);
            if (c > 0) {
                s += String.fromCharCode(c);
            }
        }
        return s;
    }

    getFixedLengthTrimmedString(len) {

        var s = "";
        var i;
        var c;
        for (i = 0; i < len; i++) {
            c = this.view.getUint8(this.position++);
            if (c > 32) {
                s += String.fromCharCode(c);
            }
        }
        return s;
    }

    getFloat() {

        var retValue = this.view.getFloat32(this.position, this.littleEndian);
        this.position += 4;
        return retValue;


    }

    getDouble() {

        var retValue = this.view.getFloat64(this.position, this.littleEndian);
        this.position += 8;
        return retValue;
    }

    skip(n) {

        this.position += n;
        return this.position;
    }


    /**
     * Return a bgzip (bam and tabix) virtual pointer
     * TODO -- why isn't 8th byte used ?
     * @returns {*}
     */
    getVPointer() {

        var position = this.position,
            offset = (this.view.getUint8(position + 1) << 8) | (this.view.getUint8(position)),
            byte6 = ((this.view.getUint8(position + 6) & 0xff) * 0x100000000),
            byte5 = ((this.view.getUint8(position + 5) & 0xff) * 0x1000000),
            byte4 = ((this.view.getUint8(position + 4) & 0xff) * 0x10000),
            byte3 = ((this.view.getUint8(position + 3) & 0xff) * 0x100),
            byte2 = ((this.view.getUint8(position + 2) & 0xff)),
            block = byte6 + byte5 + byte4 + byte3 + byte2;
        this.position += 8;

        //       if (block == 0 && offset == 0) {
        //           return null;
        //       } else {
        return new VPointer(block, offset);
        //       }
    }
}

class VPointer{
    constructor(block, offset) {
        this.block = block;
        this.offset = offset;
    }

    isLessThan(vp) {
        return this.block < vp.block ||
            (this.block === vp.block && this.offset < vp.offset);
    }

    isGreaterThan(vp) {
        return this.block > vp.block ||
            (this.block === vp.block && this.offset > vp.offset);
    }

    print() {
        return "" + this.block + ":" + this.offset;
    }
}


//*******js/bigwig/bwTotalSummary.js*************



class BWTotalSummary{
    constructor(byteBuffer) {

        if (byteBuffer) {

            this.basesCovered = byteBuffer.getLong();
            this.minVal = byteBuffer.getDouble();
            this.maxVal = byteBuffer.getDouble();
            this.sumData = byteBuffer.getDouble();
            this.sumSquares = byteBuffer.getDouble();

            this.computeStats();
        }
        else {
            this.basesCovered = 0;
            this.minVal = 0;
            this.maxVal = 0;
            this.sumData = 0;
            this.sumSquares = 0;
            this.mean = 0;
            this.stddev = 0;
        }
    }


     computeStats() {
        var n = this.basesCovered;
        if (n > 0) {
            this.mean = this.sumData / n;
            this.stddev = Math.sqrt(this.sumSquares / (n - 1));

            var min = this.minVal < 0 ? this.mean - 2 * this.stddev : 0,
                max = this.maxVal > 0 ? this.mean + 2 * this.stddev : 0;

            this.defaultRange = {
                min: 0,
                max: this.mean + 3 * this.stddev
            }
        }
    }

    updateStats(stats) {

        this.basesCovered += stats.count;
        this.sumData += status.sumData;
        this.sumSquares += sumSquares;
        this.minVal = MIN(_minVal, min);
        this.maxVal = MAX(_maxVal, max);

        computeStats.call(this);

    }


}


//***************js/bigwig/bwBPTree.js**************


const BPTREE_MAGIC_LTH = 0x78CA8C91;
const BPTREE_MAGIC_HTL = 0x918CCA78;
const BPTREE_HEADER_SIZE = 32;


 class BPTree{
     constructor(binaryParser, startOffset) {

        var self = this,
            genome =  null;

        this.header = {};
        this.header.magic = binaryParser.getInt();
        this.header.blockSize = binaryParser.getInt();
        this.header.keySize = binaryParser.getInt();
        this.header.valSize = binaryParser.getInt();
        this.header.itemCount = binaryParser.getLong();
        this.header.reserved = binaryParser.getLong();

        this.dictionary = {};

        // Recursively walk tree to populate dictionary
        readTreeNode(binaryParser, -1, this.header.keySize, this.dictionary);

        var itemSize = 8 + this.header.keySize;
        var minSize = 4 + itemSize;   // Bytes for a node with 1 item

        function readTreeNode(byteBuffer, offset, keySize, dictionary) {

            if (offset >= 0) byteBuffer.position = offset;

            var type = byteBuffer.getByte(),
                reserved = byteBuffer.getByte(),
                count = byteBuffer.getUShort(),
                i,
                key,
                chromId,
                chromSize,
                childOffset,
                bufferOffset,
                currOffset;


            if (type == 1) {

                for (i = 0; i < count; i++) {

                    key = byteBuffer.getFixedLengthTrimmedString(keySize);
                    chromId = byteBuffer.getInt();
                    chromSize = byteBuffer.getInt();

                    if(genome) key = genome.getChromosomeName(key);  // Translate to canonical chr name
                    dictionary[key] = chromId;

                }
            }
            else { // non-leaf

                for (i = 0; i < count; i++) {

                    key = byteBuffer.getFixedLengthTrimmedString(keySize);
                    childOffset = byteBuffer.getLong();
                    bufferOffset = childOffset - startOffset;
                    currOffset = byteBuffer.position;
                    readTreeNode(byteBuffer, bufferOffset, keySize, dictionary);
                    byteBuffer.position = currOffset;
                }
            }

        }
    }


}


export {BWSource,BinaryParser};


