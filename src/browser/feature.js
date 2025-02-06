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


import {Utils} from "./utils.js";
import {igvxhr,unbgzf} from "./igvxhr.js";
import {loadBamIndex} from "./bam.js";
import {BWSource} from "./bigwig.js";
import { getProjectURL } from "../dataloaders/DataLoaderUtil";

const MAX_GZIP_BLOCK_SIZE = (1 << 16);

    /**
     * feature source for "bed like" files (tab delimited files with 1 feature per line: bed, gff, vcf, etc)
     *
     * @class
     * @param config
     */
class FeatureSource{
    constructor(config) {
        this.config = config || {};
    }


    getFileHeader() {
        this.is_indexed=true;  
        return new Promise(function (fulfill, reject) {
            fulfill();        
        });
    }   


    getFeatures(chr,start,end,force,data){
      
        //no need to get header
        if (this.is_indexed){
            return this._getFeatures(chr,start,end,force,data);
        }
        //get headers/index then get features
        var self = this;
        return new Promise(function(fulfill,reject){
            self.getFileHeader().then(function(){
                self._getFeatures(chr,start,end,force,data).then(function(features){
                    fulfill(features);
                }).catch(reject)
            }).catch(reject);
        });
        
    }

    /**
     * Required function fo all data source objects.  Fetches features for the
     * range requested and passes them on to the success function.  Usually this is
     * a function that renders the features on the canvas
     *
     * @param chr
     * @param bpStart
     * @param bpEnd
     */
    _getFeatures(chr, bpStart, bpEnd,force,data) {
        if (bpStart===0){
            bpStart=1;
        }
        var self = this;
        self.time=Date.now();
        return new Promise(function (fulfill, reject) {
            if (self.featureCache && chr !== self.featureCache.range.chr){
                self.featureCache=null;
            }
            var genomicInterval = new GenomicInterval(chr, bpStart, bpEnd),
                featureCache = self.featureCache,
                maxRows = self.config.maxRows || 500;
            let ranges_to_get=false;
            if (!featureCache){
                ranges_to_get={all:[bpStart,bpEnd]};
            }
            else{
                if (featureCache.range !== undefined){
                    ranges_to_get=featureCache.range.rangesToGet(genomicInterval)
                }
            }
            if (!ranges_to_get) {
                fulfill(self.featureCache.queryFeatures(chr, bpStart, bpEnd));
            }
            else{
                let promises=[];
                let p_types=[];
                for (let type in ranges_to_get){
                    let range= ranges_to_get[type];
                    promises.push(self.retrieveFeatures(chr, range[0], range[1],force,data));
                    p_types.push([type,range]);
                }
                Promise.all(promises).then(
                    function (all_features) {
                        let existing_features=[];
                        if (self.featureCache){
                            existing_features=self.featureCache.allFeatures();//featureCache.allFeatures(chr,self.featureCache.range.start,self.featureCache.range.end);
                        }
                        let index=0;
                        for (let featureList of all_features){
                            if (featureList === null){
                                featureList=[];
                            }                        
                            if (p_types[index][0]==="left"){
                                let end = p_types[index][1][1]; 
                                //remove any already retrieved  
                                let splice=0;
                                for (let n=featureList.length-1;n>=0;n--){
                                    if (featureList[n].end>= end){
                                        featureList[n].remove=true;
                                        splice++;
                                    }
                                    
                                }
                                if (splice!==0){
                                    featureList = featureList.filter(x=>!x.remove)
                                }

                            }
                            if (p_types[index][0]==="right"){
                                let start=p_types[index][1][0];
                                //remove any already retrieved
                                let i=0
                                for (i=0;i<featureList.length;i++){
                                    if (featureList[i].start> start){
                                        break;
                                    }
                                }
                                if (i!==0){
                                    featureList.splice(0,i)
                                }
                            }
                            index++;
                            existing_features=existing_features.concat(featureList);
                        }
                        let gi = self.featureCache?self.featureCache.range:genomicInterval;
                        self.featureCache = new FeatureCache(existing_features, gi);
                        FeatureSource.packFeatures(existing_features, maxRows); 
                        fulfill(self.featureCache.queryFeatures(chr, bpStart, bpEnd));
                }).catch(function(error){
                    reject(error);
                });
            }
        });
    }




    static packFeatures(features, maxRows) {

        if (features == null || features.length === 0) {
            return;
        }

        // Segregate by chromosome

        var chrFeatureMap = {},
            chrs = [];
        features.forEach(function (feature) {

            var chr = feature.chr,
                flist = chrFeatureMap[chr];

            if (!flist) {
                flist = [];
                chrFeatureMap[chr] = flist;
                chrs.push(chr);
            }

            flist.push(feature);
        });

        // Loop through chrosomosomes and pack features;

        chrs.forEach(function (chr) {

            pack(chrFeatureMap[chr], maxRows);
        });


        // Assigns a row # to each feature.  If the feature does not fit in any row and #rows == maxRows no
        // row number is assigned.
        function pack(featureList, maxRows) {

            var rows = [];

            featureList.sort(function (a, b) {
                return a.start - b.start;
            })


            rows.push(-1000);
            featureList.forEach(function (feature) {

                var i,
                    r,
                    len = Math.min(rows.length, maxRows),
                    start = feature.start;

                for (r = 0; r < len; r++) {
                    if (start >= rows[r]) {
                        feature.row = r;
                        rows[r] = feature.end;
                        return;
                    }
                }
                feature.row = r;
                rows[r] = feature.end;


            });
        }
    }

}

class TabixBedFeatureSource extends FeatureSource{
    constructor(config,decode_function){
        super(config);
        this.reader = new FeatureFileReader(config,decode_function);
    }

    retrieveFeatures(chr,start,end){
        console.log("getting features",chr,start,end);
        return this.reader.readFeatures(chr,start,end);
    }

    getFileHeader() {
        var self = this;
        return new Promise(function (fulfill, reject) {
            self.reader.readHeader().then(function () {
                self.is_indexed=true;  
                fulfill();
             }).catch(reject);
                
        });
    }   
}

class BigBedFeatureSource extends FeatureSource{
    constructor(config,decode_function){
		config.sourceType="gtex";
		super(config);
		this.header=true;
		this.feature_source=new BWSource(config,decode_function);
	}

	retrieveFeatures(chr,bpStart,bpEnd,force,data){
		return this.feature_source.getFeatures(chr,bpStart,bpEnd,false,data);     	
	}
}

//********js/FeatureFileReader.js*****



const F_MAX_GZIP_BLOCK_SIZE = (1 << 16);

    /**
     * Reader for "bed like" files (tab delimited files with 1 feature per line: bed, gff, vcf, etc)
     *
     * @class
     * @param config
     */
class FeatureFileReader{
    constructor(config,dec_function) {

        this.config = config || {};

        if (config.localFile) {
            this.localFile = config.localFile;
            this.filename = config.localFile.name;
        }
        else {
            //this should be tested in any build configs not used by Peter...
            //probably will change how we do this at some point
            this.url = getProjectURL(config.url, false);
            this.indexURL = config.indexURL;
            this.headURL = config.headURL || this.filename;

            var uriParts = Utils.parseUri(config.url);
            this.filename = uriParts.file;
            this.path = uriParts.path;
        }

        this.format = config.format;

        this.parser = this.getParser(this.format, dec_function);
    };


    getParser(format, decode) {
        switch (format) {
            case "vcf":
                return new VcfParser();
            case "seg" :
                return new SegParser();
            default:
                return new FeatureParser(format, decode, this.config);
        }

    }

    // seg files don't have an index
    isIndexable() {
        var configIndexURL = this.config.indexURL,
            type = this.type,
            configIndexed = this.config.indexed;

        return configIndexURL || (type != "wig" && configIndexed != false);
    }


    /**
     * Return a Promise for the async loaded index
     */
    loadIndex() {
        var idxFile = this.indexURL;
        if (this.filename.endsWith(".gz")) {
            if (!idxFile) idxFile = this.url + ".tbi";
            // idxFile = getProjectURL(idxFile);
            return loadBamIndex(idxFile, this.config, true);
        }
        else {
            if (!idxFile) idxFile = this.url + ".idx";
            return loadTribbleIndex(idxFile, this.config);
        }
    }

    loadFeaturesNoIndex() {

        var self = this;

        return new Promise(function (fulfill, reject) {
            var parser = self.parser,
                options = {
                    headers: self.config.headers,           // http headers, not file header
                    withCredentials: self.config.withCredentials
                };

            if (self.localFile) {
                igvxhr.loadStringFromFile(self.localFile, options).then(parseData).catch(reject);
            }
            else {
                igvxhr.loadString(self.url, options).then(parseData).catch(reject);
            }


            function parseData(data) {
                self.header = parser.parseHeader(data);
                if (self.header instanceof String && self.header.startsWith("##gff-version 3")) {
                    self.format = 'gff3';
                }
                fulfill(parser.parseFeatures(data));   // <= PARSING DONE HERE
            };
        });
    }

 

    loadFeaturesWithIndex(chr, start, end) {
        var self = this;

        return new Promise(function (fulfill, reject) {

            var blocks,
                index = self.index,
                tabix = index && index.tabix,
                refId = tabix ? index.sequenceIndexMap[chr] : chr,
                promises = [];

            //blocks = index.getBlocksForRange(refId, start, end);
            blocks= index.blocksForRange(refId, start, end);

            if (!blocks || blocks.length === 0) {
                fulfill(null);       // TODO -- is this correct?  Should it return an empty array?
            }
            else {
                blocks.forEach(function (block) {
                    if (!block){
                        return;
                    }

                    promises.push(new Promise(function (fulfill, reject) {
                        var startPos = block.minv.block,
                            startOffset = block.minv.offset,
                            endPos = block.maxv.block + (index.tabix ? F_MAX_GZIP_BLOCK_SIZE : 0),
                            options = {
                                headers: self.config.headers,           // http headers, not file header
                                range: {start: startPos, size: endPos - startPos + 1},
                                withCredentials: self.config.withCredentials
                            },
                            success;

                        success = function (data) {
                           

                            var inflated, slicedData;

                            if (index.tabix) {

                                inflated = igvxhr.arrayBufferToString(unbgzf(data));
                                // need to decompress data
                            }
                            else {
                                inflated = data;
                            }
                           

                            slicedData = startOffset ? inflated.slice(startOffset) : inflated;
                            var f = self.parser.parseFeatures(slicedData,chr,start,end);
                            fulfill(f);
                        };


                        // Async load
                        if (self.localFile) {
                            igvxhr.loadStringFromFile(self.localFile, options).then(success).catch(reject);
                        }
                        else {
                            if (index.tabix) {
                            
                                igvxhr.loadArrayBuffer(self.url, options).then(success).catch(reject);
                            }
                            else {
                                igvxhr.loadString(self.url, options).then(success).catch(reject);
                            }
                        }
                    }))
                });

                Promise.all(promises).then(function (featureArrays) {

                    var i, allFeatures;

                    if (featureArrays.length === 1) {
                        allFeatures = featureArrays[0];
                    } else {
                        allFeatures = featureArrays[0];

                        for (i = 1; i < featureArrays.length; i++) {
                            allFeatures = allFeatures.concat(featureArrays[i]);
                        }

                        allFeatures.sort(function (a, b) {
                            return a.start - b.start;
                        });
                    }

                    fulfill(allFeatures)
                }).catch((e)=>{
                    console.log(e);
                    reject(e)}
                    );
            }
        });

    }
    


    getIndex() {

        var self = this,
        isIndeedIndexible = this.isIndexable();
        return new Promise(function (fulfill, reject) {

            if (self.indexed === undefined && isIndeedIndexible) {
                self.loadIndex().then(function (index) {
                    if (index) {
                        self.index = index;
                        self.indexed = true;
                    }
                    else {
                        self.indexed = false;
                    }
                    fulfill(self.index);
                }).catch(reject);
            }
            else {
                fulfill(self.index);   // Is either already loaded, or there isn't one
            }

        });
    }

    readHeader() {

        var self = this;

        return new Promise(function (fulfill, reject) {


            if (self.header) {
                fulfill(self.header);
            }

            else {

                // We force a load of the index first

               self. getIndex().then(function (index) {

                    if (index) {
                        // Load the file header (not HTTP header) for an indexed file.
                        // TODO -- note this will fail if the file header is > 65kb in size
                        var options = {
                                headers: self.config.headers,           // http headers, not file header
                                bgz: index.tabix,
                                range: {start: 0, size: 65000},
                                withCredentials: self.config.withCredentials
                            },
                            success = function (data) {
                                self.header = self.parser.parseHeader(data);
                                fulfill(self.header);
                            };

                        if (self.localFile) {
                            igvxhr.loadStringFromFile(self.localFile, options).then(success);
                        }
                        else {
                            const url = self.url; //getProjectURL(self.url);
                            igvxhr.loadString(url, options).then(success).catch(reject);
                        }
                    }
                    else {
                        self.loadFeaturesNoIndex(undefined).then(function (features) {
                            var header = self.header || {};
                            header.features = features;
                            fulfill(header);
                        }).catch(error);
                    }
                }).catch(reject);
            }
        });

    }

    /**
     *
     * @param chr
     * @param start 
     * @param end 
     */
    readFeatures(chr, start, end) {

        var self = this;

        return new Promise(function (fulfill, reject) {

            if (self.index) {
                self.loadFeaturesWithIndex(chr, start, end).then(packFeatures).catch(reject);
            }
            else {
                self.loadFeaturesNoIndex().then(packFeatures);
            }

            function packFeatures(features) {
                // TODO pack
                fulfill(features);
            }

        });
    }


}


//***js/feature/FeatureParsers.js**************


const maxFeatureCount = Number.MAX_VALUE;    // For future use,  controls downsampling

const gffNameFields = ["Name", "gene_name", "gene", "gene_id", "alias", "locus"];

    /**
     * A factory function.  Return a parser for the given file format.
     */
class FeatureParser{
    constructor(format, decode_func, config) {

        var customFormat;

        this.format = format;
        this.nameField = config ? config.nameField : undefined;
        this.skipRows = 0;   // The number of fixed header rows to skip.  Override for specific types as neede

        if (decode_func) {
            this.delimiter = /\s+/;

            this.decode =function(tokens,ignore){
                let feature={chr:tokens[0],start:parseInt(tokens[1]),end:parseInt(tokens[2])}
                decode_func(tokens.slice(3),feature);
                return feature;
            } ;
        }
        else if (config.decode_function==="generic"){
             this.decode =function(tokens,ignore){
                let feature={chr:tokens[0],start:parseInt(tokens[1]),end:parseInt(tokens[2])}
                feature.data=tokens.slice(3);
                return feature;
            } ;
        }
        else{
            this.decode =config.decode_function === "decodeRefflat"?FeatureParser.decodeRefflat:FeatureParser.decodeBed;
            this.delimiter = /\s+/;

        }

                                       


       /* switch (format) {
            case "narrowpeak":
            case "broadpeak":
            case "peaks":
                this.decode = this.decodePeak;
                this.delimiter = /\s+/;
                break;
            case "bedgraph":
                this.decode = this.decodeBedGraph;
                this.delimiter = /\s+/;
                break;
            case "wig":
                this.decode = this.decodeWig;
                this.delimiter = /\s+/;
                break;
            case "gff3" :
            case "gff" :
            case "gtf" :
                this.decode = this.decodeGFF;
                this.delimiter = "\t";
                break;
            case "aneu":
                this.decode = this.decodeAneu;
                this.delimiter = "\t";
                break;
            case "fusionjuncspan":
                // bhaas, needed for FusionInspector view
                this.decode = this.decodeFusionJuncSpan;
                this.delimiter = /\s+/;
                break;
            case "gtexgwas":
                this.skipRows = 1;
                this.decode = this.decodeGtexGWAS;
                this.delimiter = "\t";
                break;
            case "refflat":
                this.decode = this.decodeRefflat;
                this.delimiter = "\t";
                break;
            default:

               customFormat = igv.browser.getFormat(format);
                if (customFormat !== undefined) {
                    this.decode = decodeCustom;
                    this.format = customFormat;
                    this.delimiter = customFormat.delimiter || "\t";
                }

                else {
                
               // }

        }*/

    };

    parseHeader(data) {

        var lines = data.split("\n"),
            len = lines.length,
            line,
            i,
            header;

        for (i = 0; i < len; i++) {
            line = lines[i];
            if (line.startsWith("track") || line.startsWith("#") || line.startsWith("browser")) {
                if (line.startsWith("track")) {
                    header = this.parseTrackLine(line);
                }
                else if (line.startsWith("##gff-version 3")) {
                    this.format = "gff3";
                    if (!header) header = {};
                    header["format"] = "gff3";
                }
            }
            else {
                header={};
                break;
            }
        }
        return header;
    };

    parseFeatures(data,chr,start,end) {

        if (!data) return null;

        var wig,
            feature,
            lines = data.split("\n"),
            len = lines.length,
            tokens,
            allFeatures = [],
            line,
            i,
            cnt = 0,
            j,
            decode = this.decode,
            format = this.format,
           delimiter = this.delimiter || "\t";
       

        for (i = this.skipRows; i < len; i++) {
            line = lines[i];
            if (line.startsWith("track") || line.startsWith("#") || line.startsWith("browser")) {
                continue;
            }
            else if (format === "wig" && line.startsWith("fixedStep")) {
                wig = this.parseFixedStep(line);
                continue;
            }
            else if (format === "wig" && line.startsWith("variableStep")) {
                wig = this.parseVariableStep(line);
                continue;
            }

            tokens = lines[i].split(delimiter);
            if (tokens.length < 1 || (format!=="wig" && tokens.length<3)) continue;

            feature = this.decode(tokens, wig);


            if (feature) {
                if (feature.chr !==chr || feature.end<start || feature.start>end){
                    continue;
                }
                if (allFeatures.length < maxFeatureCount) {
                    allFeatures.push(feature);
                }
                else {
                    // Reservoir sampling,  conditionally replace existing feature with new one.
                    j = Math.floor(Math.random() * cnt);
                    if (j < maxFeatureCount) {
                        allFeatures[j] = feature;
                    }
                }
                cnt++;
            }   
        }

        return allFeatures;
    };


    static parseFixedStep(line) {

        var tokens = line.split(/\s+/),
            cc = tokens[1].split("=")[1],
            ss = parseInt(tokens[2].split("=")[1], 10),
            step = parseInt(tokens[3].split("=")[1], 10),
            span = (tokens.length > 4) ? parseInt(tokens[4].split("=")[1], 10) : 1;

        return {format: "fixedStep", chrom: cc, start: ss, step: step, span: span, index: 0};

    }

    static parseVariableStep(line) {

        var tokens = line.split(/\s+/),
            cc = tokens[1].split("=")[1],
            span = tokens.length > 2 ? parseInt(tokens[2].split("=")[1], 10) : 1;
        return {format: "variableStep", chrom: cc, span: span}

    }

    static parseTrackLine(line) {
        var properties = {},
            tokens = line.split(/(?:")([^"]+)(?:")|([^\s"]+)(?=\s+|$)/g),
            tmp = [],
            i, tk, curr;

        // Clean up tokens array
        for (i = 1; i < tokens.length; i++) {
            if (!tokens[i] || tokens[i].trim().length === 0) continue;

            tk = tokens[i].trim();

            if (tk.endsWith("=") > 0) {
                curr = tk;
            }
            else if (curr) {
                tmp.push(curr + tk);
                curr = undefined;
            }
            else {
                tmp.push(tk);
            }

        }


        tmp.forEach(function (str) {
            if (!str) return;
            var kv = str.split('=', 2);
            if (kv.length == 2) {
                properties[kv[0]] = kv[1];
            }

        });

        return properties;
    }

    /**
     * Decode the "standard" UCSC bed format
     * @param tokens
     * @param ignore
     * @returns decoded feature, or null if this is not a valid record
     */
    static decodeBed(tokens, ignore) {

        var chr, start, end, id, name, tmp, idName, exonCount, exonSizes, exonStarts, exons, exon, feature,
            eStart, eEnd;

        if (tokens.length < 3) return null;

        chr = tokens[0];
        start = parseInt(tokens[1]);
        end = tokens.length > 2 ? parseInt(tokens[2]) : start + 1;

        feature = {chr: chr, start: start, end: end, score: 1000};

        if (tokens.length > 3) {
            // Note: these are very special rules for the gencode gene files.
            tmp = tokens[3].replace(/"/g, '');
            idName = tmp.split(';');
            for (var i = 0; i < idName.length; i++) {
                var kv = idName[i].split('=');
                if (kv[0] == "gene_id") {
                    id = kv[1];
                }
                if (kv[0] == "gene_name") {
                    name = kv[1];
                }
            }
            feature.id = id ? id : tmp;
            feature.name = name ? name : tmp;
        }

        if (tokens.length > 4) {
            feature.score = parseFloat(tokens[4]);
        }
        if (tokens.length > 5) {
            feature.strand = tokens[5];
        }
        if (tokens.length > 6) {
            feature.cdStart = parseInt(tokens[6]);
        }
        if (tokens.length > 7) {
            feature.cdEnd = parseInt(tokens[7]);
        }
        if (tokens.length > 8) {
           // if (tokens[8] !== "." && tokens[8] !== "0")
                //feature.color = igv.createColorString(tokens[8]);
        }
        if (tokens.length > 11) {
            exonCount = parseInt(tokens[9]);
            exonSizes = tokens[10].split(',');
            exonStarts = tokens[11].split(',');
            exons = [];

            for (var i = 0; i < exonCount; i++) {
                eStart = start + parseInt(exonStarts[i]);
                eEnd = eStart + parseInt(exonSizes[i]);
                var exon = {start: eStart, end: eEnd};

                if (feature.cdStart > eEnd || feature.cdEnd < feature.cdStart) exon.utr = true;   // Entire exon is UTR
                if (feature.cdStart >= eStart && feature.cdStart <= eEnd) exon.cdStart = feature.cdStart;
                if (feature.cdEnd >= eStart && feature.cdEnd <= eEnd) exon.cdEnd = feature.cdEnd;

                exons.push(exon);
            }

            feature.exons = exons;
        }

       /* feature.popupData = function () {
            var data = [];
            if (feature.name) data.push({name: "Name", value: feature.name});
            if ("+" === feature.strand || "-" === feature.strand) data.push({name: "Strand", value: feature.strand});
            return data;
        };*/

        return feature;

    }

    /**
     * Decode a UCSC "refflat" record
     * @param tokens
     * @param ignore
     * @returns {*}
     */
    static decodeRefflat(tokens, ignore) {

        if (tokens.length < 10) return null;
        if (tokens[10] === undefined){
            console.log(tokens[1]);
            return null;
        }

        var feature = {
                chr: tokens[2],
                start: parseInt(tokens[4]),
                end: parseInt(tokens[5]),
                id: tokens[1],
                name: tokens[0],
                strand: tokens[3],
                cdStart: parseInt(tokens[6]),
                cdEnd: parseInt(tokens[7])
            },
            exonCount = parseInt(tokens[8]),
            exonStarts = tokens[9].split(','),
            exonEnds = tokens[10].split(','),
            exons = [];

        for (var i = 0; i < exonCount; i++) {
          
            let exon =  {start:parseInt(exonStarts[i]), end: parseInt(exonEnds[i])}
            if (feature.cdStart > exon.end || feature.cdEnd < feature.cdStart) exon.utr = true;   // Entire exon is UTR
            if (feature.cdStart >= exon.start && feature.cdStart <= exon.end) exon.cdStart = feature.cdStart;
            if (feature.cdEnd >= exon.start && feature.cdEnd <= exon.end) exon.cdEnd = feature.cdEnd;
            exons.push(exon);
        }

        feature.exons = exons;

        return feature;

    }

    static decodePeak(tokens, ignore) {

        var tokenCount, chr, start, end, strand, name, score, qValue, signal, pValue;

        tokenCount = tokens.length;
        if (tokenCount < 9) {
            return null;
        }

        chr = tokens[0];
        start = parseInt(tokens[1]);
        end = parseInt(tokens[2]);
        name = tokens[3];
        score = parseFloat(tokens[4]);
        strand = tokens[5].trim();
        signal = parseFloat(tokens[6]);
        pValue = parseFloat(tokens[7]);
        qValue = parseFloat(tokens[8]);

        if (score === 0) score = signal;

        return {
            chr: chr, start: start, end: end, name: name, score: score, strand: strand, signal: signal,
            pValue: pValue, qValue: qValue
        };
    }

    static decodeBedGraph(tokens, ignore) {

        var chr, start, end, value;

        if (tokens.length < 3) return null;

        chr = tokens[0];
        start = parseInt(tokens[1]);
        end = parseInt(tokens[2]);

        value = parseFloat(tokens[3]);

        return {chr: chr, start: start, end: end, value: value};
    }

    static decodeWig(tokens, wig) {

        var ss,
            ee,
            value;

        if (wig.format === "fixedStep") {

            ss = (wig.index * wig.step) + wig.start;
            ee = ss + wig.span;
            value = parseFloat(tokens[0]);
            ++(wig.index);
            return isNaN(value) ? null : {chr: wig.chrom, start: ss, end: ee, value: value};
        }
        else if (wig.format === "variableStep") {

            if (tokens.length < 2) return null;

            ss = parseInt(tokens[0], 10);
            ee = ss + wig.span;
            value = parseFloat(tokens[1]);
            return isNaN(value) ? null : {chr: wig.chrom, start: ss, end: ee, value: value};

        }
        else {
            return decodeBedGraph(tokens);
        }
    }

    static decodeAneu(tokens, ignore) {

        var chr, start, end, feature;


        if (tokens.length < 4) return null;

        chr = tokens[1];
        start = parseInt(tokens[2]);
        end = tokens.length > 3 ? parseInt(tokens[3]) : start + 1;

        feature = {chr: chr, start: start, end: end};

        if (tokens.length > 4) {
            feature.score = parseFloat(tokens[4]);
            feature.value = feature.score;
        }


        feature.popupData = function () {
            return [{name: "Name", value: feature.name}];
        };

        return feature;

    }

    static decodeFusionJuncSpan(tokens, ignore) {

        /*
         Format:

         0       #scaffold
         1       fusion_break_name
         2       break_left
         3       break_right
         4       num_junction_reads
         5       num_spanning_frags
         6       spanning_frag_coords

         0       B3GNT1--NPSR1
         1       B3GNT1--NPSR1|2203-10182
         2       2203
         3       10182
         4       189
         5       1138
         6       1860-13757,1798-13819,1391-18127,1443-17174,...

         */


       

        var chr = tokens[0];
        var fusion_name = tokens[1];
        var junction_left = parseInt(tokens[2]);
        var junction_right = parseInt(tokens[3]);
        var num_junction_reads = parseInt(tokens[4]);
        var num_spanning_frags = parseInt(tokens[5]);

        var spanning_frag_coords_text = tokens[6];

        var feature = {
            chr: chr,
            name: fusion_name,
            junction_left: junction_left,
            junction_right: junction_right,
            num_junction_reads: num_junction_reads,
            num_spanning_frags: num_spanning_frags,
            spanning_frag_coords: [],

            start: -1,
            end: -1
        }; // set start and end later based on min/max of span coords

        var min_coord = junction_left;
        var max_coord = junction_right;

        if (num_spanning_frags > 0) {

            var coord_pairs = spanning_frag_coords_text.split(',');

            for (var i = 0; i < coord_pairs.length; i++) {
                var split_coords = coord_pairs[i].split('-');

                var span_left = split_coords[0];
                var span_right = split_coords[1];

                if (span_left < min_coord) {
                    min_coord = span_left;
                }
                if (span_right > max_coord) {
                    max_coord = span_right;
                }
                feature.spanning_frag_coords.push({left: span_left, right: span_right});

            }
        }

        feature.start = min_coord;
        feature.end = max_coord;


        feature.popupData = function () {
            return [{name: "Name", value: feature.name}];
        };

        return feature;

    }

    static decodeGtexGWAS(tokens, ignore) {


        var tokenCount, chr, start, end, strand, name, score, qValue, signal, pValue;

        tokenCount = tokens.length;
        if (tokenCount < 8) {
            return null;
        }

        chr = tokens[0];
        start = parseInt(tokens[1]) - 1;
        end = parseInt(tokens[3].split(':')[1]);
        //name = tokens[3];
        //score = parseFloat(tokens[4]);
        //strand = tokens[5].trim();
        //signal = parseFloat(tokens[6]);
        pValue = parseFloat(tokens[5]);
        //qValue = parseFloat(tokens[8]);

        //return {chr: chr, start: start, end: end, name: name, score: score, strand: strand, signal: signal,
        //    pValue: pValue, qValue: qValue};
        return {chr: chr, start: start, end: end, pvalue: pValue};
    }

    /**
     * Decode a single gff record (1 line in file).  Aggregations such as gene models are constructed at a higher level.
     *      ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00001;Parent=gene00001
     * @param tokens
     * @param ignore
     * @returns {*}
     */
    static decodeGFF(tokens, ignore) {

        var tokenCount, chr, start, end, strand, type, score, phase, attributeString, id, parent, color, name,
            transcript_id, i,
            format = this.format;

        tokenCount = tokens.length;
        if (tokenCount < 9) {
            return null;      // Not a valid gff record
        }

        chr = tokens[0];
        type = tokens[2];
        start = parseInt(tokens[3]) - 1;
        end = parseInt(tokens[4]);
        score = "." === tokens[5] ? 0 : parseFloat(tokens[5]);
        strand = tokens[6];
        phase = "." === tokens[7] ? 0 : parseInt(tokens[7]);
        attributeString = tokens[8];

        // Find ID and Parent, or transcript_id
        var delim = ('gff3' === format) ? '=' : /\s+/;
        var attributes = {};
        attributeString.split(';').forEach(function (kv) {
            var t = kv.trim().split(delim, 2), key, value;
            if (t.length == 2) {
                key = t[0].trim();
                value = t[1].trim();
                //Strip off quotes, if any
                if (value.startsWith('"') && value.endsWith('"')) {
                    value = value.substr(1, value.length - 2);
                }
                if ("ID" === t[0]) id = t[1];
                else if ("Parent" === t[0]) parent = t[1];
                else if ("color" === t[0].toLowerCase()) color = igv.createColorString(t[1]);
                else if ("transcript_id" === t[0]) id = t[1];     // gtf format
                attributes[key] = value;
            }
        });

        // Find name (label) property
        if (this.nameField) {
            name = attributes[this.nameField];
        }
        else {
            for (i = 0; i < gffNameFields.length; i++) {
                if (attributes.hasOwnProperty(gffNameFields[i])) {
                    this.nameField = gffNameFields[i];
                    name = attributes[this.nameField];


                    break;
                }
            }
        }


        return {
            id: id,
            parent: parent,
            name: name,
            type: type,
            chr: chr,
            start: start,
            end: end,
            score: score,
            strand: strand,
            color: color,
            attributeString: attributeString,
            popupData: function () {
                var kvs = this.attributeString.split(';'),
                    pd = [],
                    key, value;
                kvs.forEach(function (kv) {
                    var t = kv.trim().split(delim, 2);
                    if (t.length === 2 && t[1] !== undefined) {
                        key = t[0].trim();
                        value = t[1].trim();
                        //Strip off quotes, if any
                        if (value.startsWith('"') && value.endsWith('"')) {
                            value = value.substr(1, value.length - 2);
                        }
                        pd.push({name: key, value: value});
                    }
                });
                return pd;
            }

        };
    }

    /**
     * Decode the "standard" UCSC bed format
     * @param tokens
     * @param ignore
     * @returns decoded feature, or null if this is not a valid record
     */
    decodeCustom(tokens, ignore) {

        var feature,
            chr, start, end,
            format = this.format,         // "this" refers to FeatureParser instance
            coords = format.coords || 0;

        if (tokens.length < 3) return null;

        chr = tokens[format.chr];
        start = parseInt(tokens[format.start]) - coords;
        end = format.end !== undefined ? parseInt(tokens[format.end]) : start + 1;

        feature = {chr: chr, start: start, end: end};

        if (format.fields) {
            format.fields.forEach(function (field, index) {
                if (index != format.chr && index != format.start && index != format.end) {
                    feature[field] = tokens[index];
                }
            });
        }

        return feature;

    }


}


//*******js/feature/featureCache.js**********************


    /**
     * Object for caching lists of features.  Supports effecient queries for sub-range  (chr, start, end)
     *
     * @param featureList
     * @param The genomic range spanned by featureList (optional)
     * @constructor
     */

class FeatureCache{
    constructor(featureList, range) {
        this.treeMap = FeatureCache.buildTreeMap(featureList);
        this.range = range;
    }

    queryFeatures(chr, start, end) {
         

        var featureList, intervalFeatures, feature, len, i, tree, intervals;

        tree = this.treeMap[chr];

        if (!tree) return [];

        intervals = tree.findOverlapping(start, end);

        if (intervals.length == 0) {
            return [];
        }
        else {
            // Trim the list of features in the intervals to those
            // overlapping the requested range.
            // Assumption: features are sorted by start position

            featureList = [];

            intervals.forEach(function (interval) {
                intervalFeatures = interval.value;
                len = intervalFeatures.length;
                for (i = 0; i < len; i++) {
                    feature = intervalFeatures[i];
                    if (feature.start > end) break;
                    else if (feature.end >= start) {
                        featureList.push(feature)
                    }
                }
            });
            return featureList;
        }

    };

    allFeatures() {

        var allFeatures = [];
        var treeMap = this.treeMap;
        if (treeMap) {
            for (var key in treeMap) {
                if (treeMap.hasOwnProperty(key)) {

                    var tree = treeMap[key];
                    tree.mapIntervals(function (interval) {
                        allFeatures = allFeatures.concat(interval.value);
                    });
                }
            }
        }
        return allFeatures;

    }

    static buildTreeMap(featureList) {

        var featureCache = {},
            chromosomes = [],
            treeMap = {},
            genome = null;

        if (featureList) {

            featureList.forEach(function (feature) {

                var chr = feature.chr,
                    geneList;

                // Translate to "official" name
                if(genome) chr = genome.getChromosomeName(chr);

                geneList = featureCache[chr];

                if (!geneList) {
                    chromosomes.push(chr);
                    geneList = [];
                    featureCache[chr] = geneList;
                }

                geneList.push(feature);

            });


            // Now build interval tree for each chromosome

            for (let i = 0; i < chromosomes.length; i++) {
                let chr = chromosomes[i];
                treeMap[chr] =FeatureCache.buildIntervalTree(featureCache[chr]);
            }
        }

        return treeMap;
    };

    /**
     * Build an interval tree from the feature list for fast interval based queries.   We lump features in groups
     * of 10, or total size / 100,   to reduce size of the tree.
     *
     * @param featureList
     */
    static buildIntervalTree(featureList) {

        var i, e, iStart, iEnd, tree, chunkSize, len, subArray;

        tree = new IntervalTree();
        len = featureList.length;

        chunkSize = Math.max(10, Math.round(len / 100));

        featureList.sort(function (f1, f2) {
            return (f1.start === f2.start ? 0 : (f1.start > f2.start ? 1 : -1));
        });

        for (i = 0; i < len; i += chunkSize) {
            e = Math.min(len, i + chunkSize);
            subArray = featureList.slice(i, e);
            iStart = subArray[0].start;
            //
            iEnd = iStart;
            subArray.forEach(function (feature) {
                iEnd = Math.max(iEnd, feature.end);
            });
            tree.insert(iStart, iEnd, subArray);
        }

        return tree;
    }


}

//*****js/intervalTree.js************


const BLACK = 1;
const RED = 2;

let NIL = {}
NIL.color = BLACK;
NIL.parent = NIL;
NIL.left = NIL;
NIL.right = NIL;

  
class IntervalTree{
    constructor() {    
        this.root = NIL;
    }


    insert(start, end, value) {

        var interval = new Interval(start, end, value);
        var x = new Node(interval);
        this.treeInsert(x);
        x.color = RED;
        while (x != this.root && x.parent.color == RED) {
            if (x.parent == x.parent.parent.left) {
                var y = x.parent.parent.right;
                if (y.color == RED) {
                    x.parent.color = BLACK;
                    y.color = BLACK;
                    x.parent.parent.color = RED;
                    x = x.parent.parent;
                } else {
                    if (x == x.parent.right) {
                        x = x.parent;
                        this.leftRotate(x);
                    }
                    x.parent.color = BLACK;
                    x.parent.parent.color = RED;
                    this.rightRotate(x.parent.parent);
                }
            } else {
                var y = x.parent.parent.left;
                if (y.color == RED) {
                    x.parent.color = BLACK;
                    y.color = BLACK;
                    x.parent.parent.color = RED;
                    x = x.parent.parent;
                } else {
                    if (x == x.parent.left) {
                        x = x.parent;
                        this.rightRotate(x);
                    }
                    x.parent.color = BLACK;
                    x.parent.parent.color = RED;
                    this.leftRotate(x.parent.parent);
                }
            }
        }
        this.root.color = BLACK;
    }


    /**
     *
     * @param start - query interval
     * @param end - query interval
     * @returns Array of all intervals overlapping the query region
     */
    findOverlapping(start, end) {


        var searchInterval = new Interval(start, end, 0);

        if (this.root === NIL) return [];

        var intervals = this.searchAll(searchInterval, this.root, []);

        if(intervals.length > 1) {
            intervals.sort(function(i1, i2) {
                 return i1.low - i2.low;
            });
        }

        return intervals;
    }

    /**
     * Dump info on intervals to console.  For debugging.
     */
    logIntervals() {

        logNode(this.root, 0);

        function logNode(node, indent) {

            var space = "";
            for(var i=0; i<indent; i++) space += " ";
            console.log(space + node.interval.low + " " + node.interval.high); // + " " + (node.interval.value ? node.interval.value : " null"));

            indent += 5;

            if(node.left != NIL) logNode(node.left, indent);
            if(node.right != NIL) logNode(node.right, indent);
        }

    }


    mapIntervals(func) {

        applyInterval(this.root);

        function applyInterval(node) {

            func(node.interval);

            if(node.left != NIL) applyInterval(node.left);
            if(node.right != NIL) applyInterval(node.right);
        }
    }

    searchAll(interval, node, results) {

        if (node.interval.overlaps(interval)) {
            results.push(node.interval);
        }

        if (node.left != NIL && node.left.max >= interval.low) {
            this.searchAll(interval, node.left, results);
        }

        if (node.right != NIL && node.right.min <= interval.high) {
            this.searchAll(interval, node.right, results);
        }

        return results;
    }

    leftRotate(x) {
        var y = x.right;
        x.right = y.left;
        if (y.left != NIL) {
            y.left.parent = x;
        }
        y.parent = x.parent;
        if (x.parent == NIL) {
            this.root = y;
        } else {
            if (x.parent.left == x) {
                x.parent.left = y;
            } else {
                x.parent.right = y;
            }
        }
        y.left = x;
        x.parent = y;

        this.applyUpdate(x);
        // no need to apply update on y, since it'll y is an ancestor
        // of x, and will be touched by applyUpdate().
    }


    rightRotate(x) {
        var y = x.left;
        x.left = y.right;
        if (y.right != NIL) {
            y.right.parent = x;
        }
        y.parent = x.parent;
        if (x.parent == NIL) {
            this.root = y;
        } else {
            if (x.parent.right == x) {
                x.parent.right = y;
            } else {
                x.parent.left = y;
            }
        }
        y.right = x;
        x.parent = y;


        this.applyUpdate(x);
        // no need to apply update on y, since it'll y is an ancestor
        // of x, and will be touched by applyUpdate().
    }


    /**
     * Note:  Does not maintain RB constraints,  this is done post insert
     *
     * @param x  a Node
     */
   treeInsert(x) {
        var node = this.root;
        var y = NIL;
        while (node != NIL) {
            y = node;
            if (x.interval.low <= node.interval.low) {
                node = node.left;
            } else {
                node = node.right;
            }
        }
        x.parent = y;

        if (y == NIL) {
            this.root = x;
            x.left = x.right = NIL;
        } else {
            if (x.interval.low <= y.interval.low) {
                y.left = x;
            } else {
                y.right = x;
            }
        }

        this.applyUpdate(x);
    }


    // Applies the statistic update on the node and its ancestors.
    applyUpdate (node) {
        while (node != NIL) {
            var nodeMax = node.left.max > node.right.max ? node.left.max : node.right.max;
            var intervalHigh = node.interval.high;
            node.max = nodeMax > intervalHigh ? nodeMax : intervalHigh;

            var nodeMin = node.left.min < node.right.min ? node.left.min : node.right.min;
            var intervalLow = node.interval.low;
            node.min = nodeMin < intervalLow ? nodeMin : intervalLow;

            node = node.parent;
        }
    }

}


class Interval {
    constructor(low, high, value) {
        this.low = low;
        this.high = high;
        this.value = value;
    }


    equals(other) {
        if (!other) {
            return false;
        }
        if (this == other) {
            return true;
        }
        return (this.low == otherInterval.low &&
            this.high == otherInterval.high);

    }


    compareTo(other) {
        if (this.low < other.low)
            return -1;
        if (this.low > other.low)
            return 1;

        if (this.high < other.high)
            return -1;
        if (this.high > other.high)
            return 1;

        return 0;
    }

    /**
     * Returns true if this interval overlaps the other.
     */
    overlaps(other) {
        try {
            return (this.low <= other.high && other.low <= this.high);
        } catch (e) {
            //alert(e);
            igv.presentAlert(e);
        }
    }
}

class Node{
    constructor(interval) {
        this.parent = NIL;
        this.left = NIL;
        this.right = NIL;
        this.interval = interval;
        this.color = RED;
    }
}

class GenomicInterval{

	constructor(chr, start, end, features) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.features = features;
    }

    contains (chr, start, end) {
        return this.chr == chr &&
            this.start <= start &&
            this.end >= end;
    }

    containsRange(range) {
        return this.chr === range.chr &&
            this.start <= range.start &&
            this.end >= range.end;
    }
    rangesToGet(range){
        let needs_range=false;
        let ranges={};
        if (this.chr !== range.chr){
            ranges.all=[range.start,range.end];
            needs_range=true;
            this.start=range.start;
            this.end=range.end;

        }
        else{   
            if (range.start<this.start){
                ranges.left=[range.start,this.start];
                needs_range=true;
                this.start=range.start;
              
            }
            if (range.end>this.end){
                ranges.right=[this.end,range.end];
                needs_range=true;
                this.end=range.end;
            }
        }
        if (!needs_range){
            return false;
        }
        return ranges;
    }
}



class FastaSequence{

    constructor(url) {

        this.file = url;
        this.indexed = true;
        if (this.indexed) {
            this.indexFile = this.file + ".fai";
        }
    

    }

   init(){

        var self = this;

        if (self.indexed) {

            return new Promise(function (fulfill, reject) {

                self.getIndex().then(function (index) {
                    var order = 0;
                    self.chromosomes = {};
                    self.chromosomeNames.forEach(function (chrName) {
                        var bpLength = self.index[chrName].size;
                        self.chromosomes[chrName] = new igv.Chromosome(chrName, order++, bpLength);
                    });


                    // Ignore index, getting chr names as a side effect.  Really bad practice
                    fulfill();
                }).catch(reject);
            });
        }
        else {
            return self.loadAll();
        }

    }

    getSequence(chr, start, end) {

        if (this.indexed) {
            return this.getSequenceIndexed(chr, start, end);
        }
        else {
            return getSequenceNonIndexed.this(chr, start, end);

        }

    }

    getSequenceIndexed(chr, start, end) {

        var self = this;

        return new Promise(function (fulfill, reject) {
            var interval = self.interval;

            if (interval && interval.contains(chr, start, end)) {

                fulfill(getSequenceFromInterval(interval, start, end));
            }
            else {

                //console.log("Cache miss: " + (interval === undefined ? "nil" : interval.chr + ":" + interval.start + "-" + interval.end));

                // Expand query, to minimum of 100kb
                var qstart = start;
                var qend = end;
                if ((end - start) < 100000) {
                    var w = (end - start);
                    var center = Math.round(start + w / 2);
                    qstart = Math.max(0, center - 50000);
                    qend = center + 50000;
                }


                self.readSequence(chr, qstart, qend).then(function (seqBytes) {
                    self.interval = new GenomicInterval(chr, qstart, qend, seqBytes);
                    fulfill(getSequenceFromInterval(self.interval, start, end));
                }).catch(reject);
            }

            function getSequenceFromInterval(interval, start, end) {
                var offset = start - interval.start;
                var n = end - start;
                var seq = interval.features ? interval.features.substr(offset, n) : null;
                return seq;
            }
        });
    }


    getSequenceNonIndexed(chr, start, end) {

        var self = this;

        return new Promise(function (fulfill, reject) {
            var seq = self.sequences[chr];
            if (seq && seq.length > end) {
                fulfill(seq.substring(start, end));
            }
        });

    }

    getIndex() {

        var self = this;

        return new Promise(function (fulfill, reject) {

            if (self.index) {
                fulfill(self.index);
            } else {
                igvxhr.load(self.indexFile,{}).then(function (data) {
                    var lines = data.split("\n");
                    var len = lines.length;
                    var lineNo = 0;

                    self.chromosomeNames = [];     // TODO -- eliminate this side effect !!!!
                    self.index = {};               // TODO -- ditto
                    while (lineNo < len) {

                        var tokens = lines[lineNo++].split("\t");
                        var nTokens = tokens.length;
                        if (nTokens == 5) {
                            // Parse the index line.
                            var chr = tokens[0];
                            var size = parseInt(tokens[1]);
                            var position = parseInt(tokens[2]);
                            var basesPerLine = parseInt(tokens[3]);
                            var bytesPerLine = parseInt(tokens[4]);

                            var indexEntry = {
                                size: size, position: position, basesPerLine: basesPerLine, bytesPerLine: bytesPerLine
                            };

                            self.chromosomeNames.push(chr);
                            self.index[chr] = indexEntry;
                        }
                    }

                    if (fulfill) {
                        fulfill(self.index);
                    }
                }).catch(reject);
            }
        });
    }

    loadAll(){

        var self = this;

        return new Promise(function (fulfill, reject) {
            self.chromosomeNames = [];
            self.chromosomes = {};
            self.sequences = {};

            igvxhr.load(self.file, {
                withCredentials: self.withCredentials

            }).then(function (data) {

                var lines = data.splitLines(),
                    len = lines.length,
                    lineNo = 0,
                    nextLine,
                    currentSeq = "",
                    currentChr,
                    order = 0;


                while (lineNo < len) {
                    nextLine = lines[lineNo++].trim();
                    if (nextLine.startsWith("#") || nextLine.length === 0) {
                        continue;
                    }
                    else if (nextLine.startsWith(">")) {
                        if (currentSeq) {
                            self.chromosomeNames.push(currentChr);
                            self.sequences[currentChr] = currentSeq;
                            self.chromosomes[currentChr] = new igv.Chromosome(currentChr, order++, currentSeq.length);
                        }
                        currentChr = nextLine.substr(1).split("\\s+")[0];
                        currentSeq = "";
                    }
                    else {
                        currentSeq += nextLine;
                    }
                }

                fulfill();

            });
        });
    }

    readSequence(chr, qstart, qend) {

        //console.log("Read sequence " + chr + ":" + qstart + "-" + qend);
        var self = this;

        return new Promise(function (fulfill, reject) {
            self.getIndex().then(function () {

                var idxEntry = self.index[chr];
                if (!idxEntry) {
                    console.log("No index entry for chr: " + chr);

                    // Tag interval with null so we don't try again
                    self.interval = new GenomicInterval(chr, qstart, qend, null);
                    fulfill(null);

                } else {

                    var start = Math.max(0, qstart);    // qstart should never be < 0
                    var end = Math.min(idxEntry.size, qend);
                    var bytesPerLine = idxEntry.bytesPerLine;
                    var basesPerLine = idxEntry.basesPerLine;
                    var position = idxEntry.position;
                    var nEndBytes = bytesPerLine - basesPerLine;

                    var startLine = Math.floor(start / basesPerLine);
                    var endLine = Math.floor(end / basesPerLine);

                    var base0 = startLine * basesPerLine;   // Base at beginning of start line

                    var offset = start - base0;

                    var startByte = position + startLine * bytesPerLine + offset;

                    var base1 = endLine * basesPerLine;
                    var offset1 = end - base1;
                    var endByte = position + endLine * bytesPerLine + offset1 - 1;
                    var byteCount = endByte - startByte + 1;
                    if (byteCount <= 0) {
                        fulfill(null);
                    }

                    igvxhr.load(self.file, {
                        range: {start: startByte, size: byteCount}
                    }).then(function (allBytes) {

                        var nBases,
                            seqBytes = "",
                            srcPos = 0,
                            desPos = 0,
                            allBytesLength = allBytes.length;

                        if (offset > 0) {
                            nBases = Math.min(end - start, basesPerLine - offset);
                            seqBytes += allBytes.substr(srcPos, nBases);
                            srcPos += (nBases + nEndBytes);
                            desPos += nBases;
                        }

                        while (srcPos < allBytesLength) {
                            nBases = Math.min(basesPerLine, allBytesLength - srcPos);
                            seqBytes += allBytes.substr(srcPos, nBases);
                            srcPos += (nBases + nEndBytes);
                            desPos += nBases;
                        }

                        fulfill(seqBytes);
                    }).catch(reject)
                }
            }).catch(reject)
        });
    }
}


export {FeatureSource,FastaSequence,BigBedFeatureSource,TabixBedFeatureSource, FeatureFileReader}