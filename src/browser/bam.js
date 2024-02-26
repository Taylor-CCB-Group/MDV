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
import {BinaryParser} from "./bigwig.js";
import {jszlib_inflate_buffer,arrayCopy} from "./vendor/inflate.js";
import {FastaSequence} from "./feature.js";
import{Zlib} from "./vendor/zlib_and_gzip.js";


const BAI_MAGIC = 21578050;
const TABIX_MAGIC = 21578324;
const MAX_HEADER_SIZE = 100000000;   // IF the header is larger than this we can't read it !
const B_MAX_GZIP_BLOCK_SIZE = (1 << 16);


    /**
     * @param indexURL
     * @param config
     * @param tabix
     *
     * @returns a Promised for the bam or tabix index.  The fulfill function takes the index as an argument.
     */
let loadBamIndex = function (indexURL, config, tabix) {

        return new Promise(function (fulfill, reject) {

            var genome = null;

            igvxhr.loadArrayBuffer(indexURL,
                {
                    headers: config.headers,
                    withCredentials: config.withCredentials
                }).then(function (arrayBuffer) {

                var indices = [],
                    magic, nbin, nintv, nref, parser,
                    blockMin = Number.MAX_VALUE,
                    blockMax = 0,
                    binIndex, linearIndex, binNumber, cs, ce, b, i, ref, sequenceIndexMap;

                if (!arrayBuffer) {
                    fulfill(null);
                    return;
                }

                if (tabix) {
                    var inflate = new Zlib.Gunzip(new Uint8Array(arrayBuffer));
                    arrayBuffer = inflate.decompress().buffer;               
                }

                parser = new BinaryParser(new DataView(arrayBuffer));

                magic = parser.getInt();

                if (magic === BAI_MAGIC || (tabix && magic === TABIX_MAGIC)) {

                    nref = parser.getInt();


                    if (tabix) {
                        // Tabix header parameters aren't used, but they must be read to advance the pointer
                        var format = parser.getInt();
                        var col_seq = parser.getInt();
                        var col_beg = parser.getInt();
                        var col_end = parser.getInt();
                        var meta = parser.getInt();
                        var skip = parser.getInt();
                        var l_nm = parser.getInt();

                        sequenceIndexMap = {};
                        for (i = 0; i < nref; i++) {
                            var seq_name = parser.getString();

                            // Translate to "official" chr name.
                            if (genome) seq_name = genome.getChromosomeName(seq_name);

                            sequenceIndexMap[seq_name] = i;
                        }
                    }

                    for (ref = 0; ref < nref; ++ref) {

                        binIndex = {};
                        linearIndex = [];

                        nbin = parser.getInt();

                        for (b = 0; b < nbin; ++b) {

                            binNumber = parser.getInt();

                            if (binNumber == 37450) {
                                // This is a psuedo bin, not used but we have to consume the bytes
                                nchnk = parser.getInt(); // # of chunks for this bin
                                cs = parser.getVPointer();   // unmapped beg
                                ce = parser.getVPointer();   // unmapped end
                                var n_maped = parser.getLong();
                                var nUnmapped = parser.getLong();

                            }
                            else {
                                
                                binIndex[binNumber] = [];
                                var nchnk = parser.getInt(); // # of chunks for this bin

                                for (i = 0; i < nchnk; i++) {
                                    cs = parser.getVPointer();
                                    ce = parser.getVPointer();
                                    if (cs && ce) {
                                        if (cs.block < blockMin) {
                                            blockMin = cs.block;    // Block containing first alignment
                                        }
                                        if (ce.block > blockMax) {
                                            blockMax = ce.block;
                                        }
                                        binIndex[binNumber].push([cs, ce]);
                                    }
                                }
                            }
                        }


                        nintv = parser.getInt();
                        for (i = 0; i < nintv; i++) {
                            cs = parser.getVPointer();
                            linearIndex.push(cs);   // Might be null
                        }

                        if (nbin > 0) {
                            indices[ref] = {
                                binIndex: binIndex,
                                linearIndex: linearIndex
                            }
                        }
                    }

                } else {
                    throw new Error(indexURL + " is not a " + (tabix ? "tabix" : "bai") + " file");
                }
                fulfill(new BamIndex(indices, blockMin, blockMax, sequenceIndexMap, tabix));
            }).catch(reject);
        })
    }


class BamIndex{
    constructor (indices, blockMin,blockMax, sequenceIndexMap, tabix) {
        this.firstAlignmentBlock = blockMin;
        this.indices = indices;
        this.lastAlignmentBlock = blockMax;
        this.sequenceIndexMap = sequenceIndexMap;
        this.tabix = tabix;
    }

    //this seems to work just as well as  blocksForRange
    //but has not been extensively tested
    getBlocksForRange(refId,min,max){
        const index = this.indices[refId];
        if (!index){
            return [];
        }
        const lIndex =index.linearIndex;
        const nintv = lIndex.length;
        const minLin = Math.min(min >> 14, nintv - 1);
        const  maxLin = Math.min(max >> 14, nintv - 1);
        return [
            {
                minv:lIndex[minLin],
                maxv:lIndex[maxLin],
                bin:1
            }
        ];
    }

    /**
     * Fetch blocks for a particular genomic range.  This method is public so it can be unit-tested.
     *
     * @param refId  the sequence dictionary index of the chromosome
     * @param min  genomic start position
     * @param max  genomic end position
     * @param return an array of {minv: {filePointer, offset}, {maxv: {filePointer, offset}}
     */
    blocksForRange(refId, min, max) {

        var bam = this,
            ba = bam.indices[refId],
            overlappingBins,
            chunks,
            nintv,
            lowest,
            minLin,
            maxLin,
            vp,
            i

        if (!ba) {
            return [];
        }
        else {
            overlappingBins = BamIndex.reg2bins(min, max);
            chunks=[]        // List of bin #s that might overlap min, max
           


            overlappingBins.forEach(function (bin) {

                if (ba.binIndex[bin]) {
                    var binChunks = ba.binIndex[bin],
                        nchnk = binChunks.length;

                    for (var c = 0; c < nchnk; ++c) {
                        var cs = binChunks[c][0];
                        var ce = binChunks[c][1];
                        chunks.push({minv: cs, maxv: ce, bin: bin});
                    }

                }
            });

            // Use the linear index to find the lowest chunk that could contain alignments in the region
            nintv = ba.linearIndex.length;
            lowest = null;
            minLin = Math.min(min >> 14, nintv - 1), maxLin = Math.min(max >> 14, nintv - 1);
            for (i = minLin; i <= maxLin; ++i) {
                vp = ba.linearIndex[i];
                if (vp) {
                
                    if (!lowest || vp.isLessThan(lowest)) {
                        lowest = vp;
                    }
                }
            }
            // Prune chunks that end before the lowest chunk
            return optimizeChunks(chunks, lowest); 
        }
    };


    /**
     * Calculate the list of bins that may overlap with region [beg, end]
     *
     */
    static reg2bins(beg, end) {
        var i = 0, k, list = [];
        if (end >= 1 << 29) {
            end = 1 << 29;
        }
        --end;
        list.push(0);
        for (k = 1 + (beg >> 26); k <= 1 + (end >> 26); ++k) list.push(k);
        for (k = 9 + (beg >> 23); k <= 9 + (end >> 23); ++k) list.push(k);
        for (k = 73 + (beg >> 20); k <= 73 + (end >> 20); ++k) list.push(k);
        for (k = 585 + (beg >> 17); k <= 585 + (end >> 17); ++k) list.push(k);
        for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k) list.push(k);
        return list;
    }


}

function optimizeChunks(chunks, lowest) {

    var mergedChunks = [],
        lastChunk = null;

    if (chunks.length === 0) return chunks;

    chunks.sort(function (c0, c1) {
        var dif = c0.minv.block - c1.minv.block;
        if (dif != 0) {
            return dif;
        } else {
            return c0.minv.offset - c1.minv.offset;
        }
    });

    chunks.forEach(function (chunk) {

        if (!lowest || chunk.maxv.isGreaterThan(lowest)) {
            if (lastChunk === null) {
                mergedChunks.push(chunk);
                lastChunk = chunk;
            } else {
                if ((chunk.minv.block - lastChunk.maxv.block) < 65000) { // Merge chunks that are withing 65k of each other
                    if (chunk.maxv.isGreaterThan(lastChunk.maxv)) {
                        lastChunk.maxv = chunk.maxv;
                    }
                } else {
                    mergedChunks.push(chunk);
                    lastChunk = chunk;
                }
            }
        }
    });

    return mergedChunks;
}




class BGZFile{

    constructor (config) {
        this.filePosition = 0;
        this.config = config;
    }

    nextBlock() {

        var self = this;

        return new Promise(function (fulfill, reject) {

            igvxhr.loadArrayBuffer(self.path,
                {
                    headers: self.config.headers,
                    range: {start: self.filePosition, size: BLOCK_HEADER_LENGTH},
                    withCredentials: self.config.withCredentials

                }).then(function (arrayBuffer) {

                var ba = new Uint8Array(arrayBuffer);
                var xlen = (ba[11] << 8) | (ba[10]);
                var si1 = ba[12];
                var si2 = ba[13];
                var slen = (ba[15] << 8) | (ba[14]);
                var bsize = (ba[17] << 8) | (ba[16]) + 1;

                self.filePosition += BLOCK_HEADER_LENGTH;

                igvxhr.loadArrayBuffer(self.path, {
                    headers: self.config.headers,
                    range: {start: self.filePosition, size: bsize},
                    withCredentials: self.config.withCredentials

                }).then(function (arrayBuffer) {

                    var unc = jszlib_inflate_buffer(arrayBuffer);

                    self.filePosition += (bsize + 8);  // "8" for CRC-32 and size of uncompressed data

                    fulfill(unc);

                }).catch(reject)
            }).catch(reject);
        })

    }
}



    var BAM_MAGIC = 21840194;
    var SECRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
    var CIGAR_DECODER = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];
    var READ_STRAND_FLAG = 0x10;
    var MATE_STRAND_FLAG = 0x20;
    var FIRST_OF_PAIR_FLAG = 0x40;
    var SECOND_OF_PAIR_FLAG = 0x80;
    var READ_UNMAPPED_FLAG = 0x4;
    var MATE_UNMAPPED_FLAG = 0x8;
    var READ_PAIRED_FLAG = 0x1;
    var PROPER_PAIR_FLAG = 0x2;
    var SECONDARY_ALIGNMNET_FLAG = 0x100;
    var SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;




    var NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    var READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    var DUPLICATE_READ_FLAG = 0x400;
    var SUPPLEMENTARY_FLAG = 0x800;

    const MAX_GZIP_BLOCK_SIZE = 65536;   //  APPARENTLY.  Where is this documented???
    const DEFAULT_SAMPLING_WINDOW_SIZE = 100;
    const DEFAULT_SAMPLING_DEPTH = 2500;//50;
    const MAXIMUM_SAMPLING_DEPTH = 2500;



class BamFeatureReader{
    constructor(config){
        this.bamPath= config.url;
        this.baiPath = config.indexURL || this.bamPath + ".bai"; // If there is an indexURL provided, use it!
        this.headPath = config.headURL || this.bamPath;
        this.alignments=[];
        this.features=[];
        this.alignmentBuffer=new SharedArrayBuffer(20000000);
        this.alignments=new Int32Array(this.alignmentBuffer);
        this.pointer=0;
        this.config=config;
        this.fragmentThreshold=50;
    }

    async getAlignments(chr,bpStart,bpEnd){
        if (this.store){
            if (chr === this.store[0] && bpStart >= this.store[1] && bpEnd <= this.store[2]){
                return this.features;
            }
        }
        
        if (!this.index){
            await this.readHeader();
        }
        this.pointer=0;
        const chrId = this.chrToIndex[chr];
        const chunks = this.index.blocksForRange(chrId, bpStart, bpEnd);

        if (!chunks) {
            return [];       
        }
        if (chunks.length === 0) {     
            return [];
        }
        

        await Promise.all(chunks.map(async (c)=>{
            var fetchMin = c.minv.block,
            fetchMax = c.maxv.block + MAX_GZIP_BLOCK_SIZE,   // Make sure we get the whole block.
            range = {start: fetchMin, size: fetchMax - fetchMin + 1};

            const compressed= await igvxhr.loadArrayBuffer(this.bamPath,
            {
                headers: this.config.headers,
                range: range,
                withCredentials: this.config.withCredentials
            })
            var ba = new Uint8Array(new unbgzf(compressed));
            this.addAlignments(ba, c.minv.offset,bpStart, bpEnd, chrId);

        }));
        this.store= [chr,bpStart,bpEnd];
        //this.calculateFragmentSize();
        this.calculateCategories();
        
       
        return this.features;
      
    }

    setCategories(data,names){
        this.catData=data;
        this.catNames=names;
    }

    setFilter(arr,arr2){
        this.filter=arr;
        this.localFilter=arr2;
    }

    calculateFragmentSize(){
        const bpStart=this.store[1];
        const al = this.store[2]-this.store[1];
        this.features=[];
        for (let i=0;i<3;i++){
            this.features.push(new Uint16Array(al))
        }
        for (let i=0;i<this.pointer;i+=4){
            let st =  this.alignments[i]-bpStart;
            const en = st+ this.alignments[i+1];
            if (st<0){
                st=0;
            }
            let arr= this.features[0]

            if (this.alignments[i+3]>this.fragmentThreshold){
                arr=this.features[1];
            }
          
            if (en>st){
                for (let n = st;n<en;n++){
                    arr[n]++;
                    this.features[2][n]++;
                }

            }
            else{
                for (let n = en;n<st;n++){
                    arr[n]++;
                    this.features[2][n]++;
                }
            }
          
        }
        this.features.store=this.store;

        
    }


    calculateCategories(){
        const len =this.catNames.length;
        const bpStart=this.store[1];
        const al = this.store[2]-this.store[1];
        this.features=[];
        for (let i=0;i<len;i++){
            this.features.push(new Uint16Array(al))
        }
        this.cut_window=0;
        const t= performance.now();
        if (this.cut_window){
            const win= this.cut_window;
            for (let i=0;i<this.pointer;i+=4){
              
                let init =  this.alignments[i]-bpStart;
                let st = init-win;
                const en = init+win;
                if (st<0){
                    st=0;
                }
                const id = this.alignments[i+2];
                if (this.filter && this.filter[id]!==0 && (this.filter[id] !==this.localFilter[id])){
                   continue;
                }
                const arr = this.features[this.catData[id]];
                for (let n = st;n<en;n++){
                    arr[n]++;
                }
            
        }

        }
        else{
       
            for (let i=0;i<this.pointer;i+=4){
                if (this.size_filter_on && this.alignments[i+3]>this.size_filter){
                    continue;
                }
                let st =  this.alignments[i]-bpStart;
                const en = st+ this.alignments[i+1];
                if (st<0){
                    st=0;
                }
                const id = this.alignments[i+2];
                if (this.filter && this.filter[id]!==0 && (this.filter[id] !==this.localFilter[id])){
                continue;
                }
                const arr = this.features[this.catData[id]];
                if (en>st){
                    for (let n = st;n<en;n++){
                        arr[n]++;
                    }

                }
                else{
                    for (let n = en;n<st;n++){
                        arr[n]++;
                    }
                }
              
        }
        
    }
    this.features.store=this.store;

   
   


    }

    getIdsInRange(range){
        const st = Math.round(range.start);
        const en = Math.round(range.end);
        const a = this.alignments
        const iset = new Set();
        for (let i=0;i<this.alignments.length;i+=3){
            if (a[i]<st){
                continue
            }
            if(a[i]>en){
                break;
            }
            iset.add(a[i+2])

        }
        return iset;
        
    }



    //alignments is an int32 array with 
    // 1-start position 2-end position 3-id 4-length of PE fragment
    addAlignments(ba,offset,min,max,chrId){
    
        var blockSize,
            blockEnd,
            tags,
            refID,
            pos,
            bmn,
            nl,
            flag_nc,
            c,
            nc,
            lseq,
            lengthOnRef,
            p,
            seqBytes;
        let unknown=0;
        let num=0;
        while (true) {
            num++
            blockSize = readInt(ba, offset);
            blockEnd = offset + blockSize + 4;

            if (blockEnd > ba.length) {
                return ;
            }

            
            refID = readInt(ba, offset + 4);
            pos = readInt(ba, offset + 8);
           
            if(refID < 0) {
                return ;   // unmapped reads
            }
            else if (refID > chrId || pos > max) {
                return ;    // off right edge, we're done
            }
            else if (refID < chrId) {
                continue;   // to left of start, not sure this is possible
            }
           
           
            bmn = readInt(ba, offset + 12);
         
            nl = bmn & 0xff;

            flag_nc = readInt(ba, offset + 16);
            const flag = (flag_nc & 0xffff0000) >> 16;
            nc = flag_nc & 0xffff;

          
            const strand = !(flag & READ_STRAND_FLAG);

            lseq = readInt(ba, offset + 20);

            const  fragmentLength = Math.abs(readInt(ba, offset + 32));
          

            p = offset + 36 + nl;

            lengthOnRef = 0;
          
            for (c = 0; c < nc; ++c) {
                var cigop = readInt(ba, p);
                var opLen = (cigop >> 4);
                var opLtr = CIGAR_DECODER[cigop & 0xf];
                if (opLtr == 'M' || opLtr == 'EQ' || opLtr == 'X' || opLtr == 'D' || opLtr == 'N' || opLtr == '=')
                    lengthOnRef += opLen;
                
                p += 4;

                
            }
           
            if (pos + lengthOnRef < min) {
                offset = blockEnd;
                continue;  // Record out-of-range "to the left", skip to next one

            }
            if (!strand){
                pos = pos+lengthOnRef;
                lengthOnRef=-lengthOnRef
            }
            this.alignments[this.pointer++]=pos;
            this.alignments[this.pointer++] = lengthOnRef;
            

            
            
            
            seqBytes = (lseq + 1) >> 1;
            
            p += seqBytes;
           
            p += lseq;


           
            tags = decodeTags(new Uint8Array(ba.buffer.slice(p, blockEnd)));
            if (this.catIndex[tags.CB]===undefined){
                unknown++;
            }
            this.alignments[this.pointer++]=this.catIndex[tags.CB];  // decode these on demand
            this.alignments[this.pointer++] = Math.abs(fragmentLength);

            p += blockEnd;

          
            offset = blockEnd;
          
        }
     

    }

 

    async readHeader() {

        var self = this;
   
        const index = await getIndex(this);

        var len = index.firstAlignmentBlock + MAX_GZIP_BLOCK_SIZE;   // Insure we get the complete compressed block containing the header

        const compressedBuffer = await igvxhr.loadArrayBuffer(self.bamPath,
                    {
                        headers: self.config.headers,

                        range: {start: 0, size: len},

                        withCredentials: self.config.withCredentials
        });

        var unc = new unbgzf(compressedBuffer, len),
            uncba = new Uint8Array(unc),
            magic = readInt(uncba, 0),
            samHeaderLen = readInt(uncba, 4),
            samHeader = '',
            genome = null;

        for (var i = 0; i < samHeaderLen; ++i) {
            samHeader += String.fromCharCode(uncba[i + 8]);
        }

        var nRef = readInt(uncba, samHeaderLen + 8);
        var p = samHeaderLen + 12;

        self.chrToIndex = {};
        self.indexToChr = [];
        for (var i = 0; i < nRef; ++i) {
            var lName = readInt(uncba, p);
            var name = '';
            for (var j = 0; j < lName - 1; ++j) {
                name += String.fromCharCode(uncba[p + 4 + j]);
            }
            var lRef = readInt(uncba, p + lName + 4);
            //dlog(name + ': ' + lRef);

            if (genome && genome.getChromosomeName) {
                name = genome.getChromosomeName(name);
            }

            self.chrToIndex[name] = i;
            self.indexToChr.push(name);

            p = p + 8 + lName;
        }

                  
    }
}
/**
* Class for reading a bam file
*
* @param config
* @constructor
*/
class BamReader{
     constructor(config,parent) {

        this.config = config;
        this.parent=parent;
        this.ac_class=AlignmentContainer
        if (parent.ac_class){
            this.ac_class= parent_class;
        }

        this.filter = config.filter || new BamFilter();

        this.bamPath = config.url;
        // Todo - deal with Picard convention.  WHY DOES THERE HAVE TO BE 2?
        this.baiPath = config.indexURL || this.bamPath + ".bai"; // If there is an indexURL provided, use it!
        this.headPath = config.headURL || this.bamPath;


        this.samplingWindowSize = config.samplingWindowSize === undefined ? DEFAULT_SAMPLING_WINDOW_SIZE : config.samplingWindowSize;
        this.samplingDepth = config.samplingDepth === undefined ? DEFAULT_SAMPLING_DEPTH : config.samplingDepth;
        if(this.samplingDepth > MAXIMUM_SAMPLING_DEPTH) {
            igv.log("Warning: attempt to set sampling depth > maximum value of 2500");
            this.samplingDepth = MAXIMUM_SAMPLING_DEPTH;
        }

        if (config.viewAsPairs) {
            this.pairsSupported = true;
        }
        else {
            this.pairsSupported = config.pairsSupported === undefined ? true : config.pairsSupported;
        }

    }

    readAlignments(chr, bpStart, bpEnd) {

        var self = this;

        return new Promise(function (fulfill, reject) {


            getChrIndex(self).then(function (chrToIndex) {

                var chrId = chrToIndex[chr],

                    alignmentContainer = new self.ac_class(chr, bpStart, bpEnd, self.samplingWindowSize, self.samplingDepth, self.pairsSupported,self.parent);
                    
                if (chrId === undefined) {
                    fulfill(alignmentContainer);
                } else {

                    getIndex(self).then(function (bamIndex) {

                        var chunks = bamIndex.blocksForRange(chrId, bpStart, bpEnd),
                            promises = [];


                        if (!chunks) {
                            fulfill(null);
                            reject("Error reading bam index");
                            return;
                        }
                        if (chunks.length === 0) {
                            fulfill(alignmentContainer);
                            return;
                        }

                        chunks.forEach(function (c) {

                            promises.push(new Promise(function (fulfill, reject) {

                                var fetchMin = c.minv.block,
                                    fetchMax = c.maxv.block + MAX_GZIP_BLOCK_SIZE,   // Make sure we get the whole block.
                                    range = {start: fetchMin, size: fetchMax - fetchMin + 1};

                                igvxhr.loadArrayBuffer(self.bamPath,
                                    {
                                        headers: self.config.headers,
                                        range: range,
                                        withCredentials: self.config.withCredentials
                                    }).then(function (compressed) {

                                    var ba = new Uint8Array(new unbgzf(compressed)); //new Uint8Array(igv.unbgzf(compressed)); //, c.maxv.block - c.minv.block + 1));
                                    decodeBamRecords(ba, c.minv.offset, alignmentContainer, bpStart, bpEnd, chrId, self.filter);

                                    fulfill(alignmentContainer);

                                }).catch(function (obj) {
                                    reject(obj);
                                });

                            }))
                        });


                        Promise.all(promises).then(function (ignored) {
                            alignmentContainer.finish();
                            fulfill(alignmentContainer);
                        }).catch(function (obj) {
                            reject(obj);
                        });
                    }).catch(reject);
                }
            }).catch(reject);
        });


        function decodeBamRecords(ba, offset, alignmentContainer, min, max, chrId, filter) {

            var blockSize,
                blockEnd,
                alignment,
                blocks,
                refID,
                pos,
                bmn,
                bin,
                mq,
                nl,
                flag_nc,
                flag,
                nc,
                lseq,
                mateRefID,
                matePos,
                readName,
                j,
                p,
                lengthOnRef,
                cigar,
                c,
                cigarArray,
                seq,
                seqBytes;
            //let reads =0;
            while (true) {
                //reads++;
                blockSize = readInt(ba, offset);
                blockEnd = offset + blockSize + 4;

                if (blockEnd > ba.length) {
                    return;
                }

                alignment = new BamAlignment();

                refID = readInt(ba, offset + 4);
                pos = readInt(ba, offset + 8);

                if(refID < 0) {
                    return;   // unmapped reads
                }
                else if (refID > chrId || pos > max) {
                    return;    // off right edge, we're done
                }
                else if (refID < chrId) {
                    continue;   // to left of start, not sure this is possible
                }
               
                bmn = readInt(ba, offset + 12);
                bin = (bmn & 0xffff0000) >> 16;
                mq = (bmn & 0xff00) >> 8;
                nl = bmn & 0xff;

                flag_nc = readInt(ba, offset + 16);
                flag = (flag_nc & 0xffff0000) >> 16;
                nc = flag_nc & 0xffff;

                alignment.flags = flag;
                alignment.strand = !(flag & READ_STRAND_FLAG);

                lseq = readInt(ba, offset + 20);

                mateRefID = readInt(ba, offset + 24);
                matePos = readInt(ba, offset + 28);
                alignment.fragmentLength = readInt(ba, offset + 32);

                readName = '';
                for (j = 0; j < nl - 1; ++j) {
                    readName += String.fromCharCode(ba[offset + 36 + j]);
                }
                //console.log(reads+readName);

                p = offset + 36 + nl;

                lengthOnRef = 0;
                cigar = '';


                cigarArray = [];
                for (c = 0; c < nc; ++c) {
                    var cigop = readInt(ba, p);
                    var opLen = (cigop >> 4);
                    var opLtr = CIGAR_DECODER[cigop & 0xf];
                    if (opLtr == 'M' || opLtr == 'EQ' || opLtr == 'X' || opLtr == 'D' || opLtr == 'N' || opLtr == '=')
                        lengthOnRef += opLen;
                    cigar = cigar + opLen + opLtr;
                    p += 4;

                    cigarArray.push({len: opLen, ltr: opLtr});
                }
                alignment.cigar = cigar;
                alignment.lengthOnRef = lengthOnRef;

                if (alignment.start + alignment.lengthOnRef < min) continue;  // Record out-of-range "to the left", skip to next one


                seq = '';
                seqBytes = (lseq + 1) >> 1;
                for (j = 0; j < seqBytes; ++j) {
                    var sb = ba[p + j];
                    seq += SECRET_DECODER[(sb & 0xf0) >> 4];
                    seq += SECRET_DECODER[(sb & 0x0f)];
                }
                seq = seq.substring(0, lseq);  // seq might have one extra character (if lseq is an odd number)

                p += seqBytes;
                alignment.seq = seq;


                if (lseq === 1 && String.fromCharCode(ba[p + j] + 33) === "*") {
                    // TODO == how to represent this?
                }
                else {
                    alignment.qual = [];
                    for (j = 0; j < lseq; ++j) {
                        alignment.qual.push(ba[p + j]);
                    }
                }
                p += lseq;


                alignment.start = pos;
                alignment.mq = mq;
                alignment.readName = readName;
                alignment.chr = self.indexToChr[refID];

                if (mateRefID >= 0) {
                    alignment.mate = {
                        chr: self.indexToChr[mateRefID],
                        position: matePos,
                        strand: !(flag & MATE_STRAND_FLAG)
                    };
                }


                alignment.tagBA = new Uint8Array(ba.buffer.slice(p, blockEnd));  // decode thiese on demand
                p += blockEnd;

                if (!min || alignment.start <= max &&
                    alignment.start + alignment.lengthOnRef >= min &&
                    filter.pass(alignment)) {
                    if (chrId === undefined || refID == chrId) {
                        blocks = makeBlocks(alignment, cigarArray);
                        alignment.blocks = blocks.blocks;
                        alignment.insertions = blocks.insertions;
                        alignmentContainer.push(alignment);
                    }
                }
                offset = blockEnd;
            }
            
            // Exits via top of loop.
        }

        /**
         * Split the alignment record into blocks as specified in the cigarArray.  Each aligned block contains
         * its portion of the read sequence and base quality strings.  A read sequence or base quality string
         * of "*" indicates the value is not recorded.  In all other cases the length of the block sequence (block.seq)
         * and quality string (block.qual) must == the block length.
         *
         * NOTE: Insertions are not yet treated // TODO
         *
         * @param record
         * @param cigarArray
         * @returns array of blocks
         */
        function makeBlocks(record, cigarArray) {

            var blocks = [],
                insertions,
                seqOffset = 0,
                pos = record.start,
                len = cigarArray.length,
                blockSeq,
                blockQuals,
                gapType,
                minQ = 5,  //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN)
                maxQ = 20; //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX)

            for (var i = 0; i < len; i++) {

                var c = cigarArray[i];

                switch (c.ltr) {
                    case 'H' :
                        break; // ignore hard clips
                    case 'P' :
                        break; // ignore pads
                    case 'S' :
                        seqOffset += c.len;
                        gapType = 'S';
                        break; // soft clip read bases
                    case 'N' :
                        pos += c.len;
                        gapType = 'N';
                        break;  // reference skip
                    case 'D' :
                        pos += c.len;
                        gapType = 'D';
                        break;
                    case 'I' :
                        blockSeq = record.seq === "*" ? "*" : record.seq.substr(seqOffset, c.len);
                        blockQuals = record.qual ? record.qual.slice(seqOffset, c.len) : undefined;
                        if (insertions === undefined) insertions = [];
                        insertions.push({start: pos, len: c.len, seq: blockSeq, qual: blockQuals});
                        seqOffset += c.len;
                        break;
                    case 'M' :
                    case 'EQ' :
                    case '=' :
                    case 'X' :

                        blockSeq = record.seq === "*" ? "*" : record.seq.substr(seqOffset, c.len);
                        blockQuals = record.qual ? record.qual.slice(seqOffset, c.len) : undefined;
                        blocks.push({start: pos, len: c.len, seq: blockSeq, qual: blockQuals, gapType: gapType});
                        seqOffset += c.len;
                        pos += c.len;

                        break;

                    default :
                        console.log("Error processing cigar element: " + c.len + c.ltr);
                }
            }

            return {blocks: blocks, insertions: insertions};

        }
    }

    readHeader() {

        var self = this;

        return new Promise(function (fulfill, reject) {

            getIndex(self).then(function (index) {

                var len = index.firstAlignmentBlock + MAX_GZIP_BLOCK_SIZE;   // Insure we get the complete compressed block containing the header

                igvxhr.loadArrayBuffer(self.bamPath,
                    {
                        headers: self.config.headers,

                        range: {start: 0, size: len},

                        withCredentials: self.config.withCredentials
                    }).then(function (compressedBuffer) {

                    var unc = new unbgzf(compressedBuffer, len),
                        uncba = new Uint8Array(unc),
                        magic = readInt(uncba, 0),
                        samHeaderLen = readInt(uncba, 4),
                        samHeader = '',
                        genome = null;

                    for (var i = 0; i < samHeaderLen; ++i) {
                        samHeader += String.fromCharCode(uncba[i + 8]);
                    }

                    var nRef = readInt(uncba, samHeaderLen + 8);
                    var p = samHeaderLen + 12;

                    self.chrToIndex = {};
                    self.indexToChr = [];
                    for (var i = 0; i < nRef; ++i) {
                        var lName = readInt(uncba, p);
                        var name = '';
                        for (var j = 0; j < lName - 1; ++j) {
                            name += String.fromCharCode(uncba[p + 4 + j]);
                        }
                        var lRef = readInt(uncba, p + lName + 4);
                        //dlog(name + ': ' + lRef);

                        if (genome && genome.getChromosomeName) {
                            name = genome.getChromosomeName(name);
                        }

                        self.chrToIndex[name] = i;
                        self.indexToChr.push(name);

                        p = p + 8 + lName;
                    }

                    fulfill();

                }).catch(reject);
            }).catch(reject);
        });
    }

}
    function getIndex(bam) {

        return new Promise(function (fulfill, reject) {

            if (bam.index) {
                fulfill(bam.index);
            }
            else {
                loadBamIndex(bam.baiPath, bam.config).then(function (index) {
                    bam.index = index;

                    fulfill(bam.index);
                }).catch(reject);
            }
        });
    }


    function getChrIndex(bam) {

        return new Promise(function (fulfill, reject) {

            if (bam.chrToIndex) {
                fulfill(bam.chrToIndex);
            }
            else {
                bam.readHeader().then(function () {
                    fulfill(bam.chrToIndex);
                }).catch(reject);
            }
        });
    }

    function readInt(ba, offset) {
        return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
    }

    function readShort(ba, offset) {
        return (ba[offset + 1] << 8) | (ba[offset]);
    }




class BamSource{
    constructor(config,parent) {

        this.config = config;
        this.alignmentContainer = undefined;
        this.maxRows = config.maxRows || 1000;
        if (config.seq_url){
            this.sequence_source=new FastaSequence(config.seq_url);
        }
        this.parent=parent;

        this.pack_alignments=true;
       

        if (config.sourceType === "ga4gh") {
            this.bamReader = new igv.Ga4ghAlignmentReader(config);
        }
        else {
            this.bamReader = new BamReader(config,parent);
        }

       this.viewAsPairs = true;
    };

    setViewAsPairs(bool) {
        var self = this;

        if (this.viewAsPairs !== bool) {
            this.viewAsPairs = bool;
            // TODO -- repair alignments
            if (this.alignmentContainer) {
                var alignmentContainer = this.alignmentContainer,
                    alignments;

                if (bool) {
                    alignments = pairAlignments(alignmentContainer.packedAlignmentRows);
                }
                else {
                    alignments = unpairAlignments(alignmentContainer.packedAlignmentRows);
                }
                alignmentContainer.packedAlignmentRows = packAlignmentRows(alignments, alignmentContainer.start, alignmentContainer.end, self.maxRows);

            }
        }
    }

    getAlignments(chr, bpStart, bpEnd) {

        var self = this;
        return new Promise(function (fulfill, reject) {

            if (self.alignmentContainer && self.alignmentContainer.contains(chr, bpStart, bpEnd)) {
                fulfill(self.alignmentContainer);
            } else {

                self.bamReader.readAlignments(chr, bpStart, bpEnd).then(function (alignmentContainer) {

                    var maxRows = self.config.maxRows || 500,
                        alignments = alignmentContainer.alignments;

                    if (!self.viewAsPairs) {
                        alignments = unpairAlignments([{alignments: alignments}]);
                    }
                    if (self.parent.display_alignments){
                         alignmentContainer.packedAlignmentRows = packAlignmentRows(alignments, alignmentContainer.start, alignmentContainer.end, maxRows);
                    }
                   


                    alignmentContainer.alignments = undefined;  // Don't need to hold onto these anymore
                    self.alignmentContainer = alignmentContainer;
                    if (self.sequence_source){
                         self.sequence_source.getSequence(alignmentContainer.chr, alignmentContainer.start, alignmentContainer.end).then(
                        function (sequence) {


                            if (sequence) {

                                alignmentContainer.coverageMap.refSeq = sequence;    // TODO -- fix this
                                alignmentContainer.sequence = sequence;           // TODO -- fix this


                                fulfill(alignmentContainer);
                            }
                        }).catch(reject);

                    }
                    else{
                        fulfill(alignmentContainer);
                    }
                  

                }).catch(reject);
            }
        });
    }
}

    function pairAlignments(rows) {

        var pairCache = {},
            result = [];

        rows.forEach(function (row) {

            row.alignments.forEach(function (alignment) {

                var pairedAlignment;

                if (canBePaired(alignment)) {

                    pairedAlignment = pairCache[alignment.readName];
                    if (pairedAlignment) {
                        pairedAlignment.setSecondAlignment(alignment);
                        pairCache[alignment.readName] = undefined;   // Don't need to track this anymore.
                    }
                    else {
                        pairedAlignment = new igv.PairedAlignment(alignment);
                        pairCache[alignment.readName] = pairedAlignment;
                        result.push(pairedAlignment);
                    }
                }

                else {
                    result.push(alignment);
                }
            });
        });
        return result;
    }

    function unpairAlignments(rows) {
        var result = [];
        rows.forEach(function (row) {
            row.alignments.forEach(function (alignment) {
                if (alignment instanceof PairedAlignment) {
                    if (alignment.firstAlignment) result.push(alignment.firstAlignment);  // shouldn't need the null test
                    if (alignment.secondAlignment) result.push(alignment.secondAlignment);

                }
                else {
                    result.push(alignment);
                }
            });
        });
        return result;
    }

    function canBePaired(alignment) {
        return alignment.isPaired() &&
            alignment.isMateMapped() &&
            alignment.chr === alignment.mate.chr &&
            (alignment.isFirstOfPair() || alignment.isSecondOfPair()) && !(alignment.isSecondary() || alignment.isSupplementary());
    }


    function packAlignmentRows(alignments, start, end, maxRows) {

        if (!alignments) return;

        alignments.sort(function (a, b) {
            return a.start - b.start;
        });

        if (alignments.length === 0) {

            return [];

        } else {

            var bucketList = [],
                allocatedCount = 0,
                lastAllocatedCount = 0,
                nextStart = start,
                alignmentRow,
                index,
                bucket,
                alignment,
                alignmentSpace = 4 * 2,
                packedAlignmentRows = [],
                bucketStart = Math.max(start, alignments[0].start);

            alignments.forEach(function (alignment) {

                var buckListIndex = Math.max(0, alignment.start - bucketStart);
                if (bucketList[buckListIndex] === undefined) {
                    bucketList[buckListIndex] = [];
                }
                bucketList[buckListIndex].push(alignment);
            });


            while (allocatedCount < alignments.length && packedAlignmentRows.length < maxRows) {

                alignmentRow = new BamAlignmentRow();

                while (nextStart <= end) {

                    bucket = undefined;

                    while (!bucket && nextStart <= end) {

                        index = nextStart - bucketStart;
                        if (bucketList[index] === undefined) {
                            ++nextStart;                     // No alignments at this index
                        } else {
                            bucket = bucketList[index];
                        }

                    } // while (bucket)

                    if (!bucket) {
                        break;
                    }
                    alignment = bucket.pop();
                    if (0 === bucket.length) {
                        bucketList[index] = undefined;
                    }

                    alignmentRow.alignments.push(alignment);
                    nextStart = alignment.start + alignment.lengthOnRef + alignmentSpace;
                    ++allocatedCount;

                } // while (nextStart)

                if (alignmentRow.alignments.length > 0) {
                    packedAlignmentRows.push(alignmentRow);
                }

                nextStart = bucketStart;

                if (allocatedCount === lastAllocatedCount) break;   // Protect from infinite loops

                lastAllocatedCount = allocatedCount;

            } // while (allocatedCount)

            return packedAlignmentRows;
        }
    }



class BamAlignment{
    constructor(){
        this.hidden = false;
    }

    isMapped() {
        return (this.flags & READ_UNMAPPED_FLAG) == 0;
    }

    isPaired () {
        return (this.flags & READ_PAIRED_FLAG) != 0;
    }

    isProperPair () {
        return (this.flags & PROPER_PAIR_FLAG) != 0;
    }

    isFirstOfPair() {
        return (this.flags & FIRST_OF_PAIR_FLAG) != 0;
    }

    isSecondOfPair() {
        return (this.flags & SECOND_OF_PAIR_FLAG) != 0;
    }

    isSecondary() {
        return (this.flags & SECONDARY_ALIGNMNET_FLAG) != 0;
    }

    isSupplementary() {
        return (this.flags & SUPPLEMENTARY_ALIGNMENT_FLAG) != 0;
    }

    isFailsVendorQualityCheck() {
        return (this.flags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }

    isDuplicate() {
        return (this.flags & DUPLICATE_READ_FLAG) != 0;
    }

    isMateMapped() {
        return (this.flags & MATE_UNMAPPED_FLAG) == 0;
    }

    isNegativeStrand() {
        return (this.flags & READ_STRAND_FLAG) != 0;
    }

    isMateNegativeStrand() {
        return (this.flags & MATE_STRAND_FLAG) != 0;
    }

    tags() {

        function decodeTags(ba) {

            var p = 0,
                len = ba.length,
                tags = {};

            while (p < len) {
                var tag = String.fromCharCode(ba[p]) + String.fromCharCode(ba[p + 1]);
                var type = String.fromCharCode(ba[p + 2]);
                var value;

                if (type == 'A') {
                    value = String.fromCharCode(ba[p + 3]);
                    p += 4;
                } else if (type === 'i' || type === 'I') {
                    value = readInt(ba, p + 3);
                    p += 7;
                } else if (type === 'c' || type === 'C') {
                    value = ba[p + 3];
                    p += 4;
                } else if (type === 's' || type === 'S') {
                    value = readShort(ba, p + 3);
                    p += 5;
                } else if (type === 'f') {
                    // TODO 'FIXME need floats';
                    value = readFloat(ba, p + 3);
                    p += 7;
                } else if (type === 'Z') {
                    p += 3;
                    value = '';
                    for (; ;) {
                        var cc = ba[p++];
                        if (cc === 0) {
                            break;
                        } else {
                            value += String.fromCharCode(cc);
                        }
                    }
                } else {
                    //'Unknown type ' + type;
                    value = 'Error unknown type: ' + type;
                    tags[tag] = value;
                    break;
                }
                tags[tag] = value;
            }
            return tags;
        }

        if (!this.tagDict) {
            if (this.tagBA) {
                this.tagDict = decodeTags(this.tagBA);
                this.tagBA = undefined;
            } else {
                this.tagDict = {};  // Mark so we don't try again.  The record has not tags
            }
        }
        return this.tagDict;

    }

    popupData(genomicLocation) {

        // if the user clicks on a base next to an insertion, show just the
        // inserted bases in a popup (like in desktop IGV).
        var nameValues = [], isFirst, tagDict;

        if(this.insertions) {
            for(var i = 0; i < this.insertions.length; i += 1) {
                var ins_start = this.insertions[i].start;
                if(genomicLocation == ins_start || genomicLocation == ins_start - 1) {
                    nameValues.push({name: 'Insertion', value: this.insertions[i].seq });
                    nameValues.push({name: 'Location', value: ins_start });
                    return nameValues;
                }
            }
        }

        nameValues.push({ name: 'Read Name', value: this.readName });

        // Sample
        // Read group
        nameValues.push("<hr>");

        // Add 1 to genomic location to map from 0-based computer units to user-based units
        nameValues.push({ name: 'Alignment Start', value: igv.numberFormatter(1 + this.start), borderTop: true });

        nameValues.push({ name: 'Read Strand', value: (true === this.strand ? '(+)' : '(-)'), borderTop: true });
        nameValues.push({ name: 'Cigar', value: this.cigar });
        nameValues.push({ name: 'Mapped', value: yesNo(this.isMapped()) });
        nameValues.push({ name: 'Mapping Quality', value: this.mq });
        nameValues.push({ name: 'Secondary', value: yesNo(this.isSecondary()) });
        nameValues.push({ name: 'Supplementary', value: yesNo(this.isSupplementary()) });
        nameValues.push({ name: 'Duplicate', value: yesNo(this.isDuplicate()) });
        nameValues.push({ name: 'Failed QC', value: yesNo(this.isFailsVendorQualityCheck()) });

        if (this.isPaired()) {
            nameValues.push("<hr>");
            nameValues.push({ name: 'First in Pair', value: !this.isSecondOfPair(), borderTop: true });
            nameValues.push({ name: 'Mate is Mapped', value: yesNo(this.isMateMapped()) });
            if (this.isMateMapped()) {
                nameValues.push({ name: 'Mate Chromosome', value: this.mate.chr });
                nameValues.push({ name: 'Mate Start', value: (this.mate.position + 1)});
                nameValues.push({ name: 'Mate Strand', value: (true === this.mate.strand ? '(+)' : '(-)')});
                nameValues.push({ name: 'Insert Size', value: this.fragmentLength });
                // Mate Start
                // Mate Strand
                // Insert Size
            }
            // First in Pair
            // Pair Orientation

        }

        nameValues.push("<hr>");
        tagDict = this.tags();
        isFirst = true;
        for (var key in tagDict) {

            if (tagDict.hasOwnProperty(key)) {

                if (isFirst) {
                    nameValues.push({ name: key, value: tagDict[key], borderTop: true });
                    isFirst = false;
                } else {
                    nameValues.push({ name: key, value: tagDict[key] });
                }

            }
        }

        return nameValues;


        function yesNo(bool) {
            return bool ? 'Yes' : 'No';
        }
    }
}


   

   

    function readFloat(ba, offset) {

        var dataView = new DataView(ba.buffer),
            littleEndian = true;

        return dataView.getFloat32(offset, littleEndian);

    }




class BamFilter{

  constructor (options) {
        if (!options) options = {};
        this.vendorFailed = options.vendorFailed === undefined ? true : options.vendorFailed;
        this.duplicates = options.duplicates === undefined ? true : options.duplicates;
        this.secondary = options.secondary || false;
        this.supplementary = options.supplementary || false;
        this.mqThreshold = options.mqThreshold === undefined ? 0 : options.mqThreshold;
    }

    pass(alignment) {

        if (this.vendorFailed && alignment.isFailsVendorQualityCheck()) return false;
        if (this.duplicates && alignment.isDuplicate()) return false;
        if (this.secondary && alignment.isSecondary()) return false;
        if (this.supplementary && alignment.isSupplementary()) return false;
        if (alignment.mq < this.mqThreshold) return false;

        return true;


    }

}



    var BLOCK_HEADER_LENGTH = 18;
    var BLOCK_LENGTH_OFFSET = 16;  // Location in the gzip block of the total block size (actually total block size - 1)
    var BLOCK_FOOTER_LENGTH = 8; // Number of bytes that follow the deflated data
    var MAX_COMPRESSED_BLOCK_SIZE = 64 * 1024; // We require that a compressed block (including header and footer, be <= this)
    var GZIP_OVERHEAD = BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH + 2; // Gzip overhead is the header, the footer, and the block size (encoded as a short).
    var GZIP_ID1 = 31;   // Magic number
    var GZIP_ID2 = 139;  // Magic number
    var GZIP_FLG = 4; // FEXTRA flag means there are optional fields


    // Uncompress data,  assumed to be series of bgzipped blocks
    // Code is based heavily on bam.js, part of the Dalliance Genome Explorer,  (c) Thomas Down 2006-2001.
class unbgzf{
    constructor(data, lim) {
        for (let a in data){
           console.log(a);
         }
        var oBlockList = [],
            ptr = [0],
            totalSize = 0;

        lim = lim || data.byteLength - 18;

        while (ptr[0] < lim) {

            var ba = new Uint8Array(data, ptr[0], 18);

            var xlen = (ba[11] << 8) | (ba[10]);
            var si1 = ba[12];
            var si2 = ba[13];
            var slen = (ba[15] << 8) | (ba[14]);
            var bsize = (ba[17] << 8) | (ba[16]) + 1;

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
            return oBlockList[0];
        } else {
            var out = new Uint8Array(totalSize);
            var cursor = 0;
            for (var i = 0; i < oBlockList.length; ++i) {
                var b = new Uint8Array(oBlockList[i]);
                arrayCopy(b, 0, out, cursor, b.length);
                cursor += b.length;
            }
            return out.buffer;
        }
    }
}



class AlignmentContainer{
    constructor(chr, start, end, samplingWindowSize, samplingDepth, pairsSupported,parent) {
        this.parent=parent;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.length = (end - start);
        if (!this.parent.cm_class){
            this.coverageMap= new CoverageMap(chr,start,end,parent)
        }
        else{
            this.coverageMap= new this.parent.cm_class(chr,start,end,parent);
        }

       
        this.alignments = [];
        this.raw_alignments=[];
        this.downsampledIntervals = [];

        this.samplingWindowSize = samplingWindowSize === undefined ? 100 : samplingWindowSize;
        this.samplingDepth = samplingDepth === undefined ? 50 : samplingDepth;

        this.pairsSupported = pairsSupported;
        this.paired = false;  // false until proven otherwise
        this.pairsCache = {};  // working cache of paired alignments by read name

        this.downsampledReads = new Set();

        this.currentBucket = new DownsampleBucket(this.start, this.start + this.samplingWindowSize, this);

        this.filter = function filter(alignment) {         // TODO -- pass this in
            return alignment.isMapped() && !alignment.isFailsVendorQualityCheck();
        }

    }

    push(alignment) {

        if (this.filter(alignment) === false) return;
        if (alignment.tagBA){
            alignment.tagBA=decodeTags(alignment.tagBA);
        }
        if (this.parent.keep_raw_alignments){
            this.raw_alignments.push(alignment);
        }

        this.coverageMap.incCounts(alignment);
        if (!this.parent.display_alignments){
            return;
        }

         // Count coverage before any downsampling

        if (this.pairsSupported && this.downsampledReads.has(alignment.readName)) {
            return;   // Mate already downsampled -- pairs are treated as a single alignment for downsampling
        }

        if (alignment.start >= this.currentBucket.end) {
            finishBucket.call(this);
            this.currentBucket = new DownsampleBucket(alignment.start, alignment.start + this.samplingWindowSize, this);
        }

        this.currentBucket.addAlignment(alignment);

    }

    forEach(callback) {
        this.alignments.forEach(callback);
    }

    finish() {

        if (this.currentBucket !== undefined) {
            finishBucket.call(this);
        }

        // Need to remove partial pairs whose mate was downsampled
        if(this.pairsSupported) {
            var tmp = [], ds = this.downsampledReads;

            this.alignments.forEach(function (a) {
                if (!ds.has(a.readName)) {
                    tmp.push(a);
                }
            })
            this.alignments = tmp;
        }

        this.alignments.sort(function (a, b) {
            return a.start - b.start
        });

        this.pairsCache = undefined;
        this.downsampledReads = undefined;
    }

    contains(chr, start, end) {
        return this.chr == chr &&
            this.start <= start &&
            this.end >= end;
    }

    hasDownsampledIntervals() {
        return this.downsampledIntervals && this.downsampledIntervals.length > 0;
    }
}

    function finishBucket() {
        this.alignments = this.alignments.concat(this.currentBucket.alignments);
        if (this.currentBucket.downsampledCount > 0) {
            this.downsampledIntervals.push(new DownsampledInterval(
                this.currentBucket.start,
                this.currentBucket.end,
                this.currentBucket.downsampledCount));
        }
        this.paired = this.paired || this.currentBucket.paired;
    }


   
    



class DownsampleBucket{
    constructor(start, end, alignmentContainer) {

        this.start = start;
        this.end = end;
        this.alignments = [];
        this.downsampledCount = 0;
        this.samplingDepth = alignmentContainer.samplingDepth;
        this.pairsSupported = alignmentContainer.pairsSupported;
        this.downsampledReads = alignmentContainer.downsampledReads;
        this.pairsCache = alignmentContainer.pairsCache;
    }

    addAlignment(alignment) {

        var samplingProb, idx, replacedAlignment, pairedAlignment;

        if (this.alignments.length < this.samplingDepth) {

            if (this.pairsSupported && canBePaired(alignment)) {
                pairedAlignment = this.pairsCache[alignment.readName];
                if (pairedAlignment) {
                    //Not subject to downsampling, just update the existing alignment
                    pairedAlignment.setSecondAlignment(alignment);
                    this.pairsCache[alignment.readName] = undefined;   // Don't need to track this anymore. NOTE: Don't "delete", causes runtime performance issues
                }
                else {
                    // First alignment in a pair
                    pairedAlignment = new PairedAlignment(alignment);
                    this.paired = true;
                    this.pairsCache[alignment.readName] = pairedAlignment;
                    this.alignments.push(pairedAlignment);
                }
            }
            else {
                this.alignments.push(alignment);
            }

        } else {

            samplingProb = this.samplingDepth / (this.samplingDepth + this.downsampledCount + 1);

            if (Math.random() < samplingProb) {

                idx = Math.floor(Math.random() * (this.alignments.length - 1));
                replacedAlignment = this.alignments[idx];   // To be replaced

                if (this.pairsSupported && canBePaired(alignment)) {

                    if(this.pairsCache[replacedAlignment.readName] !== undefined) {
                        this.pairsCache[replacedAlignment.readName] = undefined;
                    }

                    pairedAlignment = new PairedAlignment(alignment);
                    this.paired = true;
                    this.pairsCache[alignment.readName] = pairedAlignment;
                    this.alignments[idx] = pairedAlignment;

                }
                else {
                    this.alignments[idx] = alignment;
                }
                this.downsampledReads.add(replacedAlignment.readName);

            }
            else {
                this.downsampledReads.add(alignment.readName);
            }

            this.downsampledCount++;
        }

    }
}


    // TODO -- refactor this to use an object, rather than an array,  if end-start is > some threshold
class CoverageMap{
    constructor(chr, start, end,parent) {

        this.chr = chr;
        this.bpStart = start;
        this.length = (end - start);
        this.parent=parent;
        this.coverage = new Array(this.length);

        this.maximum = 0;

       
    }

    incCounts(alignment) {

        var self = this;

        if (alignment.blocks === undefined) {

            incBlockCount(alignment);
        }
        else {
            alignment.blocks.forEach(function (block) {
                incBlockCount(block);
            });
        }

        function incBlockCount(block) {

            var key,
                base,
                i,
                j,
                q;

            for (i = block.start - self.bpStart, j = 0; j < block.len; i++, j++) {

                if (!self.coverage[i]) {
                    self.coverage[i] = new Coverage();
                }

                base = block.seq.charAt(j);
                key = (alignment.strand) ? "pos" + base : "neg" + base;
                q = block.qual[j];

                self.coverage[i][key] += 1;
                self.coverage[i]["qual" + base] += q;

                self.coverage[i].total += 1;
                self.coverage[i].qual += q;

                self.maximum = Math.max(self.coverage[i].total, self.maximum);

            }
        }
    }
}

CoverageMap.threshold = 0.2;
CoverageMap.qualityWeight = true;

class Coverage{
    constructor() {
        this.posA = 0;
        this.negA = 0;

        this.posT = 0;
        this.negT = 0;

        this.posC = 0;
        this.negC = 0;
        this.posG = 0;

        this.negG = 0;

        this.posN = 0;
        this.negN = 0;

        this.pos = 0;
        this.neg = 0;

        this.qualA = 0;
        this.qualT = 0;
        this.qualC = 0;
        this.qualG = 0;
        this.qualN = 0;

        this.qual = 0;

        this.total = 0;
    }

    isMismatch(refBase) {

        var myself = this,
            mismatchQualitySum,
            threshold = CoverageMap.threshold * ((CoverageMap.qualityWeight && this.qual) ? this.qual : this.total);

        mismatchQualitySum = 0;
        ["A", "T", "C", "G"].forEach(function (base) {

            if (base !== refBase) {
                mismatchQualitySum += ((CoverageMap.qualityWeight && myself.qual) ? myself["qual" + base] : (myself["pos" + base] + myself["neg" + base]));
            }
        });

        return mismatchQualitySum >= threshold;

    }
}

class DownsampledInterval{
    constructor(start, end, counts) {
        this.start = start;
        this.end = end;
        this.counts = counts;
    }

    popupData(genomicLocation) {
        return [
            {name: "start", value: this.start + 1},
            {name: "end", value: this.end},
            {name: "# downsampled:", value: this.counts}]
    }
}


class PairedAlignment{
    constructor(firstAlignment) {

        this.firstAlignment = firstAlignment;
        this.chr = firstAlignment.chr;
        this.readName = firstAlignment.readName;

        if (firstAlignment.start < firstAlignment.mate.position) {
            this.start = firstAlignment.start;
            this.end = Math.max(firstAlignment.mate.position, firstAlignment.start + firstAlignment.lengthOnRef);  // Approximate
            this.connectingStart = firstAlignment.start + firstAlignment.lengthOnRef;
            this.connectingEnd = firstAlignment.mate.position;
        }
        else {
            this.start = firstAlignment.mate.position;
            this.end = firstAlignment.start + firstAlignment.lengthOnRef;
            this.connectingStart = firstAlignment.mate.position;
            this.connectingEnd = firstAlignment.start;
        }
        this.lengthOnRef = this.end - this.start;

    }

    setSecondAlignment(alignment) {

        // TODO -- check the chrs are equal,  error otherwise
        this.secondAlignment = alignment;

        if (alignment.start > this.firstAlignment.start) {
            this.end = alignment.start + alignment.lengthOnRef;
            this.connectingEnd = alignment.start;
        }
        else {
            this.start = alignment.start;
            this.connectingStart = alignment.start + alignment.lengthOnRef;
        }
        this.lengthOnRef = this.end - this.start;


    }

    popupData(genomicLocation) {

        var nameValues = [];

        nameValues = nameValues.concat(this.firstAlignment.popupData(genomicLocation));

        if (this.secondAlignment) {
            nameValues.push("-------------------------------");
            nameValues = nameValues.concat(this.secondAlignment.popupData(genomicLocation));
        }
        return nameValues;
    }

    isPaired () {
        return true; // By definition
    }

    firstOfPairStrand () {
        if (this.firstAlignment.isFirstOfPair()) {
            return this.firstAlignment.strand;
        }
        else if (this.secondAlignment) {
            return this.secondAlignment.strand;
        }
        else {
            return this.firstAlignment.strand;          // This assumes inward pointing pairs
        }
    }

}

class BamAlignmentRow {
    constructor(){
        this.alignments = [];
        this.score = undefined;
    }

    findCenterAlignment(bpStart, bpEnd) {

        var centerAlignment = undefined;

        // find single alignment that overlaps sort location
        this.alignments.forEach(function(a){

            if (undefined === centerAlignment) {

                if ((a.start + a.lengthOnRef) < bpStart || a.start > bpEnd) {
                    // do nothing
                } else {
                    centerAlignment = a;
                }

            }

        });

        return centerAlignment;
    }

    updateScore(genomicLocation, genomicInterval, sortOption) {

        this.score = this.caculateScore(genomicLocation, (1 + genomicLocation), genomicInterval, sortOption);

    };

    caculateScore(bpStart, bpEnd, genomicInterval, sortOption) {

        var baseScore,
            alignment;

        alignment = this.findCenterAlignment(bpStart, bpEnd);
        if (undefined === alignment) {
            return Number.MAX_VALUE;
        }

        if ("NUCLEOTIDE" === sortOption.sort) {

            baseScore = undefined;

            alignment.blocks.forEach(function (block) {

                var sequence = genomicInterval.sequence,
                    coverageMap = genomicInterval.coverageMap,
                    reference,
                    base,
                    coverage,
                    count,
                    phred;

                if ("*" !== block.seq) {

                    for (var i = 0, indexReferenceSequence = block.start - genomicInterval.start, bpBlockSequence = block.start, lengthBlockSequence = block.seq.length;
                         i < lengthBlockSequence;
                         i++, indexReferenceSequence++, bpBlockSequence++) {

                        if (bpStart === bpBlockSequence) {

                            reference = sequence.charAt(indexReferenceSequence);
                            base = block.seq.charAt(i);

                            if (base === "=") {
                                base = reference;
                            }

                            if (base === 'N') {
                                baseScore = 2;
                            }
                            else if (base === reference) {
                                baseScore = 3;
                            }
                            else if (base === "X" || base !== reference){

                                coverage = coverageMap.coverage[ (bpBlockSequence - coverageMap.bpStart) ];
                                count = coverage[ "pos" + base ] + coverage[ "neg" + base ];
                                phred = (coverage.qual) ? coverage.qual : 0;
                                baseScore = -(count + (phred / 1000.0));
                            } else {
                                console.log("BamAlignmentRow.caculateScore - huh?");
                            }

                        } // bpStart === bpBlockSequence

                    } // block.seq.length

                }
                else {
                    baseScore = 3;
                }

            });

            return (undefined === baseScore) ? Number.MAX_VALUE : baseScore;
        }
        else if ("STRAND" === sortOption.sort) {

            return alignment.strand ? 1 : -1;
        }
        else if ("START" === sortOption.sort) {

            return alignment.start;
        }

        return Number.MAX_VALUE;

    }
}

let bgzBlockSize =function(data) {
    const ba = new Uint8Array(data);
    const bsize = (ba[17] << 8) | (ba[16]) + 1;
    return bsize;
}


function decodeTags(ba) {

            var p = 0,
                len = ba.length,
                tags = {};

            while (p < len) {
                var tag = String.fromCharCode(ba[p]) + String.fromCharCode(ba[p + 1]);
                var type = String.fromCharCode(ba[p + 2]);
                var value;

                if (type == 'A') {
                    value = String.fromCharCode(ba[p + 3]);
                    p += 4;
                } else if (type === 'i' || type === 'I') {
                    value = readInt(ba, p + 3);
                    p += 7;
                } else if (type === 'c' || type === 'C') {
                    value = ba[p + 3];
                    p += 4;
                } else if (type === 's' || type === 'S') {
                    value = readShort(ba, p + 3);
                    p += 5;
                } else if (type === 'f') {
                    // TODO 'FIXME need floats';
                    value = readFloat(ba, p + 3);
                    p += 7;
                } else if (type === 'Z') {
                    p += 3;
                    value = '';
                    for (; ;) {
                        var cc = ba[p++];
                        if (cc === 0) {
                            break;
                        } else {
                            value += String.fromCharCode(cc);
                        }
                    }
                } else {
                    //'Unknown type ' + type;
                    value = 'Error unknown type: ' + type;
                    tags[tag] = value;
                    break;
                }
                tags[tag] = value;
            }
            return tags;
        }

export {loadBamIndex,BamReader,BamSource,BamFilter,BamAlignment,AlignmentContainer,PairedAlignment,bgzBlockSize,CoverageMap,Coverage,BamFeatureReader};