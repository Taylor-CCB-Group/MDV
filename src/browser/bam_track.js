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

import {MLVTrack,Graphics,MLVWigTrack} from "./tracks.js";
import {BamSource,PairedAlignment} from "./bam.js";
import {Utils} from "./utils.js";

var alignmentRowYInset = 0;
var alignmentStartGap = 5;
var downsampleRowHeight = 5;
const DEFAULT_COVERAGE_TRACK_HEIGHT = 50;







class ATACBamTrack extends MLVTrack{
    constructor(config) {
        if (!config.color){
            config.color="#D3D3D3";
        }
        super(config);
        this._setFeatureSource();
        if (this.config.bigwig){
            this.bigwig_track = new MLVWigTrack({
                url:this.config.bigwig,
                type:"bigwig",
                discrete:true,
                color:"#D3D3D3",
                scale:"dynamic",
                height:100
            })
        }
    }


    _setFeatureSource(){
        this.feature_source=new BamSource(this.config,this);
        //this.feature_source.setViewAsPairs(true);
	    //this.coverageTrack = new CoverageTrack(this.config, this);
        //this.alignmentTrack = new AlignmentTrack(this.config, this);
    }

    getFeatures(chr, bpStart, bpEnd,force,data) {
        /*if (bpEnd- bpStart>500000 && this.bigwig_track){
            this.draw_bigwig=true;
            return this.bigwig_track.getFeatures(chr, bpStart, bpEnd,force,data);
        }
        this.draw_bigwig=false;*/
        return this.feature_source.getAlignments(chr, bpStart, bpEnd);
    }


}

class BAMTrack extends MLVTrack{

    constructor(config) {
        if (!config.color){
            config.color="#D3D3D3";
        }
        super(config);
        this._setFeatureSource();
        if (this.config.bigwig){
            this.bigwig_track = new MLVWigTrack({
                url:this.config.bigwig,
                type:"bigwig",
                discrete:true,
                color:"#D3D3D3",
                scale:"dynamic",
                height:100
            })
        }

        

        this.config=config;

        if(config.coverageTrackHeight === undefined) {
            config.coverageTrackHeight = DEFAULT_COVERAGE_TRACK_HEIGHT;
        }

        

        this.visibilityWindow = config.visibilityWindow || 30000;     // 30kb default

        this.viewAsPairs = true;//config.viewAsPairs;

        this.pairsSupported = config.pairsSupported === undefined ? true : false;

        this.color = config.color || "rgb(185, 185, 185)";

        // sort alignment rows
        this.sortOption = config.sortOption || {sort: "NUCLEOTIDE"};
        this.sortDirection = true;

        // filter alignments
        this.filterOption = config.filterOption || {name: "mappingQuality", params: [30, undefined]};
        this.display_alignments=true;

    }

    _setFeatureSource(){
        this.feature_source=new BamSource(this.config,this);
        this.feature_source.setViewAsPairs(true);

	   this.coverageTrack = new CoverageTrack(this.config, this);
        this.alignmentTrack = new AlignmentTrack(this.config, this);
    }

    getFeatures(chr, bpStart, bpEnd,force,data) {
        if (bpEnd- bpStart>500000 && this.bigwig_track){
            this.draw_bigwig=true;
            return this.bigwig_track.getFeatures(chr, bpStart, bpEnd,force,data);
        }
        this.draw_bigwig=false;
        return this.feature_source.getAlignments(chr, bpStart, bpEnd);
    }

    // Alt - Click to Sort alignment rows
    altClick(genomicLocation, event) {

        this.alignmentTrack.sortAlignmentRows(genomicLocation, this.sortOption);

        this.trackView.redrawTile(this.featureSource.alignmentContainer);
        $(this.trackView.viewportDiv).scrollTop(0);

        this.sortDirection = !this.sortDirection;
    };


    /**
     * Optional method to compute pixel height to accomodate the list of features.  The implementation below
     * has side effects (modifiying the samples hash).  This is unfortunate, but harmless.
     *
     * @param alignmentContainer
     * @returns {number}
     */
    computePixelHeight(alignmentContainer) {

        return this.coverageTrack.computePixelHeight(alignmentContainer) +
            this.alignmentTrack.computePixelHeight(alignmentContainer);

    };

    drawFeatures(options) {
        this.top=options.top;
        if (this.draw_zoom_in){
            this.bottom=this.top+this.config.height;
            return this.bottom;
        }
        if (this.draw_bigwig){
            this.bottom=this.bigwig_track.drawFeatures(options);
            return;
        }
        if(this.coverageTrack.height > 0) {
            this.coverageTrack.draw(options);

        }

         this.bottom = this.alignmentTrack.draw(options);
         return this.bottom;
    };

    paintAxis(ctx, pixelWidth, pixelHeight) {

        this.coverageTrack.paintAxis(ctx, pixelWidth, this.coverageTrackHeight);

    }
    getFeatureAt(genomicLocation,chr,offset,bpPerPixel){
        if (!this.feature_source.alignmentContainer){
            return null;
        }
    	var packedAlignmentRows = this.feature_source.alignmentContainer.packedAlignmentRows,
            downsampledIntervals = this.feature_source.alignmentContainer.downsampledIntervals,
            packedAlignmentsIndex,
            alignmentRow,
            clickedObject,
            i, len, tmp;
        if (!packedAlignmentRows){
            return null;
        }
        packedAlignmentsIndex = Math.floor(((offset.y -this.top-this.coverageTrack.height) - (alignmentRowYInset)) / this.config.featureHeight);

        if (packedAlignmentsIndex < 0) {

            for (i = 0, len = downsampledIntervals.length; i < len; i++) {


                if (downsampledIntervals[i].start <= genomicLocation && (downsampledIntervals[i].end >= genomicLocation)) {
                    clickedObject = downsampledIntervals[i];
                    break;
                }

            }
        }
        else if (packedAlignmentsIndex < packedAlignmentRows.length) {

            alignmentRow = packedAlignmentRows[packedAlignmentsIndex];

            clickedObject = undefined;

            for (i = 0, len = alignmentRow.alignments.length, tmp; i < len; i++) {

                tmp = alignmentRow.alignments[i];

                if (tmp.start <= genomicLocation && (tmp.start + tmp.lengthOnRef >= genomicLocation)) {
                    clickedObject = tmp;
                    break;
                }

            }
        }

        return clickedObject;

    }


    popupData(genomicLocation, xOffset, yOffset) {

        if (yOffset >= this.coverageTrack.top && yOffset < this.coverageTrack.height) {
            return this.coverageTrack.popupData(genomicLocation, xOffset, this.coverageTrack.top);
        }
        else {
            return this.alignmentTrack.popupData(genomicLocation, xOffset, yOffset - this.alignmentTrack.top);
        }

    }

    popupMenuItems(popover) {

        var self = this,
            html,
            menuItems = [],
            colorByMenuItems = [],
            tagLabel = 'tag' + (self.alignmentTrack.colorByTag ? ' (' + self.alignmentTrack.colorByTag + ')' : ''),
            selected;


        colorByMenuItems.push({key: 'none', label: 'track color'});

        if(!self.viewAsPairs) {
            colorByMenuItems.push({key: 'strand', label: 'read strand'});
        }
        if (self.pairsSupported && self.alignmentTrack.hasPairs) {
            colorByMenuItems.push({key: 'firstOfPairStrand', label: 'first-of-pair strand'});
        }
        colorByMenuItems.push({key: 'tag', label: tagLabel});

        menuItems.push(igv.colorPickerMenuItem(popover, this.trackView));

        menuItems.push('<div class="igv-track-menu-category igv-track-menu-border-top">Color by</div>');

        colorByMenuItems.forEach(function (item) {
            selected = self.alignmentTrack.colorBy === item.key;
            menuItems.push(colorByMarkup(item, selected));
        });


        html = [];
        if (self.pairsSupported && self.alignmentTrack.hasPairs) {
            html.push('<div class="igv-track-menu-item igv-track-menu-border-top">');
            html.push(true === self.viewAsPairs ? '<i class="fa fa-check fa-check-shim">' : '<i class="fa fa-check fa-check-shim fa-check-hidden">');
            html.push('</i>');
            html.push('View as pairs');
            html.push('</div>');
            menuItems.push({
                object: $(html.join('')),
                click: function () {
                    var $fa = $(this).find('i');

                    popover.hide();

                    self.viewAsPairs = !self.viewAsPairs;

                    if (true === self.viewAsPairs) {
                        $fa.removeClass('fa-check-hidden');
                    } else {
                        $fa.addClass('fa-check-hidden');
                    }

                    self.featureSource.setViewAsPairs(self.viewAsPairs);
                    self.trackView.update();
                }
            });
        }

        return menuItems;

        function colorByMarkup(menuItem, showCheck, index) {

            var parts = [],
                item = {};


            //parts.push((0 === index) ? '<div class=\"igv-track-menu-item igv-track-menu-border-top\">' : '<div class="igv-track-menu-item">');
            parts.push('<div class="igv-track-menu-item">');

            parts.push(showCheck ? '<i class="fa fa-check fa-check-shim"></i>' : '<i class="fa fa-check fa-check-shim fa-check-hidden"></i>');

            //parts.push('<span>');
            //parts.push('Color by: ');
            //parts.push('</span>');

            if (menuItem.key === 'tag') {
                parts.push('<span id="color-by-tag">');
            } else {
                parts.push('<span>');
            }
            parts.push(menuItem.label);
            parts.push('</span>');

            parts.push('</div>');

            item.object = $(parts.join(''));

            item.click = function () {

                igv.popover.hide();

                if ('tag' === menuItem.key) {

                    igv.dialog.configure(function () {
                            return "Tag Name"
                        },

                        self.alignmentTrack.colorByTag ? self.alignmentTrack.colorByTag : '',

                        function () {
                            var tag = igv.dialog.$dialogInput.val().trim();
                            self.alignmentTrack.colorBy = 'tag';

                            if(tag !== self.alignmentTrack.colorByTag) {
                                self.alignmentTrack.colorByTag = igv.dialog.$dialogInput.val().trim();
                                self.alignmentTrack.tagColors = new igv.PaletteColorTable("Set1");
                                $('#color-by-tag').text(self.alignmentTrack.colorByTag);
                            }

                            self.trackView.update();
                        });

                    igv.dialog.show($(self.trackView.trackDiv));

                } else {
                    self.alignmentTrack.colorBy = menuItem.key;
                    self.trackView.update();
                }
            };

            return item;
        }
    }


   
}

 function shadedBaseColor(qual, nucleotide, genomicLocation) {

        var color,
            alpha,
            minQ = 5,   //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN),
            maxQ = 20,  //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX);
            foregroundColor =[0,255,255], //nucleotideColorComponents[nucleotide],
            backgroundColor = [255, 255, 255];   // White


        //if (171167156 === genomicLocation) {
        //    // NOTE: Add 1 when presenting genomic location
        //    console.log("shadedBaseColor - locus " + igv.numberFormatter(1 + genomicLocation) + " qual " + qual);
        //}

        if (!foregroundColor) return;

        if (qual < minQ) {
            alpha = 0.1;
        } else {
            alpha = Math.max(0.1, Math.min(1.0, 0.1 + 0.9 * (qual - minQ) / (maxQ - minQ)));
        }
        // Round alpha to nearest 0.1
        alpha = Math.round(alpha * 10) / 10.0;

        if (alpha >= 1) {
            color = Graphics.nucleotideColors[nucleotide];
        }
        else {
            color = "rgba(" + foregroundColor[0] + "," + foregroundColor[1] + "," + foregroundColor[2] + "," + alpha + ")";    //igv.getCompositeColor(backgroundColor, foregroundColor, alpha);
        }
        return color;
    }


class CoverageTrack{
	constructor(config, parent) {

        this.parent = parent;
        this.featureSource = parent.featureSource;
       


        this.height = 20;//config.coverageTrackHeight;
        this.dataRange = {min: 0};   // Leav max undefined
        
    }

    computePixelHeight(alignmentContainer) {
        return this.height;
    }

    draw(options) {

        var self = this,
            alignmentContainer = options.features,
            ctx = options.context,
            bpPerPixel = options.bpPerPixel,
            bpStart = options.bpStart,
            pixelWidth = options.pixelWidth,
            bpEnd = bpStart + pixelWidth * bpPerPixel + 1,
            coverageMap = alignmentContainer.coverageMap,
            bp,
            x,
            y,
            w,
            h,
            refBase,
            i,
            len,
            item,
            accumulatedHeight,
            sequence;

        

        if (coverageMap.refSeq) sequence = coverageMap.refSeq.toUpperCase();

        this.dataRange.max = coverageMap.maximum;
        var top =options.top;
        options.top+=this.height;
        // paint backdrop color for all coverage buckets
        w = Math.max(1, Math.ceil(1.0 / bpPerPixel));
        for (i = 0, len = coverageMap.coverage.length; i < len; i++) {

            bp = (coverageMap.bpStart + i);
            if (bp < bpStart) continue;
            if (bp > bpEnd) break;

            item = coverageMap.coverage[i];
            if (!item) continue;

            h = Math.round((item.total / this.dataRange.max) * this.height);
            y = this.height - h;
            x = Math.floor((bp - bpStart) / bpPerPixel);

            Graphics.setProperties(ctx, {fillStyle: this.parent.config.color, strokeStyle: this.color});
            Graphics.fillRect(ctx, x, y+top, w, h);
        }

        // coverage mismatch coloring -- don't try to do this in above loop, color bar will be overwritten when w<1
        if (sequence) {
            for (i = 0, len = coverageMap.coverage.length; i < len; i++) {

                bp = (coverageMap.bpStart + i);
                if (bp < bpStart) continue;
                if (bp > bpEnd) break;

                item = coverageMap.coverage[i];
                if (!item) continue;

                h = (item.total / this.dataRange.max) * this.height;
                y = this.height - h;
                x = Math.floor((bp - bpStart) / bpPerPixel);

                refBase = sequence[i];
                if (item.isMismatch(refBase)) {

                    Graphics.setProperties(ctx, {fillStyle: Graphics.nucleotideColors[refBase]});
                    Graphics.fillRect(ctx, x, y+top, w, h);

                    accumulatedHeight = 0.0;
                    ["A", "C", "T", "G"].forEach(function (nucleotide) {

                        var count,
                            hh;

                        count = item["pos" + nucleotide] + item["neg" + nucleotide];


                        // non-logoritmic
                        hh = (count / self.dataRange.max) * self.height;

                        y = (self.height - hh) - accumulatedHeight;
                        accumulatedHeight += hh;

                        Graphics.setProperties(ctx, {fillStyle: Graphics.nucleotideColors[nucleotide]});
                        Graphics.fillRect(ctx, x, y+top, w, hh);
                    });
                }
            }
        }

    }


    popupData(genomicLocation, xOffset, yOffset) {

        var coverageMap = this.featureSource.alignmentContainer.coverageMap,
            coverageMapIndex,
            coverage,
            nameValues = [];


        coverageMapIndex = genomicLocation - coverageMap.bpStart;
        coverage = coverageMap.coverage[coverageMapIndex];

        if (coverage) {


            nameValues.push(igv.browser.referenceFrame.chr + ":" + igv.numberFormatter(1 + genomicLocation));

            nameValues.push({name: 'Total Count', value: coverage.total});

            // A
            tmp = coverage.posA + coverage.negA;
            if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor(((coverage.posA + coverage.negA) / coverage.total) * 100.0) + "%)";
            nameValues.push({name: 'A', value: tmp});


            // C
            tmp = coverage.posC + coverage.negC;
            if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor((tmp / coverage.total) * 100.0) + "%)";
            nameValues.push({name: 'C', value: tmp});

            // G
            tmp = coverage.posG + coverage.negG;
            if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor((tmp / coverage.total) * 100.0) + "%)";
            nameValues.push({name: 'G', value: tmp});

            // T
            tmp = coverage.posT + coverage.negT;
            if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor((tmp / coverage.total) * 100.0) + "%)";
            nameValues.push({name: 'T', value: tmp});

            // N
            tmp = coverage.posN + coverage.negN;
            if (tmp > 0)  tmp = tmp.toString() + " (" + Math.floor((tmp / coverage.total) * 100.0) + "%)";
            nameValues.push({name: 'N', value: tmp});

        }


        return nameValues;

    }
}

class AlignmentTrack{
	constructor(config, parent) {
	    this.config=config;

        this.parent = parent;
        this.featureSource = parent.feature_source;
     
        this.alignmentRowHeight = config.alignmentRowHeight || 6;

        this.negStrandColor = config.negStrandColor || "rgba(150, 150, 230, 0.75)";
        this.posStrandColor = config.posStrandColor || "rgba(230, 150, 150, 0.75)";
        this.insertionColor = config.insertionColor || "rgb(138, 94, 161)";
        this.deletionColor = config.deletionColor || "black";
        this.skippedColor = config.skippedColor || "rgb(150, 170, 170)";

        this.colorBy = config.colorBy || "none";
        this.colorByTag = config.colorByTag;
        this.bamColorTag = config.bamColorTag === undefined ? "YC" : config.bamColorTag;

        // sort alignment rows
        this.sortOption = config.sortOption || {sort: "NUCLEOTIDE"};

        this.sortDirection = true;

        this.hasPairs = false;   // Until proven otherwise

    }

    computePixelHeight(alignmentContainer) {

        if (alignmentContainer.packedAlignmentRows) {
            var h = 0;
            if (alignmentContainer.hasDownsampledIntervals()) {
                h += downsampleRowHeight + alignmentStartGap;
            }
            return h + (this.config.featureHeight * alignmentContainer.packedAlignmentRows.length) + 5;
        }
        else {
            return this.height;
        }

    }

    draw(options) {

        var self = this,
            alignmentContainer = options.features,
            ctx = options.context,
            bpPerPixel = options.bpPerPixel,
            bpStart = options.bpStart,
            pixelWidth = options.pixelWidth,
            bpEnd = bpStart + pixelWidth * bpPerPixel + 1,
            packedAlignmentRows = alignmentContainer.packedAlignmentRows,
            sequence = alignmentContainer.sequence,
            base_text="normal "+(this.config.featureHeight-2)+"px Arial";
        if (!packedAlignmentRows){
            return;
        }
        let top = options.top;
        if (sequence) {
            sequence = sequence.toUpperCase();
        }

        if (alignmentContainer.hasDownsampledIntervals()) {
            alignmentRowYInset = downsampleRowHeight + alignmentStartGap;

            alignmentContainer.downsampledIntervals.forEach(function (interval) {
                var xBlockStart = (interval.start - bpStart) / bpPerPixel,
                    xBlockEnd = (interval.end - bpStart) / bpPerPixel;

                if (xBlockEnd - xBlockStart > 5) {
                    xBlockStart += 1;
                    xBlockEnd -= 1;
                }
                Graphics.fillRect(ctx, xBlockStart, 2, (xBlockEnd - xBlockStart), downsampleRowHeight - 2, {fillStyle: "black"});
            })

        }
        else {
            alignmentRowYInset = 3;
        }

        if (packedAlignmentRows) {

            packedAlignmentRows.forEach(function renderAlignmentRow(alignmentRow, i) {

                var yRect = alignmentRowYInset + (self.config.featureHeight * i)+top,
                    alignmentHeight = self.config.featureHeight - 2,
                    i,
                    b,
                    alignment;

                for (i = 0; i < alignmentRow.alignments.length; i++) {

                    alignment = alignmentRow.alignments[i];

                    self.hasPairs = self.hasPairs || alignment.isPaired();

                    if ((alignment.start + alignment.lengthOnRef) < bpStart) continue;
                    if (alignment.start > bpEnd) break;


                    if (true === alignment.hidden) {
                        continue;
                    }

                    if (alignment instanceof PairedAlignment) {

                        drawPairConnector(alignment, yRect, alignmentHeight);

                        drawSingleAlignment(alignment.firstAlignment, yRect, alignmentHeight);

                        if (alignment.secondAlignment) {
                            drawSingleAlignment(alignment.secondAlignment, yRect, alignmentHeight);
                        }

                    }
                    else {
                        drawSingleAlignment(alignment, yRect, alignmentHeight);
                    }

                }
            });
        }


        // alignment is a PairedAlignment
        function drawPairConnector(alignment, yRect, alignmentHeight) {


            var alignmentColor = self.getAlignmentColor(alignment.firstAlignment),
                outlineColor = alignmentColor,
                xBlockStart = (alignment.connectingStart - bpStart) / bpPerPixel,
                xBlockEnd = (alignment.connectingEnd - bpStart) / bpPerPixel,
                yStrokedLine = yRect + alignmentHeight / 2;
            
            if ((alignment.connectingEnd) < bpStart || alignment.connectingStart > bpEnd) return;

            if (alignment.mq <= 0) {
                alignmentColor = Utils.addAlphaToRGB(alignmentColor, "0.15");
            }

            Graphics.setProperties(ctx, {fillStyle: alignmentColor, strokeStyle: outlineColor});

            Graphics.strokeLine(ctx, xBlockStart, yStrokedLine, xBlockEnd, yStrokedLine);

        }


        function drawSingleAlignment(alignment, yRect, alignmentHeight) {
        

            var alignmentColor = self.getAlignmentColor(alignment),
                outlineColor = alignmentColor,
                lastBlockEnd,
                blocks = alignment.blocks,
                block,
                b;

                 
           if (!alignmentColor){
            alignmentColor="gray";
        }
      
            if ((alignment.start + alignment.lengthOnRef) < bpStart || alignment.start > bpEnd) return;

            if (alignment.mq <= 0) {
                alignmentColor = Utils.addAlphaToRGB(alignmentColor, "0.15");
            }

            Graphics.setProperties(ctx, {fillStyle: alignmentColor, strokeStyle: outlineColor});

            for (b = 0; b < blocks.length; b++) {   // Can't use forEach here -- we need ability to break

                block = blocks[b];

                if ((block.start + block.len) < bpStart) continue;

                drawBlock(block);

                if ((block.start + block.len) > bpEnd) break;  // Do this after drawBlock to insure gaps are drawn


                if (alignment.insertions) {
                    alignment.insertions.forEach(function (block) {
                        var refOffset = block.start - bpStart,
                            xBlockStart = refOffset / bpPerPixel - 1,
                            widthBlock = 3;
                        Graphics.fillRect(ctx, xBlockStart, yRect - 1, widthBlock, alignmentHeight + 2, {fillStyle: self.insertionColor});
                    });
                }

            }

            function drawBlock(block) {
                var seqOffset = block.start - alignmentContainer.start,
                    xBlockStart = (block.start - bpStart) / bpPerPixel,
                    xBlockEnd = ((block.start + block.len) - bpStart) / bpPerPixel,
                    widthBlock = Math.max(1, xBlockEnd - xBlockStart),
                    widthArrowHead = self.config.featureHeight / 2.0,
                    blockSeq = block.seq.toUpperCase(),
                    skippedColor = self.skippedColor,
                    deletionColor = self.deletionColor,
                    refChar,
                    readChar,
                    readQual,
                    xBase,
                    widthBase,
                    colorBase,
                    x,
                    y,
                    i,
                    yStrokedLine = yRect + alignmentHeight / 2;

                if (block.gapType !== undefined && xBlockEnd !== undefined && lastBlockEnd !== undefined) {
                    if ("D" === block.gapType) {
                        Graphics.strokeLine(ctx, lastBlockEnd, yStrokedLine, xBlockStart, yStrokedLine, {strokeStyle: deletionColor});
                    }
                    else {
                        Graphics.strokeLine(ctx, lastBlockEnd, yStrokedLine, xBlockStart, yStrokedLine, {strokeStyle: skippedColor});
                    }
                }
                lastBlockEnd = xBlockEnd;

                if (true === alignment.strand && b === blocks.length - 1) {
                    // Last block on + strand
                    x = [
                        xBlockStart,
                        xBlockEnd,
                        xBlockEnd + widthArrowHead,
                        xBlockEnd,
                        xBlockStart,
                        xBlockStart];
                    y = [
                        yRect,
                        yRect,
                        yRect + (alignmentHeight / 2.0),
                        yRect + alignmentHeight,
                        yRect + alignmentHeight,
                        yRect];
                    Graphics.fillPolygon(ctx, x, y, {fillStyle: alignmentColor});
                    if (alignment.mq <= 0) {
                        Graphics.strokePolygon(ctx, x, y, {strokeStyle: outlineColor});
                    }
                }
                else if (false === alignment.strand && b === 0) {
                    // First block on - strand
                    x = [
                        xBlockEnd,
                        xBlockStart,
                        xBlockStart - widthArrowHead,
                        xBlockStart,
                        xBlockEnd,
                        xBlockEnd];
                    y = [
                        yRect,
                        yRect,
                        yRect + (alignmentHeight / 2.0),
                        yRect + alignmentHeight,
                        yRect + alignmentHeight,
                        yRect];
                    Graphics.fillPolygon(ctx, x, y, {fillStyle: alignmentColor});
                    if (alignment.mq <= 0) {
                        Graphics.strokePolygon(ctx, x, y, {strokeStyle: outlineColor});
                    }
                }
                else {
                    //      igv.graphics.fillRect(ctx, xBlockStart, yRect, widthBlock, height, {fillStyle: "white"});
                    Graphics.fillRect(ctx, xBlockStart, yRect, widthBlock, alignmentHeight, {fillStyle: alignmentColor});
                    if (alignment.mq <= 0) {
                        ctx.save();
                        ctx.strokeStyle = outlineColor;
                        ctx.strokeRect(xBlockStart, yRect, widthBlock, alignmentHeight);
                        ctx.restore();
                    }
                }
                // Only do mismatch coloring if a refseq exists to do the comparison
                if (sequence && blockSeq !== "*") {
                    for (var i = 0, len = blockSeq.length; i < len; i++) {
                        readChar = blockSeq.charAt(i);
                        refChar = sequence.charAt(seqOffset + i);
                        if (readChar === "=") {
                            readChar = refChar;
                        }
                        if (readChar === "X" || refChar !== readChar) {
                            if (block.qual && block.qual.length > i) {
                                readQual = block.qual[i];
                                colorBase = shadedBaseColor(readQual, readChar, i + block.start);
                            }
                            else {
                                colorBase = Graphics.nucleotideColors[readChar];
                            }
                            if (colorBase) {
                                xBase = ((block.start + i) - bpStart) / bpPerPixel;
                                widthBase = Math.max(1, 1 / bpPerPixel);
                                if (bpPerPixel<0.15){
                                  Graphics.strokeText(ctx, readChar, xBase+widthBase/2, yRect+2, {
                                    strokeStyle: colorBase,
                                    font: base_text,
                                    textAlign: 'center',
                                    textBaseline:'hanging'
                                });
                                }
                                else{
                                Graphics.fillRect(ctx, xBase, yRect, widthBase, alignmentHeight, {fillStyle: colorBase});
                                }
                            }
                        }
                    }
                }
            }
        }
        return alignmentRowYInset + (self.config.featureHeight * packedAlignmentRows.length)+top;

    }

    sortAlignmentRows(genomicLocation, sortOption) {

        var self = this,
            alignmentContainer = this.featureSource.alignmentContainer,
            alignmentRows = alignmentContainer.packedAlignmentRows;

        alignmentRows.forEach(function (alignmentRow) {
            alignmentRow.updateScore(genomicLocation, alignmentContainer, sortOption);
        });

        alignmentRows.sort(function (a, b) {
            return self.sortDirection ? a.score - b.score : b.score - a.score;
        });

    }

    static doSortAlignmentRows(genomicLocation, genomicInterval, sortOption, sortDirection) {

        var alignmentRows = genomicInterval.packedAlignmentRows,
            sequence = genomicInterval.sequence;

        if (sequence) {
            sequence = sequence.toUpperCase();
        } else {
            console.log("No sequence, no traversal. No discussion!");
            return;
        }

        alignmentRows.forEach(function (alignmentRow) {
            alignmentRow.updateScore(genomicLocation, genomicInterval, sortOption);
        });

        alignmentRows.sort(function (a, b) {
            return sortDirection ? a.score - b.score : b.score - a.score;
        });

    }

    popupData(genomicLocation, xOffset, yOffset) {

        var packedAlignmentRows = this.featureSource.alignmentContainer.packedAlignmentRows,
            downsampledIntervals = this.featureSource.alignmentContainer.downsampledIntervals,
            packedAlignmentsIndex,
            alignmentRow,
            clickedObject,
            i, len, tmp;

        packedAlignmentsIndex = Math.floor((yOffset - (alignmentRowYInset)) / this.config.featureHeight);

        if (packedAlignmentsIndex < 0) {

            for (i = 0, len = downsampledIntervals.length; i < len; i++) {


                if (downsampledIntervals[i].start <= genomicLocation && (downsampledIntervals[i].end >= genomicLocation)) {
                    clickedObject = downsampledIntervals[i];
                    break;
                }

            }
        }
        else if (packedAlignmentsIndex < packedAlignmentRows.length) {

            alignmentRow = packedAlignmentRows[packedAlignmentsIndex];

            clickedObject = undefined;

            for (i = 0, len = alignmentRow.alignments.length, tmp; i < len; i++) {

                tmp = alignmentRow.alignments[i];

                if (tmp.start <= genomicLocation && (tmp.start + tmp.lengthOnRef >= genomicLocation)) {
                    clickedObject = tmp;
                    break;
                }

            }
        }

        if (clickedObject) {
            return clickedObject.popupData(genomicLocation);
        }
        else {
            return [];
        }

    };


    getAlignmentColor(alignment) {
       
        if (this.parent.color_function){
            return this.parent.color_function(alignment);
        }

        var alignmentTrack = this,
            option = alignmentTrack.colorBy,
            tagValue, color,
            strand;

        color = alignmentTrack.parent.config.color; // default

        switch (option) {

            case "strand":
                color = alignment.strand ? alignmentTrack.posStrandColor : alignmentTrack.negStrandColor;
                break;
            case "firstOfPairStrand":

                if(alignment instanceof PairedAlignment) {
                    color = alignment.firstOfPairStrand() ? alignmentTrack.posStrandColor : alignmentTrack.negStrandColor;
                }
                else if (alignment.isPaired()) {

                    if (alignment.isFirstOfPair()) {
                        color = alignment.strand ? alignmentTrack.posStrandColor : alignmentTrack.negStrandColor;
                    }
                    else if (alignment.isSecondOfPair()) {
                        color = alignment.strand ? alignmentTrack.negStrandColor : alignmentTrack.posStrandColor;
                    }
                    else {
                        console.log("ERROR. Paired alignments are either first or second.")
                    }
                }
                break;

            case "tag":

                tagValue = alignment.tags()[alignmentTrack.colorByTag];
                if (tagValue !== undefined) {

                    if (alignmentTrack.bamColorTag === alignmentTrack.colorByTag) {
                        // UCSC style color option
                        color = "rgb(" + tagValue + ")";
                    }
                    else {
                        color = alignmentTrack.tagColors.getColor(tagValue);
                    }
                }
                break;
            default:
                color = alignmentTrack.parent.config.color;
        }
        return color;

    }

}



BAMTrack.filters = {

        noop: function () {
            return function (alignment) {
                return false;
            };
        },

        strand: function (strand) {
            return function (alignment) {
                return alignment.strand === strand;
            };
        },

        mappingQuality: function (lower, upper) {
            return function (alignment) {

                if (lower && alignment.mq < lower) {
                    return true;
                }

                if (upper && alignment.mq > upper) {
                    return true;
                }

                return false;
            }
        }
    };

MLVTrack.custom_tracks["bam"]=BAMTrack;

MLVTrack.track_types["bam"]={
	"class":BAMTrack,
	name:"BAM"
}

export {BAMTrack,CoverageTrack,AlignmentTrack};
