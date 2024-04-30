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




import {BWSource} from "./bigwig.js";
import {FeatureSource,FastaSequence,BigBedFeatureSource,TabixBedFeatureSource} from "./feature.js";


class MLVTrack{
	constructor(config){
		this.config=config;
	}

	_setFeatureSource(){
		//overriden by tracks with feature sources
	}

	drawScale(ctx){
		//overidden in tracks with scale
	}

	getConfig(){
		return JSON.parse(JSON.stringify(this.config));
	}

	setConfigAttribute(attribute,value){
        if (value == null){
            delete this.config[attribute];
            return;
        }
		this.config[attribute]=value;
    	if (attribute==="url"){
    		this._setFeatureSource();
    	}
    }

    setColorFunction(func){
		this.color_function=func;
	}
    getFeatureAt(genomicLocation, chr,yOffset, bpPerPixel){
    	return null;
    }

    setConfig(config){
    	if (this.config.url !== config.url){
    		this.config=config;
    		this._setFeatureSource();
    	}
    	else{
    		this._setFeatureSource();
    	}
    }

     /**
	* Retrieves the features requested, the default is just to get the features
	* from the feature source 
	* @param {string} chr - The chromosome 
	* @param {integer} start - The start of the range from which to obtain features
	* @param {integer} end - The end of the range from which to obtain features 
	* @param {boolean} force - If true then cached features should not be used
	* but fresh features retrieved
	* @param {Object} data - contains bp  ber pixel and width of the canvas 
	*/
	getFeatures (chr, bpStart, bpEnd,force,data) {
		return this.feature_source.getFeatures(chr,bpStart,bpEnd,force,data);
	}


    getSettings(panel){
        return [{
            type:"slider",
            min:10,
            max:1000,
            label:"Height",
            current_value:this.config.height,
            func:x=>{
                this.setConfigAttribute("height",x);
                panel.update()

            }
        }]
    }


	static calculateLabel(url){
		if (typeof url !== "string"){
			url = url[0];
		}
		let arr =url.split("/");
		let label= arr[arr.length-1];
		arr= label.split(".");
		label = arr[0];
		return label;
	}

	addExtraControls(div,panel){
	}

	static getTypeFromURL(url){
		let config={}
		if (typeof url !== "string"){
			return config;
		}
			if (url.endsWith("bw")){
				config.type="bigwig";
				config.format="wig";
				
			}
			else if (url.endsWith(".bed.gz")){
				config.type="bed";
				config.format="feature"

			}
			else if (url.endsWith(".bb") && !(config.type)){
				config.type="bigbed";
				config.format="feature"

			}
			else if (url.endsWith(".fasta")){
				config.type="fasta";
				config.format="sequence";
			}
			else if (url.endsWith(".bam")){
				config.type="bam";
				config.format="alignment";
			}
		return config;

	}

	static parseConfig(con){
		let config = JSON.parse(JSON.stringify(con));
		
		if (!(config.type) && config.url){
			let info = MLVTrack.getTypeFromURL(config.url);
			if (info.type){
				config.type=info.type;
				config.format=info.format;
			}
		}
		if (config.type==="bed" || config.type==="bigbed"){
			config.format="feature";
		}
		if (!config.short_label && config.url){
			config.short_label=MLVTrack.calculateLabel(config.url);
		}
		
		if (!config.track_id){
			config.track_id=config.url;
		}

		if (config.format==="feature"){
			config.displayMode = config.displayMode || "COLLAPSED";    // COLLAPSED | EXPANDED | SQUISHED
        	config.labelDisplayMode = "SLANT";
        	config.squishedCallHeight = config.squishedCallHeight || 30;
        	config.expandedCallHeight = config.expandedCallHeight || 15;
        	config.featureHeight=config.featureHeight || 12;
		
		}

		if (config.format==="wig" || config.type==="bigwig"){
			if (!config.scale){
				config.scale="dynamic";
			}
			if (!config.min_y){
				config.min_y=0;
			}
			if (!config.max_y){
				config.max_y=100;
			}
			if (!config.height){
				config.height=100;
			}
		}
		if (config.type==="bam"){
			if (!config.featureHeight){
				config.featureHeight=12;
			}
		}
		if (!config.height){
			config.height=100;
		}

		if (!config.color){
			if (config.type==="bam"){
				config.color="#D3D3D3";
			}
			else{
				//config.color="black";
			}
		}
		if (!config.opacity){
			config.opacity=1.0;
		}
           
		return config;
	}

    drawMessage(options){
        const ctx = options.context;
        ctx.save();
        ctx.font="15px Verdana";
        const dim = ctx.measureText(options.features);
        const y= options.top+(options.height/2)+20;
        ctx.fillStyle = options.textColor;
        for (let x=0;x<options.pixelWidth;x+=dim.width+10){
            ctx.fillText(options.features, x, y);
        }
        ctx.restore();
    }




	
	static getTrack(config){
		config=MLVTrack.parseConfig(config);
		let cl= MLVTrack.track_types[config.type];
		if (!cl){
			throw new Error("Track type "+ config.type+" not recognised")
		}
        // pjt: nb, simple bed FeatureSource is not used anywhere and we have no way of loading simple bed files
        // config.type = "bed" is used for tabix indexed bed files
		return new cl["class"](config);

	}
		
}

MLVTrack.custom_tracks={};

MLVTrack.track_types={}





//*******************js/rulerTrack.js**********************


class RulerTrack extends MLVTrack{
	constructor(config){
		if (!config){
			config={"track_id":"ruler"+RulerTrack.count,format:"ruler",short_label:"Ruler",type:"ruler"};
		}
		super(config);
        this.height = 30;
        this.config.height=30;
        this.name = "";
       
        this.disableButtons = true;
        this.ignoreTrackMenu = true;
        this.order = -Number.MAX_VALUE;
      
        RulerTrack.count++;
    }

    getFeatures(chr, bpStart, bpEnd) {

        return new Promise(function (fulfill, reject) {
            fulfill([]);
        });
    }

    drawFeatures(options) {

        const ctx = options.context;
        const fontStyle = { textAlign: 'center', font: '10px PT Sans', fillStyle: options.textColor, strokeStyle: options.textColor };

        const range = Math.floor(1100 * options.bpPerPixel);
        const ts = RulerTrack.findSpacing(range);
        const spacing = ts.majorTick;

        // Find starting point closest to the current origin
        let nTick = Math.floor(options.bpStart / spacing) - 1;
        let x = 0;
		let y_pos=options.top+this.height;
        //canvas.setProperties({textAlign: 'center'});
        Graphics.setProperties(ctx, fontStyle );
        const shim = 2;
        while (x < options.pixelWidth) {

            const l = Math.floor(nTick * spacing);

            x = Math.round(((l - 1) - options.bpStart + 0.5) / options.bpPerPixel);
            const chrPosition = formatNumber(l / ts.unitMultiplier, 0) + " " + ts.majorUnit;

            if (nTick % 1 == 0) {
                Graphics.fillText(ctx, chrPosition, x, y_pos - 15);
            }

            Graphics.strokeLine(ctx, x, y_pos - 10, x, y_pos - shim);

            nTick++;
        }
        Graphics.strokeLine(ctx, 0, y_pos - shim, options.pixelWidth, y_pos - shim);


        function formatNumber(anynum, decimal) {
            //decimal  - the number of decimals after the digit from 0 to 3
            //-- Returns the passed number as a string in the xxx,xxx.xx format.
            //anynum = eval(obj.value);
            let divider = 10;
            switch (decimal) {
                case 0:
                    divider = 1;
                    break;
                case 1:
                    divider = 10;
                    break;
                case 2:
                    divider = 100;
                    break;
                default:       //for 3 decimal places
                    divider = 1000;
            }

            const workNum = Math.abs((Math.round(anynum * divider) / divider));

            let workStr = "" + workNum

            if (workStr.indexOf(".") == -1) {
                workStr += "."
            }

            let dStr = workStr.substring(0, workStr.indexOf("."));
            const dNum = dStr - 0
            let pStr = workStr.substring(workStr.indexOf("."))

            while (pStr.length - 1 < decimal) {
                pStr += "0"
            }

            if (pStr == '.') pStr = '';

            //--- Adds a comma in the thousands place.
            if (dNum >= 1000) {
                var dLen = dStr.length
                dStr = parseInt("" + (dNum / 1000)) + "," + dStr.substring(dLen - 3, dLen)
            }

            //-- Adds a comma in the millions place.
            if (dNum >= 1000000) {
                dLen = dStr.length
                dStr = parseInt("" + (dNum / 1000000)) + "," + dStr.substring(dLen - 7, dLen)
            }
            let retval = dStr + pStr
            //-- Put numbers in parentheses if negative.
            if (anynum < 0) {
                retval = "(" + retval + ")";
            }

            //You could include a dollar sign in the return value.
            //retval =  "$"+retval
            return retval;
        }
        return y_pos;


    }
    static findSpacing(maxValue) {

        if (maxValue < 10) {
            return new TickSpacing(1, "", 1);
        }


        // Now man zeroes?
        var nZeroes = Math.floor(log10(maxValue));
        var majorUnit = "";
        var unitMultiplier = 1;
        if (nZeroes > 9) {
            majorUnit = "gb";
            unitMultiplier = 1000000000;
        }
        if (nZeroes > 6) {
            majorUnit = "mb";
            unitMultiplier = 1000000;
        } else if (nZeroes > 3) {
            majorUnit = "kb";
            unitMultiplier = 1000;
        }

        var nMajorTicks = maxValue / Math.pow(10, nZeroes - 1);
        if (nMajorTicks < 25) {
            return new TickSpacing(Math.pow(10, nZeroes - 1), majorUnit, unitMultiplier);
        } else {
            return new TickSpacing(Math.pow(10, nZeroes) / 2, majorUnit, unitMultiplier);
        }

        function log10(x) {
            var dn = Math.log(10);
            return Math.log(x) / dn;
        }
    }

}

RulerTrack.count=0;

MLVTrack.track_types["ruler"]={
	"class":RulerTrack,
	name:"Ruler"

}

class TickSpacing{
	constructor(majorTick, majorUnit, unitMultiplier) {
        this.majorTick = majorTick;
        this.majorUnit = majorUnit;
        this.unitMultiplier = unitMultiplier;
    }

}

 
class MLVBedTrack extends MLVTrack{
	constructor(config){
		super(config);
		this._setFeatureSource();
		this.filter_function=null;
		this.color_function=null;
		
	}

	_setFeatureSource(){
		this.feature_source= new TabixBedFeatureSource(this.config);
	}

	setFilterFunction(func){
		this.filter_function=func;
	}
	

	getCurrentFeatures(chr,start,end){
		return this.feature_source.featureCache.queryFeatures(chr,start,end);
	}
	
	drawFeatures(options) {
		let max_y_val=0;
        var track = this,
        	py,
            featureList = options.features,
            ctx = options.context,
            bpPerPixel = options.bpPerPixel,
            bpStart = options.bpStart,
            pixelWidth = options.pixelWidth,
            pixelHeight = options.pixelHeight,
            offset=0,
            bpEnd = bpStart + pixelWidth * bpPerPixel + 1;
	 	let top=0;
       	if(options.top){
           top=options.top;
       	}
       	this.top=top;
        let conf = this.config;
        let windowX = 0;
        let windowX1 = windowX + pixelWidth;

        this.config.squishedCallHeight = this.config.featureHeight+10;
        this.config.expandedCallHeight = (this.config.featureHeight/2)+2;
        let ki=null;
        if (featureList.length>50000){
        	ki=Math.round(featureList.length/50000)+1;
        }



        if (featureList) {
        	let opacity = this.config.opacity?this.config.opacity:1;
        	ctx.globalAlpha=opacity;

            for (var gene, i = 0, len = featureList.length; i < len; i++) {
            	if (ki && i%ki!==0){
            		continue;
            	}
                gene = featureList[i];
                if (this.filter_function && !this.filter_function(gene)){
                	gene.display=false;
                	continue;
                }
                gene.display=true;
                if (gene.end < bpStart) continue;
                if (gene.start > bpEnd) break;
                let coord = this.calculateFeatureCoordinates(gene, bpStart,bpPerPixel);
                let h = conf.featureHeight
                if (conf.displayMode === "SQUISHED" && gene.row != undefined) {
            		h = conf.featureHeight / 2;
            		py = conf.expandedCallHeight * gene.row + 2;
        		} else if (conf.displayMode === "EXPANDED" && gene.row != undefined) {
            		py = conf.squishedCallHeight * gene.row + 5;
        		} else {
             // collapsed
           			 py = 5;
        		}
        		py+=top;
        		if (py+h>max_y_val){
        			max_y_val=py+h;
        		}
        		coord.py=py;
        		coord.h=h;
        		let info={bpPerPixel:bpPerPixel,bpStart:bpStart,pixelWidth:pixelWidth}
        		if (this.color_function){
        			gene.color=this.color_function(gene);
        		}
        		
                this.renderFeature(gene,coord,ctx,info,options.textColor);
                this.renderFeatureLabel(ctx, gene, coord.px, coord.px1, py, windowX, windowX1,options.textColor);
            }
            ctx.globalAlpha=1;
        }
        else {
            console.log("No feature list");
        }
        this.bottom=max_y_val;
        if (this.config.displayMode==="COLLAPSED"){
        	max_y_val+=25;
        }	
    	return max_y_val;
    }
        
           

    
  
	 /**
     * @param ctx       the canvas 2d context
     * @param feature
     * @param featureX  feature start x-coordinate
     * @param featureX1 feature end x-coordinate
     * @param featureY  feature y-coordinate
     * @param windowX   visible window start x-coordinate
     * @param windowX1  visible window end x-coordinate
     */
  



	renderFeatureLabel(ctx, feature, featureX, featureX1, featureY, windowX, windowX1,defaultColor) {
		let info = this.config;
        var geneFontStyle, transform,
            boxX, boxX1,    // label should be centered between these two x-coordinates
            labelX, labelY,
            textFitsInBox;

        // feature outside of viewable window
        if (featureX1 < windowX || featureX > windowX1) {
            boxX = featureX;
            boxX1 = featureX1;
        } else {
            // center label within visible portion of the feature
            boxX = Math.max(featureX, windowX);
            boxX1 = Math.min(featureX1, windowX1);
        }

        let text= feature.name;
        if (this.label_function){
        	text=this.label_function(feature);
        }
       

        //if (igv.browser.selection && "genes" === this.config.type && feature.name !== undefined) {
            // TODO -- for gtex, figure out a better way to do this
            //geneColor = igv.browser.selection.colorForGene(feature.name);
      //  }

        textFitsInBox = (boxX1 - boxX) > ctx.measureText(text).width;
        //geneColor="black";

        if ((textFitsInBox) && info.displayMode != "SQUISHED" && text) {
            geneFontStyle = {
                font: '10px PT Sans',
                textAlign: 'center',
                fillStyle: info.color || defaultColor,
                strokeStyle: info.color || defaultColor
            };

            if (info.displayMode === "COLLAPSED" && info.labelDisplayMode === "SLANT") {
                transform = {rotate: {angle: 45}};
                delete geneFontStyle.textAlign;
            }

            labelX = boxX + ((boxX1 - boxX) / 2);
            labelY = getFeatureLabelY(featureY, transform);

            Graphics.fillText(ctx, text, labelX, labelY, geneFontStyle, transform);
        }
        function getFeatureLabelY(featureY, transform) {
        	return transform ? featureY + info.featureHeight+14 : featureY + info.featureHeight+7;
    	}
    }



	calculateFeatureCoordinates(feature, bpStart, xScale) {
    	var px = Math.round((feature.start-1 - bpStart) / xScale),
        px1 = Math.round((feature.end - bpStart) / xScale),
        pw = px1 - px;

        if (pw < 3) {
        	pw = 3;
            px -= 1;
        }

        return {
        	px: px,
            px1: px1,
            pw: pw
        };
	}

	/**
     * Renders the feature to the canvas
     * @param feature - The feature itself
     * @param coord An object containing information on where to draw the feature
     * px1,px2 the left and right pixels - pw - the width
     * py the top, h - the height
     * @param ctx - The context to draw the feature
     * @param info - An object containing information about the genomic location
     * bpStart.bpPerPixel and pixelWidth
     * 
     */
	
	renderFeature(feature, coord,ctx,info,defaultColor){
		var e,x, cy, direction, exon, ePx, ePx1, ePxU, ePw, py2, h2, 
            step = 20,
            color = this.config.color || defaultColor;
        if (feature.color) {
            color = feature.color;
        }
  
        ctx.fillStyle = color;
        ctx.strokeStyle = color;


        cy = coord.py + coord.h / 2;
        h2 = coord.h / 2;
        py2 = cy - h2 / 2;

		let exonCount = feature.exons ? feature.exons.length : 0;
			if (exonCount == 0) {
            	// single-exon transcript
            	ctx.fillRect(coord.px, coord.py, coord.pw, coord.h);
        	}
        else{
            // multi-exon transcript
            coord.px=Math.max(coord.px,0);
            coord.px1=Math.min(coord.px1,info.pixelWidth);
            Graphics.strokeLine(ctx, coord.px + 1, cy, coord.px1 - 1, cy); // center line for introns
            direction = feature.strand == '+' ? 1 : -1;
            
           
          
            for ( x=coord.px + step / 2; x <  coord.px1; x += step) {

                // draw arrowheads along central line indicating transcribed orientation
                Graphics.strokeLine(ctx, x - direction * 2, cy - 2, x, cy);
                Graphics.strokeLine(ctx, x - direction * 2, cy + 2, x, cy);
            }
            for (e = 0; e < exonCount; e++) {
                // draw the exons
                exon = feature.exons[e];
                ePx = Math.round((exon.start - info.bpStart) / info.bpPerPixel);
                ePx1 = Math.round((exon.end - info.bpStart) / info.bpPerPixel);
                ePw = Math.max(1, ePx1 - ePx);

                if (exon.utr) {
                    ctx.fillRect(ePx, py2, ePw, h2); // Entire exon is UTR
                }
                else {
                    if (exon.cdStart) {
                        ePxU = Math.round((exon.cdStart - info.bpStart) / info.bpPerPixel);
                        ctx.fillRect(ePx, py2, ePxU - ePx, h2); // start is UTR
                        ePw -= (ePxU - ePx);
                        ePx = ePxU;

                    }
                    if (exon.cdEnd) {
                        ePxU = Math.round((exon.cdEnd - info.bpStart) / info.bpPerPixel);
                        ctx.fillRect(ePxU, py2, ePx1 - ePxU, h2); // start is UTR
                        ePw -= (ePx1 - ePxU);
                        ePx1 = ePxU;
                    }

                    ctx.fillRect(ePx, coord.py, ePw, coord.h);

                    // Arrows
                    if (ePw > step + 5) {
                        ctx.fillStyle = "white";
                        ctx.strokeStyle = "white";
                        for (x = ePx + step / 2; x < ePx1; x += step) {
                            // draw arrowheads along central line indicating transcribed orientation
                            Graphics.strokeLine(ctx, x - direction * 2, cy - 2, x, cy);
                            Graphics.strokeLine(ctx, x - direction * 2, cy + 2, x, cy);
                        }
                        ctx.fillStyle = color;
                        ctx.strokeStyle = color;

                    }
                }
            }
        }

	}

	getFeatureAt(genomicLocation, chr, coord, bpPerPixel) {
		let yOffset=coord.y-this.top;
        // We use the featureCache property rather than method to avoid async load.  If the
        // feature is not already loaded this won't work,  but the user wouldn't be mousing over it either.
        if (this.feature_source.featureCache) {

          
               var tolerance = 2 * bpPerPixel,  // We need some tolerance around genomicLocation, start with +/- 2 pixels
                featureList = this.feature_source.featureCache.queryFeatures(chr, genomicLocation - tolerance, genomicLocation + tolerance),
                row;

            if (this.config.displayMode != "COLLAPSED") {
                row = (Math.floor)(this.config.displayMode === "SQUISHED" ? yOffset / this.config.expandedCallHeight : yOffset / this.config.squishedCallHeight);
            }
            else{
                if (yOffset<5 || yOffset>5+this.config.featureHeight){
                    return;
                }
            }


            if (featureList && featureList.length > 0) {


                var popupData = [];
                for (let feature of featureList){
                    if (feature.end >= genomicLocation - tolerance &&
                        feature.start <= genomicLocation + tolerance) {

                        // If row number is specified use it
                        if ((row === undefined || feature.row === undefined || row === feature.row)&&  feature.display) {
                           return feature;

                        }
                    }
                }

               
            }

        }

        return null;
    };

}

MLVTrack.track_types["bed"]={
	"class":MLVBedTrack,
	name:"bed(tabix)"
}


class MLVBigBedTrack extends MLVBedTrack{
	constructor(config){
		super(config);
		
	}
	_setFeatureSource(){
		this.feature_source=new BigBedFeatureSource(this.config);
	}

}

MLVTrack.track_types["bigbed"]={
	"class":MLVBigBedTrack,
	"name":"BigBed"
	
}





class MLVWigTrack extends MLVTrack{
	constructor(config){
		config.format="wig";
		super(config);
		this._setFeatureSource();
	}

	_setFeatureSource(){
		this.feature_source=new BWSource(this.config);	
	}
    getSettings(panel) {
        return [
            ...super.getSettings(panel),
            {
                type: "check",
                label: "Show data range",
                value: this.config.showDataRange,
                func: v => {
                    this.setConfigAttribute("showDataRange", v);
                    panel.update();
                }
            }
        ]
    }

	drawScale(pixel_height,ctx,defaultColor){
		if (this.config.scale_link_to && this.config.group){
			return;
		}
	
		let	top=this.top;
		let	bot = this.bottom;
		
		ctx.fillStyle=defaultColor;
        ctx.strokeStyle=defaultColor;

		ctx.beginPath();
		ctx.moveTo(0,top);
		ctx.lineTo(0,bot);
		ctx.moveTo(0,top);
		ctx.lineTo(20,top);
		ctx.moveTo(0,bot);
		ctx.lineTo(20,bot);
		ctx.font="12px Arial";
		ctx.stroke();
		ctx.textBaseline="top";
		
		if (this.max_y != null){
		    ctx.fillText(this.max_y.toFixed(2),20,top);
		    ctx.textBaseline="alphabetic";
		    ctx.fillText(this.min_y,20,bot);
		}

	}



	
	
	
	
	drawFeatures(options) {
		let self = this,
	    features = options.features,
	    color=this.config.color,
	    ctx = options.context,
	    bpPerPixel = options.bpPerPixel,
	    bpStart = options.bpStart,
	    pixelWidth = options.pixelWidth,
	    pixelHeight =options.pixelHeight,
	    y_offset=options.top,
	    bpEnd = bpStart + pixelWidth * bpPerPixel + 1,
	    featureValueRange;
	    if (this.config.group){
            pixelHeight=options.height;
	    }
	    else {
	    	pixelHeight=this.config.height;	
	    }
	    if (!color){
	    	color="black";       
	    }
	    this.prev_coords={x:0,y:0};
		
        if (features) {
            if (this.scale_link_to) {
                let t = this.scale_link_to.config;
                if (t) {
                    this.config.scale = t.scale;
                    this.max_y = this.scale_link_to.max_y;
                    this.min_y = this.scale_link_to.min_y;
                }
            }
            else if (this.set_scale) {
                this.min_y = this.set_scale.min;
                this.max_y = this.set_scale.max;
            }
            else if ((this.max_y === undefined && this.config.scale === "automatic") || this.config.scale === "dynamic") {
                var s = autoscale(features);
                this.min_y = s.min;
                this.max_y = s.max;
            }
            else if (this.config.scale === "fixed") {
                this.min_y = this.config.min_y;
                this.max_y = this.config.max_y;
            }

            featureValueRange = this.max_y - this.min_y;

            //$dataRangeTrackLabel = $(this.trackView.trackDiv).find('.igv-data-range-track-label');
            //
            //min = (Math.floor(track.dataRange.min) === track.dataRange.min) ? track.dataRange.min : track.dataRange.min.toFixed(2);
            //max = (Math.floor(track.dataRange.max) === track.dataRange.max) ? track.dataRange.max : track.dataRange.max.toFixed(2);
            //str = '[' + min + ' - ' + max + ']';
            //
            //$dataRangeTrackLabel.text(str);
            let prev_x = 0;
            let prev_y = 0;
            ctx.globalAlpha = this.config.opacity ? this.config.opacity : 1;

            if (this.is_line) {
                let y = (1.0 - this.config.value / featureValueRange) * pixelHeight;
                Graphics.strokeLine(ctx, 0, y, pixelWidth, y, { "strokeStyle": this.config.color, "lineWidth": this.config.width ? this.config.width : 1 });
            }

            else {
                // pass on any error message from feature source...
                if (!Array.isArray(features)) {
                    throw features;
                }
                features.forEach(renderFeature);
                if (this.config.showDataRange) {
                    ctx.save();
                    features.forEach(f => {
                        f.valueBak = f.value;
                        f.value = f.minVal;
                    });
                    ctx.globalAlpha *= 0.4;
                    features.forEach(renderFeature);
                    features.forEach(f => f.value = f.maxVal);
                    features.forEach(renderFeature);
                    features.forEach(f => f.value = f.valueBak);
                    ctx.restore();
                }
            }

            ctx.globalAlpha = 1
            if (this.config.threshold) {
                let y = y_offset + (1.0 - this.config.threshold / featureValueRange) * pixelHeight;
                Graphics.strokeLine(ctx, 0, y, pixelWidth, y, { "strokeStyle": "black", "lineWidth": 1 });
            }
        }
        function renderFeature(feature) {

            let yUnitless,
                heightUnitLess,
                x,
                y,
                width,
                height,
                rectEnd;

            if (feature.end < bpStart) return;
            if (feature.start > bpEnd) return;
                if (feature.end===feature.start){
                feature.start-=1;
            }
            
            
            x = Math.floor((feature.start - bpStart) / bpPerPixel);
            
            rectEnd = Math.floor((feature.end - bpStart) / bpPerPixel);
            width = Math.max(0, rectEnd - x);

            //height = ((feature.value - featureValueMinimum) / featureValueRange) * pixelHeight;
            //rectBaseline = pixelHeight - height;
            //canvas.fillRect(rectOrigin, rectBaseline, rectWidth, rectHeight, {fillStyle: track.color});
            // feature.value = .... think about doing something to anti-alias.
            if (signsDiffer(self.min_y, self.max_y)) {

                if (feature.value < 0) {
                    yUnitless = self.max_y/ featureValueRange;
                    heightUnitLess = -feature.value / featureValueRange;
                } else {
                    yUnitless = ((self.max_y - feature.value) / featureValueRange);
                    heightUnitLess = feature.value / featureValueRange;
                }

            }
            else if (self.min_y < 0) {
                yUnitless = 0;
                heightUnitLess = -feature.value / featureValueRange;
            }
            else {
                yUnitless = 1.0 - ((feature.value-self.min_y) / featureValueRange);
                heightUnitLess = (feature.value+self.min_y) / featureValueRange;
            }

            y = (yUnitless*pixelHeight)+y_offset;
            y=y<y_offset?y_offset:y;
            height=heightUnitLess * pixelHeight;
            height=height>pixelHeight?pixelHeight:height

            //canvas.fillRect(x, yUnitless * pixelHeight, width, heightUnitLess * pixelHeight, { fillStyle: igv.randomRGB(64, 255) });
            if (self.config.display==='line'){
                if (self.prev_coords.x) {
                    Graphics.strokeLine(ctx,x,y,self.prev_coords.x,self.prev_coords.y,{"strokeStyle":color,"lineWidth":3});
                }
                self.prev_coords.x = x;
                self.prev_coords.y = y;
            }
            else{
                //phase data
                const ph = self.config.phase_data;
                if(ph){
                    const pos = (bpPerPixel*x)+bpStart;
                    if (pos >ph.st && pos <ph.en){
                        const pro = ph.ratio*height;
                        Graphics.fillRect(ctx, x, y, width, pro, {fillStyle: "blue"});
                        Graphics.fillRect(ctx, x, y+pro, width, height-pro, {fillStyle: "red"});
                    }
                    else{
                        Graphics.fillRect(ctx, x, y, width, height, {fillStyle: color});
                    }
                }
                else{
                    Graphics.fillRect(ctx, x, y, width, height, {fillStyle: color});
                }
            }           
        }
        function autoscale(features) {
            let min = Number.MAX_VALUE, max = -Number.MAX_VALUE;
            if (!Array.isArray(features)) {
                console.warn('features is not an array');
                return {min: 0, max: 100};
            }
            features.forEach(function (f) {
                min = Math.min(min, f.value);
                max = Math.max(max, f.value);
            });
            return {min: min, max: max};
        }

        function signsDiffer(a, b) {
            return (a > 0 && b < 0 || a < 0 && b > 0);
        }
        this.top=y_offset;
        this.bottom=y_offset+pixelHeight;
        
        return this.bottom;
	}
}



MLVTrack.track_types["bigwig"]={
	"class":MLVWigTrack,
	name:"BigWig"
}


 class LineTrack extends MLVWigTrack{
	 constructor(config){
		 super(config);
		 this.is_line=true;
	 }
	 
	 getFeatures(chr, bpStart, bpEnd) {
		 return new Promise(function (fulfill, reject) {
			 fulfill([]);
	     });
	 }
	 
	 
 }

MLVTrack.track_types["line"]={
	"class":LineTrack,
	name:"Line Track"
}




class SequenceTrack extends MLVTrack{
    constructor(config){
		super(config) 
        this._setFeatureSource(config);
        this.sequenceType = config.sequenceType || "dna";  
        this.height = 15;
    }

    _setFeatureSource(config){
    	this.feature_source = new FastaSequence(config.url);
    }



    getFeatures(chr, bpStart, bpEnd,force,data) {
		let self = this;
        return new Promise(function (fulfill, reject) {
            if (data.bpPerPixel > 1/*igv.browser.trackViewportWidthBP() > 30000*/) {
                fulfill(null);
            }
            else {
                self.feature_source.getSequence(chr, bpStart, bpEnd).then(fulfill).catch(reject);
            }
        });
    }


    drawFeatures(options) {
        var sequence = options.features,
            ctx = options.context,
            bpPerPixel = options.bpPerPixel,
            bpStart = options.bpStart,
            pixelWidth = options.pixelWidth,
            bpEnd = bpStart + pixelWidth * bpPerPixel + 1,
            len, w, y, pos, offset, b, p0, p1, pc, c;

        let y_pos1=options.top;
        let y_pos2=y_pos1+5;

        if (sequence) {
            len = sequence.length;
            w = 1 / bpPerPixel;

            y = y_pos1+this.height / 2;
            for (pos = bpStart; pos <= bpEnd; pos++) {
                offset = pos - bpStart;
                if (offset < len) {
                    // var b = sequence.charAt(offset);
                    b = sequence[offset];
                    p0 = Math.floor(offset * w);
                    p1 = Math.floor((offset + 1) * w);
                    pc = Math.round((p0 + p1) / 2);

                    if (this.color) {
                        c = this.color;
                    }
                    else if ("dna" === this.sequenceType) {
                        c = Graphics.nucleotideColors[b];
                    }
                    else {
                        c = "rgb(0, 0, 150)";
                    }

                    if (!c) c = "gray";

                    if (bpPerPixel >0.18) {
                        Graphics.fillRect(ctx, p0, y_pos1, p1 - p0, 10, {fillStyle: c});
                    }
                    else {
                        Graphics.strokeText(ctx, b, pc, y, {
                            strokeStyle: c,
                            font: 'normal 10px Arial',
                            textAlign: 'center'
                        });
                    }
                }
            }
        }
        return y_pos1+10;
    }

}


MLVTrack.track_types["fasta"]={
	"class":SequenceTrack,
	name:"Fasta"
}



//*******js/ifv-canvas.js***********************


class Graphics{


       static setProperties(ctx, properties) {

            for (var key in properties) {
                if (properties.hasOwnProperty(key)) {
                    var value = properties[key];
                    ctx[key] = value;
                }
            }
        }

        static strokeLine (ctx, x1, y1, x2, y2, properties) {

            x1 = Math.floor(x1) + 0.5;
            y1 = Math.floor(y1) + 0.5;
            x2 = Math.floor(x2) + 0.5;
            y2 = Math.floor(y2) + 0.5;

            //log("stroke line, prop: " + properties);

            ctx.save();
            if (properties) Graphics.setProperties(ctx, properties);

            ctx.beginPath();
            ctx.moveTo(x1, y1);
            ctx.lineTo(x2, y2);
            ctx.stroke();
            ctx.restore();
        }

        static fillRect(ctx, x, y, w, h, properties) {

            var c;
            x = Math.round(x);
            y = Math.round(y);

            if (properties) {
                ctx.save();
                Graphics.setProperties(ctx, properties);
            }
            
            ctx.fillRect(x, y, w, h);

            if (properties) ctx.restore();
        }

        static fillPolygon(ctx, x, y, properties) {
            ctx.save();
            if (properties)   Graphics.setProperties(ctx, properties);
            Graphics.doPath(ctx, x, y);
            ctx.fill();
            ctx.restore();
        }

        static strokePolygon(ctx, x, y, properties) {
            ctx.save();
            if (properties)   Graphics.setProperties(ctx, properties);
            Graphics.doPath(ctx, x, y);
            ctx.stroke();
            ctx.restore();
        }

        static fillText(ctx, text, x, y, properties, transforms) {

            if (properties) {
                ctx.save();
                Graphics.setProperties(ctx, properties);
            }


            ctx.save();

            ctx.translate(x, y);
            if (transforms) {

                for (var transform in transforms) {
                    var value = transforms[transform];

                    // TODO: Add error checking for robustness
                    if (transform == 'translate') {
                        ctx.translate(value['x'], value['y']);
                    }
                    if (transform == 'rotate') {
                        ctx.rotate(value['angle'] * Math.PI / 180);
                    }
                }

            }

            ctx.fillText(text, 0, 0);
            ctx.restore();

            if (properties) ctx.restore();

        }

        static strokeText(ctx, text, x, y, properties, transforms) {


            ctx.save();
            if (properties) {
                Graphics.setProperties(ctx, properties);
            }


            ctx.translate(x, y);
            if (transforms) {

                for (var transform in transforms) {
                    var value = transforms[transform];

                    // TODO: Add error checking for robustness
                    if (transform == 'translate') {
                        ctx.translate(value['x'], value['y']);
                    }
                    if (transform == 'rotate') {
                        ctx.rotate(value['angle'] * Math.PI / 180);
                    }
                }
            }


            ctx.strokeText(text, 0, 0);
            ctx.restore();

        }

        static strokeCircle(ctx, x, y, radius) {

            ctx.beginPath();
            ctx.arc(x, y, radius, 0, 2 * Math.PI);
            ctx.stroke();
        }

        static fillCircle (ctx, x, y, radius) {

            ctx.beginPath();
            ctx.arc(x, y, radius, 0, 2 * Math.PI);
            ctx.fill();
        }

        static drawArrowhead(ctx, x, y, size, lineWidth) {

            ctx.save();
            if (!size) {
                size = 5;
            }
            if (lineWidth) {
                ctx.lineWidth = lineWidth;
            }
            ctx.beginPath();
            ctx.moveTo(x, y - size / 2);
            ctx.lineTo(x, y + size / 2);
            ctx.lineTo(x + size, y);
            ctx.lineTo(x, y - size / 2);
            ctx.closePath();
            ctx.fill();
            ctx.restore();
        }

        static dashedLine(ctx, x1, y1, x2, y2, dashLen, properties) {
            ctx.save();
            x1 = Math.round(x1);
            y1 = Math.round(y1);
            x2 = Math.round(x2);
            y2 = Math.round(y2);
            dashLen = Math.round(dashLen);
            log("dashedLine");
            if (properties) Graphics.setProperties(ctx, properties);

            if (dashLen == undefined) dashLen = 2;
            ctx.moveTo(x1, y1);

            var dX = x2 - x1;
            var dY = y2 - y1;
            var dashes = Math.floor(Math.sqrt(dX * dX + dY * dY) / dashLen);
            var dashX = dX / dashes;
            var dashY = dY / dashes;

            var q = 0;
            while (q++ < dashes) {
                x1 += dashX;
                y1 += dashY;
                ctx[q % 2 == 0 ? 'moveTo' : 'lineTo'](x1, y1);
            }
            ctx[q % 2 == 0 ? 'moveTo' : 'lineTo'](x2, y2);

            ctx.restore();
        }


    

        static doPath(ctx, x, y) {


        	var i, len = x.length;
        	for (i = 0; i < len; i++) {
        		x[i] = Math.round(x[i]);
        		y[i] = Math.round(y[i]);
        	}

        	ctx.beginPath();
        	ctx.moveTo(x[0], y[0]);
        	for (i = 1; i < len; i++) {
        		ctx.lineTo(x[i], y[i]);
        	}
        	ctx.closePath();
        }

}

Graphics.nucleotideColors={
	"A":"green",
	"T":"red",
	"G":"gray",
	"C":"blue",
	"a":"green",
	"t":"red",
	"c":"blue",
	"g":"gray"

}



export {MLVTrack,MLVWigTrack,MLVBedTrack,RulerTrack,MLVBigBedTrack,Graphics}