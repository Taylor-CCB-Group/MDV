import {BAMTrack,CoverageTrack,AlignmentTrack} from "../bam_track.js";
import {BamSource,CoverageMap,Coverage} from "../bam.js";
import {MLVTrack} from '../tracks.js';

class SCACoverageMap{
            constructor(chr, start, end,parent) {
            
            this.chr = chr;
            this.bpStart = start;
            this.length = (end - start);
            this.parent =parent;
            this.reset();
          
         
            
       
        }
        reset(){
        	this.coverages= {"selected":null};
        
            for (let c in this.parent.categories){
            	   let a = new Array(this.length);
            	   a.maximum=0;
            	   a.count = this.parent.category_counts[c]	 
                   this.coverages[c]=a;
            }
            this.resetSelected();

        }

        resetSelected(){
        	  let a = new Array(this.length);
                a.maximum=0;
                a.count=this.parent.selected_count;
                this.coverages["selected"]=a;
        }

  

        incCounts(alignment,sel_only) {
        	let self = this;
            if (!sel_only){
				
				let cat = this.parent.category_function(alignment.tagBA.CB);
			

				if (alignment.blocks === undefined) {

					incBlockCount(alignment,cat);
				}
				else {
					alignment.blocks.forEach(function (block) {
						incBlockCount(block,cat);
					});
				}
            }

            let cat = "selected"

            if (this.parent.selected){
                if (!this.parent.selected[alignment.tagBA.CB]){
                	return
                }
            }
            

            if (alignment.blocks === undefined) {

                incBlockCount(alignment,cat);
            }
            else {
                alignment.blocks.forEach(function (block) {
                    incBlockCount(block,cat);
                });
            }

            


            function incBlockCount(block,cat) {
            	let coverage = self.coverages[cat];
            	if (!coverage){
            		return;
            	}

				var //key,
					//base,
					i,
					j,
					q;
			

				for (i = block.start - self.bpStart, j = 0; j < block.len; i++, j++) {

					if (!coverage[i]) {
						coverage[i] = {total:0,maximum:0};//new Coverage();
					}

					//base = block.seq.charAt(j);
					//key = (alignment.strand) ? "pos" + base : "neg" + base;
					//q = block.qual[j];

					//coverage[i][key] += 1;
					//coverage[i]["qual" + base] += q;

					coverage[i].total += 1;
					//coverage[i].qual += q;

					coverage.maximum = Math.max(coverage[i].total, coverage.maximum);

				}
            }
        }
    }

class SCACoverageTrack extends CoverageTrack{

     		constructor(coverage,height){
     			super(coverage,height);	
     		}

     		drawScale(ctx){
                let top=this.parent.top;
	            let cms = this.parent.feature_source.alignmentContainer.coverageMap;
	            let ind_track_height= this.parent.config.ind_track_height?this.parent.config.ind_track_height:50;
	            for (let name of this.parent.cat_order){
	            	let c= cms.coverages[name];
	            	let bot = top+ind_track_height;
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
		            ctx.fillStyle="black";
		            let num = c.maximum/c.count;
		            if (isNaN(num)){
		            	num=0;
		            }
		            ctx.fillText(num.toFixed(2),20,top);
		            ctx.font="14px Arial";
		            ctx.textBaseline="middle";
		            ctx.fillText(name,2,top+(ind_track_height/2));
		            ctx.font="12px Arial";
		            top+=ind_track_height;
	            }
     		}


     		draw(options){
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
                h,
                refBase,
                i,
                len,
                item,
                accumulatedHeight,
                sequence;

                let ind_track_height= this.parent.config.ind_track_height?this.parent.config.ind_track_height:50;

        

                if (coverageMap.refSeq) sequence = coverageMap.refSeq.toUpperCase();
				let top = options.top;
				let w = Math.max(1, Math.ceil(1.0 / bpPerPixel));
				for (let cb of this.parent.cat_order){

					let coverage= coverageMap.coverages[cb];
					console.log(cb+":"+coverage.count);
					let max = coverage.maximum/coverage.count
					ctx.fillStyle=this.parent.categories[cb];
					if (cb==="selected"){
						ctx.fillStyle="red";
					}
					for (i = 0, len = coverage.length; i < len; i++) {

						bp = (coverageMap.bpStart + i);
						if (bp < bpStart) continue;
						if (bp > bpEnd) break;

						item = coverage[i];
						if (!item) continue;

						h = Math.round(((item.total/coverage.count) / max) * ind_track_height);
						y = ind_track_height - h;
						x = Math.floor((bp - bpStart) / bpPerPixel);

						
						ctx.fillRect( x, y+top, w, h);
				    }
				    top+=ind_track_height;

     		    }
     		    this.dataRange.max = coverageMap.maximum;
                options.top+=top;
     		}
     	}


  

    class SCABAMTrack extends BAMTrack{
        constructor(config){
            if (!config.featureHeight){
                config.featureHeight=12;
            }
            super(config);
            this.display_alignments=false;
            this.cm_class= SCACoverageMap;
            this.keep_raw_alignments=true;
        }

        _setFeatureSource(){
            this.feature_source=new BamSource(this.config,this);

            this.feature_source.setViewAsPairs(true);

	        this.coverageTrack = new SCACoverageTrack(this.config, this);
            this.alignmentTrack = new AlignmentTrack(this.config, this);
        }

        setSelected(selected,count){
        	this.selected = selected;
        	let cm =this.feature_source.alignmentContainer.coverageMap;
        	this.selected_count=count;
        	cm.resetSelected();
            for (let alignment of this.feature_source.alignmentContainer.raw_alignments){
            	    cm.incCounts(alignment,true);
				}
        }

        getFeatures(chr, bpStart, bpEnd,force,data) {
			if (bpEnd- bpStart>200000){
				this.draw_zoom_in=true;
				return new Promise(function (fulfill, reject) {
					fulfill([]);
				});
			}
			this.draw_zoom_in=false;
			return super.getFeatures(chr, bpStart, bpEnd,force,data);
        }

    

        setCoverageCategory(cat_to_color,cat_counts,cat_func,update){
            this.categories=cat_to_color;
        	let ind_track_height= this.config.ind_track_height?this.config.ind_track_height:50;
        	let li=[];
        	for (let val in cat_counts){
        		li.push({val:val,num:cat_counts[val]});
        	}
        	li.sort(function(a,b){
        		if (a.val==="other"){
        			return 1;
        		}
        		return b.num-a.num;
        	});
        	this.cat_order=["selected"];
        	for (let l of li){
        		this.cat_order.push(l.val);
        	};
        	

        	this.config.height=(Object.keys(cat_to_color).length+1)*ind_track_height;
        	this.category_counts = cat_counts; 
        	this.category_function=cat_func;
        	if (update){
                let cm = this.feature_source.alignmentContainer.coverageMap;
                cm.reset();
                for (let alignment of this.feature_source.alignmentContainer.raw_alignments){
            	    cm.incCounts(alignment);
				}
        	}

     	}
     drawScale(pixel_height,ctx){
     	if (this.draw_zoom_in){
     		return;
     	}
     	this.coverageTrack.drawScale(ctx);
     
	
		}

	addExtraControls(dialog){
		let self = this;
		dialog.div.empty();
	
	

	
    	dialog.scale_slider=$("<div>").slider({
       		max:100,
       		min:10,
       		value:this.config.ind_track_height?this.config.ind_track_height:50,
       		slide:function(e,ui){
           		dialog.config.y_max=ui.value
           		if (dialog.panel){
                	dialog.panel.setTrackAttribute(dialog.config.track_id,"ind_track_height",ui.value);
                	let h = Object.keys(self.category_counts).length+1;
                	dialog.panel.setTrackAttribute(dialog.config.track_id,"height",h*ui.value);
                	dialog.panel.update();
           		}
           		
       		} 
    	}).appendTo(dialog.div);
    
				
	}





     

    }

     	MLVTrack.track_types["sca_bam_track"]={
            "class":SCABAMTrack	
     	}

export {SCABAMTrack};