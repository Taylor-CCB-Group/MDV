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

/**
 * Class representing a lightweight panel the can host multiple tracks of
 * different types
 */
import {MLVTrack,RulerTrack} from "./tracks.js";
import "./bam_track.js";
import "./BamCoverageTrack.js";
import {createEl} from "../utilities/Elements.js";
import SettingsDialog from "../utilities/SettingsDialog.js";
import canvasToSvg from "canvas-to-svg"

class MLVPanel {
	/**
	* Creates a panel
	* @param {array} tracks - a list of config objects describing each track
	* @param {object} config - A config with the panel settings
	*/
	constructor (div,config={},tracks=[]) {
		if (typeof div === "string"){
			this.trackDiv=document.getElementById("#"+div)
		}
		else{
			this.trackDiv = div;
		}
		const box = this.trackDiv.getBoundingClientRect();
	
		this.height = box.height
		this.width= box.width;
	
		this.tracks={}
		this.track_order=[];
		this.trackLabels={};
		this.trackDialogs={};
		for (let t_config of tracks){
			let track=MLVTrack.getTrack(t_config);
			this.tracks[track.config.track_id]=track;
			this.track_order.push(track.config.track_id);
			this._addTrackLabel(track.config)
		}
		//check for linked scales
		this._tracksChanged();
		this.legend= null;

      
       	this.trackDiv.style.position="relative";

		this.canvas = createEl("canvas",{
			classes:["igv-content-canvas"],
			width:this.width,
			height:this.height
		},this.trackDiv);

        this.ctx = this.canvas.getContext("2d");

	

        if (config.show_scale){
        	this.addScaleCanvas(this.height);
			this.show_scale=true;
        }

       
        //this.trackDiv.append(Utils.spinner());

        /*let icon_div=$("<div>").css({"z-index":100,position:"absolute",top:"2px",right:"2px"}).appendTo(this.trackDiv)
        					   .attr("class","panel-icon-div")

		*/
        //for event handlers
       	this.is_dragging=false;
       	this.isMouseDown = false,
      	this.lastMouseX = undefined;
       	this.mouseDownX = undefined;

		//amount to show each side of view port
       	this.buffer_level=1;

       	this.groups={};

       	this.highlighted_regions={};

       	//listeners
		this.listeners={};
       	

       	if (config.allow_user_interactions){
       		this.allowUserDrag();
			this.allowUserZoom();
			this.allowUserFeatureClick();
			this.allowUserFeatureOver();
			this.allowUserRangeSelection();
			this.trackDiv.addEventListener("dragover",e=>{
				e.preventDefault();
			})

       	}
      

       	if(config.ruler_track){
       		this.addRulerTrack();
       	}
       
       	if (config.new_layout){
       		this.new_layout=true;
       	}

		this.__doc__=document;
      
       
       	this.retries=0;
       	this.yOffset=0;

		this.loadingIcon=createEl("i",{
			classes:["fas","fa-spinner","fa-spin","panel-loading-icon"]
		},this.trackDiv);
		this.messageAlert=createEl("div",{
			classes:["panel-error-message"]
		},this.trackDiv)
    }

	closeAllDialogs(){
		for (let td in this.trackDialogs){
			this.trackDialogs[td].close();
		}
		this.trackDialogs={};
	}

	_addTrackLabel(config){
		if (config.show_track_label===false){
			//return;
		}
		const trackLabel =createEl("div",{
			classes:["panel-track-label"],
			text:config.short_label,
			draggable:true
		},this.trackDiv);
		trackLabel.addEventListener("click",e=>{
			this._openTrackDialog(e,config.track_id);
		
		});
		trackLabel.addEventListener("drag",e=>{
			this.trackDiv.style.cursor="ne-risize";

		});
		trackLabel.addEventListener("pointerdown",e=>{
			e.stopImmediatePropagation();

		});


		trackLabel.addEventListener("dragend",(e)=>{
			let y = this._getCoords(e).y +this.yOffset;
			//trackLabel.style.cursor="pointer"
			
			let tafter = null;
			if (y<0){
				tafter==this.track_order[0];
			}
			else{
				for (let tc in this.tracks){
					const t = this.tracks[tc];
					
					if (y>=t.top && y<=t.bottom){
						tafter =tc;
						break;
						
					}
				}
			}
			if (tafter===config.track_id){
				return;
			}
			tafter= tafter || this.track_order[this.track_order.length-1];
			this.track_order.splice(this.track_order.indexOf(config.track_id),1);
			this.track_order.splice(this.track_order.indexOf(tafter)+1,0,config.track_id);
			this.update();
			
		});
		
		this.trackLabels[config.track_id]=trackLabel;

	}

	_openTrackDialog(e,id){
		if (this.trackDialogs[id]){
			return;
		}
		const t =  this.tracks[id];
		this.trackDialogs[id]=  new SettingsDialog({
			maxHeight:400,
			doc:this.__doc__ || document,
			width:300,
			title:t.config.short_label,
			position:[e.pageX,e.pageY],
			onclose:()=>delete this.trackDialogs[id]

		},t.getSettings(this));
	}

    _parseConfig(config){
    	//check the tracks have the right settings
    	if (this.fixed_height_mode){
    		if (!config.height){
    			config.height=150;
    		}
    	}
    }

    addScaleCanvas(height){
    	this.scale_canvas = createEl("canvas",{
			styles:{
				position:"absolute",
				top:"0px",
				left:"5px"
			},
			width:100,
			height:this.height
		},this.trackDiv)
    	
        this.scale_ctx=this.scale_canvas.getContext("2d");
    }

    /**
	* sets the extra amount of track to draw each side of the view. A value 
	* of 1 will retreive 1 x the view width each side i.e. 3 x the visible window
	* @param {integer} level - The type of listener - track_empty 
	*/
    setBufferLevel(level){
    	this.buffer_level=level;
    }
    /**
	* Returns the element that houses the panel
	* @returns {integer} level - The type of listener - track_empty 
	*/

    getDiv(){
    	return this.trackDiv;
    }
    

        /**
	* Sets the highligted region
	* @param {Object} location - An object containing chr, start and end
	* @param {name} The name(id) of the region (used to remove the region)
	* @param {String} The color to give the highligted region
	*/
    setHighlightedRegion(location,name,color){
    	this.highlighted_regions[name]={
    		chr:location.chr,
    		start:location.start,
    		end:location.end,
    		color:color
    	}
    	this.force_redraw=true;
    }

    /**
	* Removes the highlighted region from the panel
	* @param {string} name - The name of the highlighted region
	* that was given when it was created.
	*/
    removeHighlightedRegion(name){
    	delete this.highlighted_regions[name];
    	this.force_redraw=true;
    }


    addRulerTrack(){
    	let track=new RulerTrack();
    	let config = track.getConfig();
		//config.show_label=false;
		this.tracks[config.track_id]=track;
		this.track_order.unshift(config.track_id);
		this._addTrackLabel(config);
		return this;
    }

    /**
	* Adds a listener to the panel
	* @param {string} type - The type of listener - track_empty
	* @param {function} func - The function to call 
	* @param {string} id - The id of the handler (can be used to remove the handler)
	* Optional - an id will be assigned (and returned) if not supplied
	* @returns{string} The id of the handler or null if the type did not exist 
	*/
    addListener(id,func){
    	this.listeners[id]=func;
    }

    /**
	* Removes a listener to the panel
	* @param {string} type - The type of listener - track_empty 
	* @param {string} id - The id of the handler to remove
	* @returns{boolean} true if the listener was removed, otherwise false 
	*/
    removeListener(id){
    	delete this.listeners[id];
    }
    
     
    /**
	* Removes a listener to the panel
	* @param {object} config - The config of the track to addTrack
	* @param {integer} index - Optional, the vertical order of the track
	*/
    addTrack(config,index,no_propagate){
    	let track=MLVTrack.getTrack(config);
    	
		this.tracks[track.config.track_id]=track;
		if (index || index==0){
			this.track_order.splice(index,0,track.config.track_id)
		}
		else{
			this.track_order.push(track.config.track_id);
		}
		this._tracksChanged();	
		this._addTrackLabel(track.config);
		if (!no_propagate){
			this._callListeners("track_added",track.config);
		}
    }
    
    _callListeners(type,data){
    	for (let id in this.listeners){
			this.listeners[id](type,data);
		}   
    }

    removeAllTracks(){
    	let dup_array = this.track_order.slice();
    	for (let id of dup_array){
    		this.removeTrack(id,true)
    	}
    }
    
	/**
	* Removes a listener to the panel
	* @param {object} config - The config of the track to add 
	*/
    removeTrack(track_id,not_repaint,not_propagate){
    	if (!this.tracks[track_id]){
    		return null;
    	}
    	this.track_order = this.track_order.filter(e => e !== track_id);
    	if (!not_repaint){
    		this.repaint(true,true);
    	}
    
    	if (this.legend){
    		this.legend.removeTrack(track_id);
    	}
    	let config =  this.tracks[track_id].config
    	delete this.tracks[track_id];
    	if (! not_propagate){
    		this._callListeners("track_removed",config);
		}
    	if (this.track_order.length===0){
            for (let l_id in this.listeners["panel_empty"]){
                this.listeners["panel_empty"][l_id](this);
            }
        }
        return config;
    }

    getTrackConfig(track_id){
    	let track = this.tracks[track_id];
    	return track.getConfig();
    }

	getTrack(track_id){
		return this.tracks[track_id];
	}

    getAllTrackConfigs(){
    	let configs=[];
    	for (let id of this.track_order){
			configs.push(this.tracks[id].getConfig());
    	}
    	return configs;
    }
    
    setTrackAttribute(track_id,key,value){
    	let track = this.tracks[track_id];
    	if (!track){
    		return;
    	}
    	track.setConfigAttribute(key,value);
    	if (key==="scale_link_to"){
    		this.tracks[track_id].scale_link_to = this.tracks[value];
    	}
    	if ((key==="color" || key==="display") && this.legend){
    		this.legend.updateTrack(track_id);
    	}
    }
    
    setTrackAttributes(track_id,attributes){
    	let track = this.tracks[track_id];
    	for (let key in attributes){
    		track.setConfigAttribute(key,attributes[key]);
    		if (key==="color" && this.legend){
        		this.legend.updateTrack(track_id);
        	}
    	}
    }


    /**
	* Sets the filter  function for track. 
	* @param {string} track_id- The id of the track
	* @param {string} func - The filter function. It should accept the feature
	* and return true to dispaly the feature and false to hide it. Use null 
	* to cancel the filter
	*/
    setTrackFeatureFilter(track_id,func){
    	let track = this.tracks[track_id];
    	track.setFilterFunction(func);
    }

    /**
	* Sets the filter  function for track 
	* @param {string} track_id- The id of the track
	* @param {string} func - The color function. It should accept the feature
	* and return the feature color. Use null to go back to default colors 
	*/
    setTrackColorFunction(track_id,func){
    	let track = this.tracks[track_id];
    	track.setColorFunction(func);
    }

    setTrackLabelFunction(track_id,func){
    	let track = this.tracks[track_id];
    	if (track){
    		track.label_function=func;
    	}
    }


    _tracksChanged(){
    	for (let t_id of this.track_order){
    		let track = this.tracks[t_id];
    		//if this track is linked to the scale of another
    		//get pointer to the track
    		let link_to = track.config['scale_link_to'];
    		if (link_to){
				let other_track = this.tracks[link_to];
				if (other_track){
					track.scale_link_to=other_track;
				}
    		}
    	}
    }


    getCurrentTrackFeatures(track_id){
    	let track = this.tracks[track_id];
    	return track.getCurrentFeatures(this.chr,this.start,this.end);
    }



   /**
    * Updated the panel view, if chromosome start and end are supplied
    * it will go to this location. If no parameters are given the panel
    * will be redrawn at the same location e.g after the color, scale or another
    * paramter has been set
    * @param {string} force - If true then a cached image will not be used
    * @param {integer} start of the region to draw
    * @param {integer} end of the region to draw
    */

    update (chr,start,end,no_propagation) {
    	this.call_update_listener=no_propagation;
		if (start>end){
			const temp =start;
			start=end;
			end=temp;
		}
        if (chr){
            this.chr=chr;
            this.start=start;
            this.end=end;
            this.repaint();
        }
        else{
        	this.repaint(true,true);
        }


       
    };

    getTracksHeight(){
    	let h =0;
    	let groups={}
    	for (let tid of this.track_order){
    		let track=  this.tracks[tid];
    		let g = track.config.group;
			if (track.config.hide){
				continue;
			}
    		if (g){
    			if (!groups[g]){
    				h+=track.config.height;
    				groups[g]=true;
    			}
    		}
    		else{
				h+=track.config.height;
    		}
    	}
    	return h;
    }


    async getAllFeatures(bpStart, bpEnd,force,data) {
        let promises = [];
		this._display_order=[];
        for (let track_id  of this.track_order){
        	let track = this.tracks[track_id];
			if (track.config.hide){
				continue;
			}
			this._display_order.push(track_id);
        	promises.push(track.getFeatures(this.chr,bpStart,bpEnd,force,data));       
        }
        return await Promise.all(promises.map(p => p.catch(e => e)));    
    }

    /**
     * Repaint the view, using a cached image if available.
     * @param {boolean} force - If true then a cached image will not be used
     * @param {boolean} range_from_tile - Redraw the tile
	 * @param {function} [svg] - if a callback is supplied , then the
	 * actual browser will not repaint but a serialises svg string of the 
	 * browser's contents will be passed to the function
      */
    repaint(force,range_from_tile,svg) {
		this.messageAlert.style.display="none";
        var pixelWidth,
            bpWidth,
            bpStart,
            bpEnd,
            self = this,
            ctx,
            chr,
            refFrameStart,
            refFrameEnd;

        chr = this.chr;
        refFrameStart = this.start;
        refFrameEnd = this.end;
        this.bpPerPixel=(this.end-this.start)/this.canvas.width;
        let get_features=true;
        if (this.tile && this.tile.containsRange(chr, refFrameStart, refFrameEnd, this.bpPerPixel)){
            get_features=false;
        } 
        if (!get_features && !force && !this.force_redraw) {
            this.force_redraw=false;
            this.paintImage();
            if (!self.call_update_listener){
            	self._callListeners("view_changed",{
					chr:self.chr,
					start:parseInt(self.start),
					end:parseInt(self.end)
				})
            }
            self.call_update_listener=false;
            self.retries=0;
        }
        else {
            // Expand the requested range so we can pan a bit without reloading
            this.force_redraw=false;
			if (!svg){
            	pixelWidth = ((this.buffer_level*2)+1) * this.canvas.width;
				bpStart = Math.max(0, Math.round(this.start-(this.buffer_level*this.canvas.width*this.bpPerPixel)));
			}
			else{
				pixelWidth= this.canvas.width;
				bpStart = Math.max(0, Math.round(this.start));
			}
            bpWidth = Math.round(pixelWidth*this.bpPerPixel);
           
            bpEnd = bpStart + bpWidth;
            if (self.loading && !(svg)){
            	if (force && range_from_tile){
            		self.update_required=true;
            	}
            	else{
            		self.update_required="location";
            	}
            	return;
            }
            if (range_from_tile){
            	if (this.tile){
                    bpStart=this.tile.startBP;
                    bpEnd=this.tile.endBP;
            	}
            }

         
			if (!svg){
				self.loading = {start: bpStart, end: bpEnd};
			
				self.loadingIcon.style.display="block";
			}
			const textColor = svg?"black":getComputedStyle(self.trackDiv).getPropertyValue('color');

            self.getAllFeatures( bpStart, bpEnd,!get_features,{pixelWidth:pixelWidth,bpPerPixel:self.bpPerPixel})

                .then(function (all_features) {
                    
                    if (all_features) {
                        
						let buffer = null;
						if (svg){
							//create mock canvas
							buffer= {};

						}	
						else{
							buffer = document.createElement('canvas');
						}
                        
                        buffer.width = pixelWidth;
                        buffer.height = self.getTracksHeight();
						if (svg){
							ctx= new canvasToSvg(buffer.width,buffer.height);
						}
						else{
							ctx = buffer.getContext('2d');
						}

                        
                        if (self.show_scale && !(svg)){
        					self.scale_buffer= document.createElement('canvas');
        					self.scale_buffer.width = 200;
        					self.scale_buffer.height = buffer.height;
        					self.scale_buffer_ctx=self.scale_buffer.getContext("2d");
                        }
                   
                        var options ={
                             context: ctx,
                             bpStart: bpStart,
                             bpPerPixel: self.bpPerPixel,
                             pixelWidth: buffer.width,
                             pixelHeight: buffer.height,
                             chr:chr,
							 textColor:textColor
                        };
                        let top=0;
                        self.groups={};
                        self.calculateMaxScale(all_features);
						for (let id in self.tracks){
							const c= self.tracks[id].config

							if (c.hide){
								self.trackLabels[c.track_id].style.display="none";
							}
						}
						
                        for (var i in all_features){
                        	let track = self.tracks[self._display_order[i]];
							
							const l = self.trackLabels[track.config.track_id];
							if (l){
								l.style.display="block";
								l.style.top=(top+self.yOffset)+"px";

							}
							
                        	options.features=all_features[i];
						
                        	let group = track.config.group
                        	if (group){
                        		if (!self.groups[group]){
                        			self.groups[group]={top:top,height:track.config.height}
                        			//first time increase top
                        			top+=track.config.height;
                        		}
                        		options.top= track.top=self.groups[group].top;
                        		options.height=self.groups[group].height;
								track.bottom=options.top+options.height;

                        	}
							else{
                        		options.top =top;
								track.top=top;
								options.height= track.config.height;
								top+=options.height;
								track.bottom=top;
                        	}
							if (typeof options.features === "string" && track.config.type!=="fasta" ){
								track.drawMessage(options);
								continue;
							}
                        	
							
							ctx.save();
							ctx.rect(0,options.top,options.pixelWidth,options.height);
							ctx.clip();
							ctx.beginPath();
							track.drawFeatures(options);
							ctx.restore()
                                   
                            if (self.show_scale){
								const sctx = svg?ctx:self.scale_buffer_ctx;
                            	track.drawScale(options.pixelHeight,sctx,textColor);
                            }
                                      
                        }
						
                        for (let name in self.highlighted_regions){
                        	let region = self.highlighted_regions[name];
                        	if (self.chr !== region.chr){
                        		continue;
                        	}
                        	if (region.end<bpStart ||region.start>bpEnd){
                        		continue
                        	}
                        	self.drawHighlightedRegion(region,options);
                        }
                        if (svg){
							svg(ctx.getSerializedSvg());
							return;
						}
                        self.retries=0;
                        self.loading = false;
                        self.tile = new Tile(chr, bpStart, bpEnd, self.bpPerPixel, buffer);
                        self.paintImage();
                        if (!self.call_update_listener){
                        	self._callListeners("view_changed",{
								chr:self.chr,
								start:parseInt(self.start),
								end:parseInt(self.end)
							})
                        }
                        self.call_update_listener=false;
                    
                    }
                    else {
                        self.ctx.clearRect(0, 0, self.canvas.width, self.canvas.height);
                    }
                    if (self.update_required){
                    	if (self.update_required==="location"){
                    		self.update(self.chr,self.start,self.end);
                    		self.update_required=false;
                    	}
                    	else{
                    		self.update_required=false;
                    		self.update();
                    	}
                    }	
                    self.loadingIcon.style.display="none";

                })
                .catch(function (error) {
                    self.loading = false;
                  
                    console.log(error);
                    if (self.retries<3 && error!=="Timed out"){
                    	self.retries++;
                    	self.repaint(force,range_from_tile);
                    }
                   
                    else{
                       
                        self.loading=false;
                        self.force_redraw=true;
						self.loadingIcon.style.display="none";
                        self.messageAlert.textContent=error;
						self.messageAlert.style.display="block";
                    }
                });
        } 
    }


    autoScale(features,min,max){
                if (!features || typeof features == "string" ){
                	return({min:0,max:1})
                }
        		features.forEach(function (f) {
            		min = Math.min(min, f.value);
           			max = Math.max(max, f.value);
        		});
        		return {min: min, max: max};
    		
    }

    calculateMaxScale(all_features){
    	  let groups={};
    	  for (var i in all_features){
              let track = this.tracks[this.track_order[i]];
              track.set_scale=null;
              let group =track.config.group;
             if (group && track.config.scale!=="fixed" && !(track.config.scale_link_to)){
             		track.config.scale_group=group;
             }       	
             group = track.config.scale_group
             if (group){
             	let group_info= groups[group];
             	if (!group_info){
             		group_info={tracks:[track],features:[all_features[i]]}
             		groups[group]=group_info

             		
             	}
             	else{
             		group_info.features.push(all_features[i]);
             		group_info.tracks.push(track);
             	}
             
             }
    	  }
    	  for (let name in groups){
    	  	let g= groups[name];
    	  	if (!g.ignore){
    	  		let min=0;
    	  		let max = -Number.MAX_VALUE;
    	  		let scale=null;
    	  		for (let f of g.features){
    	  			 scale= this.autoScale(f,min,max)
    	  			 min= scale.min;
    	  			 max=scale.max;
    	  		}
    	  		for (let t of g.tracks){
    	  			t.set_scale=scale;
    	  		}
    	  	}
    	  } 
    }

    drawHighlightedRegion(region,options){
    	let start= (region.start-options.bpStart)/options.bpPerPixel;
    	start = start<0?0:start;

    	let width = (region.end-region.start)/options.bpPerPixel;
    	width = width<3?3:width;
    	width =width>options.pixelWidth?options.pixelWidth:width;
		options.context.globalAlpha=0.1;
		options.context.fillStyle=region.color;
    	options.context.fillRect(start,0,width,options.pixelHeight);
    	options.context.globalAlpha=1.0;
    }
    

    paintImage() {

        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        if (this.show_scale){
        	this.scale_ctx.clearRect(0, 0,100, this.canvas.height);
        	this.scale_ctx.drawImage(this.scale_buffer,0,this.yOffset)
        }

        if (this.tile) {
            this.xOffset = Math.round((this.tile.startBP - this.start)/this.bpPerPixel);
            this.ctx.drawImage(this.tile.image, this.xOffset, this.yOffset);
            this.ctx.save();
            this.ctx.restore();
        }
    };

    allowUserFeatureOver(){
    	let self = this;
		if (this.foListener){
			return;
		}
		this.foListener= (e)=>{
			
				if (this.loading){
					return;
				}
				clearTimeout(this.moto);
				this.moto=setTimeout(()=>{
					if (!this.is_dragging){
					   let info=this.getFeatureAt(e);
					   let i = this.mouse_over_feature;
					   if (info.feature ){
						   if(i && i.feature!==info.feature){
							   this._callListeners("featureout",{track:i.track,feature:i.feature,event:e});
						   }
					   
						   if ((!i) || (i.feature!==info.feature)){
							this._callListeners("featureover",{track:info.track,feature:info.feature,event:e});
							   this.mouse_over_feature=info;
						   }
					   }
					   else{
						   let i = this.mouse_over_feature
						   if (i){
							this._callListeners("featureout",{track:i.track,feature:i.feature,event:e});
							   this.mouse_over_feature=null;
						   }
					   }
					}
				},10);
			}
			this.trackDiv.addEventListener("mousemove",this.foListener)
		
    }

    removeFeatureOverHandler(){
    		this.trackDiv.removeEventListener("mouseover",this.foListener);
			this.foListener=null;
    }


	allowUserFeatureClick(){
    	if (this.fcListener){
			return;
		}
    	this.allowUserFeatureOver();
		if (!this.listeners["_fcfo"]){
			this.addListener("_fcfo",(type,data)=>{
				if (type==="featureout"){
					this.trackDiv.style.cursor="default"
				}
				if(type==="featureover"){
					this.trackDiv.style.cursor="pointer";
				}
	
			})
		}
    
		this.fcListener= (e)=> {
    	 	if (this.loading){
    	 		return;
    	 	}
    	 	clearTimeout(this.foto);
    	 	this.foto=setTimeout(()=>{
    	 		if (!this.is_dragging){
					let i=this.getFeatureAt(e);
					if (i.track){
						this._callListeners("featureclick",{track:i.track,feature:i.feature,event:e});
					}
    	 		}
    	 	},200);
    	 };
		 this.trackDiv.addEventListener("click",this.fcListener)

    }


    removeFeatureClickHandler(){
    	this.trackDiv.removeEventListener("click",this.fcListener);
		this.fcListener=null;
		this.removeListener("_fcfo")
    }



    allowUserDrag(){
		this.drmdListener= e => {
    	 	if (e.shiftKey){
    	 		return;
    	 	}
            const co = this._getCoords(e)
            this.isMouseDown = true;
            this.start_dragging=true;
            this.lastMouseX = co.x;
            this.lastMouseY= co.y;
            this.mouseDownX = this.lastMouseX;

        }
		this.trackDiv.addEventListener("pointerdown",this.drmdListener)


       this.drmmListener= e=>{
            let canvasCoords = this._getCoords(e)
            if (this.is_dragging || this.start_dragging){
                var diff = canvasCoords.x-this.lastMouseX;
                var bp_diff=this.bpPerPixel*diff;
                this.start-=bp_diff;
                this.end-=bp_diff;
                let dd = this.canvas.height-this.tile.image.height;
                let y_diff=0;
                if (dd<0 || this.yOffset !==0){
                	let y_diff=  canvasCoords.y-this.lastMouseY;
                	this.yOffset+=y_diff;
                	if (this.yOffset>0){
                		this.yOffset=0;
                	}
                	else if (this.yOffset<dd){
						this.yOffset=dd;
                	}
					
					for (let id in this.tracks){
						const tr = this.tracks[id];
						const l= this.trackLabels[id];
						if (l){
							l.style.display="block";
							l.style.top=(tr.top+this.yOffset)+"px";

						}

					}
                }
                this.repaint();
                this.lastMouseX=canvasCoords.x;
                this.lastMouseY=canvasCoords.y;
                if (this.start_dragging && (diff>5 || y_diff>5)){
                	this.is_dragging=true;
                	this.start_dragging=false;
                }
               }
        }

		this.trackDiv.addEventListener("mousemove",this.drmmListener)

        this.drmuListener = e => {   
              this.is_dragging=false;
              this.start_dragging=false;
        }
		this.trackDiv.addEventListener("pointerup",this.drmuListener)
  
    }

    removeDragHandler(){
    	this.trackDiv.removeEventListener("pointerdown",this.drmdListener);
		this.trackDiv.removeEventListener("mousemove",this.drmmListener);
		this.trackDiv.removeEventListener("pointerup",this.drmuListener);
		this.drmdListener=null;
		this.drmmListener=null;
		this.drmuListener=null;
    }

    _getCoords(e){
		const box = this.canvas.getBoundingClientRect();
    	const x = e.clientX - Math.round(box.left);
        const y = e.clientY - Math.round(box.top);
        return {x,y};   	
    }


    allowUserZoom(){
		this.zmListener= e =>{
		
			e.preventDefault();
			let deltaY= e.deltaY;
			if (deltaY === undefined){
				deltaY=e.detail
			}
    	 	if (this.loading || (this.bpPerPixel<0.05 && deltaY<0)){
    	 		return;
    	 	}
    	 	
    	 	let canvasCoords = this._getCoords(e);
            let factor = deltaY>0?2:0.5;
            let mbp= (this.start+ canvasCoords.x * this.bpPerPixel)
            let new_length = (this.end-this.start)*factor;
            let new_start = mbp-((canvasCoords.x/this.canvas.width)*new_length);
            this.start=  new_start
            this.end= new_start+new_length;        
            this.repaint();      
		}
    	this.trackDiv.addEventListener('wheel',this.zmListener);     
    }

    disableUserZoom(){
		this.trackDiv.removeEventListener("wheel",this.zmListener);
    }
    
    
    allowUserRangeSelection(){
		this.rsmdListener= e=>{
            if (e.shiftKey){
                this.start_select =this._getCoords(e).x;
                let left = this.start_select+"px";		
                this.select_div=createEl("div",{
					styles:{
						"position":"absolute","opacity":0.2,"background-color":"blue","top":"0px","height":this.height+"px",left:left,"width":"0px"
					}
				},this.trackDiv)
                e.stopPropagation();
            }
        } 
        this.trackDiv.addEventListener("pointerdown",this.rsmdListener);

		this.rsmmListener= e =>{
            if (e.shiftKey && this.start_select){
                let x=this._getCoords(e).x;
                if (x<this.start_select){
                    this.select_div.style.left=x+"px";
					this.select_div.style.width=(this.start_select-x)+"px";
                }
                else{                  
                     this.select_div.style.left=this.start_select+"px";
					 this.select_div.style.width=(x-this.start_select)+"px";
                }           
                e.stopPropagation();
            }
		}
        this.trackDiv.addEventListener("mousemove",this.rsmmListener);

		this.rsmuListener= e =>{
			if (this.start_select){
                let x=this._getCoords(e).x;
                let start = this.start + (this.start_select*this.bpPerPixel);
                let end =  this.start + (x*this.bpPerPixel);
                this.start_select=null;
                this.select_div.remove();
                if (start>end){
                	let temp=end;
                	end=start;
                	start=temp;
                }
                this._callListeners("range_selected",{
						chr:this.chr,
						start:start,
						end:end
				});
            }
		}
        this.trackDiv.addEventListener("pointerup",this.rsmuListener);
    }

    
    removeAllowSelection(){
    	this.trackDiv.removeEventListener("pointerdown",this.rsmdListener);
		this.trackDiv.removeEventListener("mousemove",this.rsmmListener);
		this.trackDiv.removeEventListener("pointerup",this.rsmuListener);
		this.rsmdListener=null;
		this.rsmmListener=null;
		this.rsmuListener=null;
    }

    getImage(){
    	 var imgURL = this.canvas[0].toDataURL(MIME_TYPE);
    }
    
    
    /**
	* Gets the feature that was clicked
	* @param {JQuery Event} e - Can be any object- all that is required is pageX and PageY
	* @returns {object} An object with track - the track config at the event position(or null) and
	* feature - the feature at the postition (or null). 
	*/
     getFeatureAt(e){
    	 let co = this._getCoords(e);
    	 co.y-=this.yOffset;
    	 let gl = Math.round(this.start+(co.x*this.bpPerPixel));
    	 for (let t in this.tracks){
    	 	let track = this.tracks[t];
    	 	if (co.y>track.top && co.y<track.bottom){
    	 		return {track:track,
    	 				feature:track.getFeatureAt(gl,this.chr,co,this.bpPerPixel,this.ctx,this.yOffset)
    	 		};
    	 	}		
    	 }
    	 return {track:null,feature:null};
    }
    



	setSize(width,height){
		if (!width){
			const b = this.trackDiv.getBoundingClientRect();
			this.width=b.width;
			this.height=b.height;
		}
		else{
			this.height=height;
			this.width=width;
			this.trackDiv.style.height= height+"px";
			this.trackDiv.style.width=width+"px";
		}
		this.canvas.setAttribute("width",this.width);
		this.canvas.setAttribute("height",this.height);
		if (this.show_scale){
        	this.scale_canvas.setAttribute("height",this.height);
        }

	}

  

    
    redrawTile(features) {

        if (!this.tile) return;

        var self = this,
            chr = self.tile.chr,
            bpStart = self.tile.startBP,
            bpEnd = self.tile.endBP,
            buffer = document.createElement('canvas'),
            bpPerPixel = self.tile.scale;

        buffer.width = self.tile.image.width;
        buffer.height = self.tile.image.height;
        var ctx = buffer.getContext('2d');

      

        self.track.draw({
            features: features,
            context: ctx,
            bpStart: bpStart,
            bpPerPixel: bpPerPixel,
            pixelWidth: buffer.width,
            pixelHeight: buffer.height
        });


        self.tile = new Tile(chr, bpStart, bpEnd, bpPerPixel, buffer);
        self.paintImage();
    }

}


class Tile{
	constructor (chr, tileStart, tileEnd, scale, image) {
		this.chr = chr;
		this.startBP = tileStart;
		this.endBP = tileEnd;
		this.scale = scale;
		this.image = image;
	}

	containsRange(chr, start, end, scale) {
		if (start<0){
			start=0;
		}
		return this.scale.toFixed(3) === scale.toFixed(3) && start >= this.startBP && end <= this.endBP && chr === this.chr;
	}

	overlapsRange(chr, start, end) {
		return this.chr === chr && this.endBP >= start && this.startBP <= end;
	}
}

export default MLVPanel;