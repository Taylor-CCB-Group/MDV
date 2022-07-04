    	import {MLVTrack} from '../tracks.js'

    	class UCSCMLVTrack extends MLVTrack{
			constructor(config){
				config.url= config.url.replace("hgTracks","hgRenderTracks");
				if (config.url.includes("/s/")){
					config.url=UCSCMLVTrack._convertSessionURL(config.url);
				}
				super(config);
				this.image=null;
				this.legend_width=74;
			}

			static _convertSessionURL(url){
				let arr=url.split("/");
				let new_url = arr[0]+"//"+arr[2]+"/cgi-bin/hgRenderTracks?hgS_doOtherUser=submit&hgS_otherUserName="+arr[4]
					+"&hgS_otherUserSessionName="+arr[5];
				return new_url;

			}

			addExtraControls(dialog){
				let self = this;
				dialog.div.append("<hr>");
				dialog.div.append($("<label>").text("Session Url"));
				let inp=$("<textarea>").val(dialog.config.url).css({width:"100%",height:"70px"});
				let but = $("<i class='fas fa-download'></i>").click(function(e){
					let url = self._convertURL(inp.val());
					dialog.config.url=url;
					dialog.panel.setTrackAttribute(self.config.track_id,"url",url);
					dialog.panel.update();
				});
				dialog.div.append(but).append(inp);
			}

			_convertURL(url){
				url=url.replace("hgTracks","hgRenderTracks");
				return url
			}

			getFeatures(chr,bpStart,bpEnd,force,data){
				let self=this;
				let width = data.pixelWidth+this.legend_width;
				width = ((bpEnd-bpStart)/data.bpPerPixel)+this.legend_width;
				return new Promise(function (fulfill, reject) {
					self.image = new Image();
					self.image.onload = function () {
						self.config.height = self.image.height;
    					fulfill([]);
					};
					let url = self.config.url;
					let args ="&position="+chr+":"+bpStart+"-"+bpEnd+"&pix="+(width)+"&hgt.labelWidth=10";
					self.image.src =url+args;        
        		});		
			}

			drawFeatures(options) {
          		let ctx = options.context;
          		this.top=options.top;
          		ctx.globalAlpha=this.config.opacity;
           		ctx.drawImage(this.image,-this.legend_width,options.top);
           		ctx.globalAlpha=1.0;
           		this.bottom = options.top+this.image.height;
           		return this.bottom;
			}

			drawScale(pixelHeight,ctx){
				ctx.globalAlpha=0.7;
				ctx.drawImage(this.image,0,0,this.legend_width,pixelHeight,0,this.top,this.legend_width,pixelHeight);
				ctx.globalAlpha=1.0;
			}

		}

		MLVTrack.track_types['ucsc_track']={
			"class":UCSCMLVTrack,
			name:"UCSC Session"
		}

export {UCSCMLVTrack};
