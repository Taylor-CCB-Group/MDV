    	import {MLVTrack} from '../tracks.js'

    	class TiledImageTrack extends MLVTrack{
			constructor(config){
				super(config);
			}


			getFeatures(chr,bpStart,bpEnd,force,data){
				let locations = this.config.images[chr];
				if (!locations){
					return new Promise(function(fulfill,reject){
						fulfill([]);
					})
				}
				let promises= [];
				for (let loc of locations){
					if (loc[0]>bpEnd){
						break;
					}
					if (loc[1]<bpStart){
						continue;
					}
					promises.push(getImagePromise(loc));
				}
				return Promise.all(promises);
				
			}

			drawScale(pixelHeight,ctx){
 
				//ctx.drawImage(this.image,0,0,this.legend_width,pixelHeight,0,0,this.legend_width,pixelHeight)
			}

			drawFeatures(options) {
				let bpEnd = options.bpStart + options.pixelWidth * options.bpPerPixel + 1;
				let height=0;
				for (let image of options.features){
					let x1= (image.start-options.bpStart)/options.bpPerPixel;
					x1=x1<0?0:x1;
					let x_clip=0;
					if (options.bpStart>image.start){
						x_clip=((options.bpStart-image.start)/(image.end-image.start))*image.width;
 					}
					let x_clip2=image.width
					if (image.end>bpEnd){
						x_clip2 = image.width-(((image.end-bpEnd)/(image.end-image.start))*image.width);
					}
                           let x2 = options.pixelWidth;
					if (image.end<bpEnd){
						x2  = options.pixelWidth-((bpEnd-image.end)/options.bpPerPixel);
					}
					let factor=this.config.y_scale_factor?options.bpPerPixel/this.config.y_scale_factor:1;
					options.context.drawImage(image,x_clip,0,x_clip2-x_clip,image.height,x1,options.top,x2-x1,image.height/factor);
					height=image.height;
                         
					
				}
				
                this.bottom = options.top+height;
           		return this.bottom;
			}

			

		}
		function getImagePromise(loc){
			return new Promise(function (fulfill, reject) {
					let image = new Image();
					image.start=loc[0];
					image.end=loc[1];
					image.onload = function () {
    						fulfill(image);
					};
					image.src =loc[2]; 
			});       
		}
		

		MLVTrack.track_types['tiled_image_track']={
			"class":TiledImageTrack,
			name:"Tiled Image Track"

		}

export {TiledImageTrack};

