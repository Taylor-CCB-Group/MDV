
import regl from  "regl";
import * as glMatrix from "gl-matrix";
import {createEl, makeDraggable} from "../utilities/Elements.js";
import {Camera} from "./Camera.js";


class WGL3DI{
	/** 
    * Creates a new wg2di instance
    * @param {string|object} div - The id of the div to house the instance or jquery element
    * @param {integer} width - The width of the instance (optional- will be the width of the parent div if not supplied)
    * @param {integer} height - The height  of the instance (optional- will be the height of the parent div if not supplied)
    * @param {object} config (optional)
    * <ul>
    * <li> in_view_only - if true then only those objects in view will be drawn after panning/zooming.This will
    * speed things up if maniluplating individual objects but slow down panning and zooming</li>
    * <li> allow_object_drag - if true individual objects can be dragged </li>
    * <li> default_brush - if true then brushing will be the default drag behavior - uses shift for panning
    *  Otherwise panning is default and shift used for brushing
	*
    * </ul>
    */
	constructor(div,width,height,config){
		width = Math.round(width);
		height= Math.round(height);
		var self = this;
		if (!config){
			config={};
		}
		this.config=config;
		if (config.circle_borders==null){
			this.circle_borders=true;
		}
		else{
				this.circle_borders=config.circle_borders;
		}


		this.draw_options=config.draw_options?config.draw_options:{depth:{},blend:{}};
		
	

    
    	this.regl=null;
		this.pick_buffer=null;
		this.objects=[];
		this.keys={};
		this.object_types=[];
		this.objects_in_view=0;
        this.camera = new Camera({distance:500,theta:0.75,phi:0.5});
		//html elements
		if (typeof div === "string"){
			this.div_container=document.getElementById(div);
		}
		else{
			this.div_container=div;
		}
		this.canvas=null;
		this.label_context;
		this.label_context=null;
		this.height=0;
		this.width=0;
		this._setUpDocument(width,height);
		this.images_to_load=0;
		this.is_filtered={};
		this.is_hidden={};
    
		//handlers
		this.handlers={
			object_clicked:{},
			object_over:{},
			object_out:{},
			brush_stopped:{},
			zoom_stopped:{},
			panning_stopped:{},
		};

		//switches
		this.draw_labels=false;

		//circle shapes
		this.circle_properties={x_pos:1,y_pos:1,z_pos:1,color:1,pick_color:1,localFilter:1,globalFilter:1};
		this.circles_to_draw={};
		this.circles={};
		for (var prop in this.circle_properties){
			this.circles[prop]=[];
		}
		this.circles.count=0;
		this.object_types.push({data:this.circles,
								properties:this.circle_properties,
								vertices:1,
							primitive:"points"});
		this.universal_circle_radius=20;

		//line shapes
		this.line_properties={'position':2,"color":2,"pick_color":2};
		this.lines_to_draw={};
		this.lines={};
		for (var prop in this.line_properties){
			var list=[];
			this.lines[prop]=list;   
			this.lines_to_draw[prop]=list;
		}
		this.lines.count=0;
		this.object_types.push({data:this.lines,
								data_in_view:this.lines_to_draw,
								properties:this.line_properties,
								vertices:2,
								primitive:"lines"});

		//rectangles
		this.rect_properties={'position':2,"color":2,"pick_color":2,"opacity":2};
		this.rects_to_draw={};
		this.rects={};
		for (var prop in this.rect_properties){     
			this.rects[prop]=[];   
			this.rects_to_draw[prop]=[];
		}
		this.rects.count=0;

		this.object_types.push({data:this.rects,
								data_in_view:this.rects_to_draw,
								properties:this.rect_properties,
								vertices:6,
								primitive:"triangles"});


		//squares
		this.square_properties={'position':1,"color":1,"pick_color":1,"size":1,"opacity":1};
		this.squares_to_draw={};
		this.squares={};
		for (var prop in this.square_properties){     
			this.squares[prop]=[];   
			this.squares_to_draw[prop]=[];
		}
		this.squares.count=0;
		this.object_types.push({data:this.squares,
								data_in_view:this.squares_to_draw,
								properties:this.square_properties,
								vertices:1,
								primitive:"points"});


		this.scale=1.0;
		this.x_scale=1.0;
		this.y_scale=1.0;
		this.offset=[0,0];


		//images
		this.image_properties={'position':2,"pick_color":2,"color":2,"opacity":2};
		this.images_to_draw={};
		this.images={};
		for (var prop in this.image_properties){     
			this.images[prop]=[];   
			this.images_to_draw[prop]=[];
		}
		
		this.images.count=0;
		this.images.props=[];
		this.images.globals=[this.offset[0],this.offset[1],this.x_scale,this.y_scale];
		this.object_types.push({data:this.images,
								data_in_view:this.images_to_draw,
								properties:this.image_properties,
								vertices:6,
								primitive:"triangles"});
	    this.images.display_as_image=true;

        

        //image tiles
	    this.imagetile_properties={'position':2,"pick_color":2,"color":2,"opacity":2,"x_y":2,"c_r":2};
		this.imagetiles_to_draw={};
		this.imagetiles={};
		for (var prop in this.imagetile_properties){     
			this.imagetiles[prop]=[];   
			this.imagetiles_to_draw[prop]=[];
		}
		
		this.imagetiles.count=0;
		this.object_types.push({data:this.imagetiles,
								data_in_view:this.imagetiles_to_draw,
								properties:this.imagetile_properties,
								vertices:6,
								primitive:"triangles"});




		//The last mouse position recorded
		this.mouse_position=null;
		//Was an object clicked
		this.object_clicked=null;
		//an object was clicked
		this.dragging=false;
		//object which mouse is over
		this.object_mouse_over=null;
		this.mouse_over_color=null;//[255,0,0];


		this.zoom_amount=0;
		if (!this.config.draw_order){
			this.config.draw_order=[0,1,2,3,4,5];
		}



		

		regl({
			onDone: function(err,regl){
				self.regl=regl;
				regl._refresh()
				self.pickbuffer = regl.framebuffer({ colorFormat: 'rgba',height:self.height,width:self.width});
				self._initDrawMethods();
				self._addHandlers();
				 
			},
			canvas:self.canvas,
			attributes:{
				antialias:false,
			}
			
			

		});

	}





	removeAllObjects(){
		for(let ob of this.object_types){
			for (let prop in ob.data){
				if (prop==="count"){
					ob.data[prop]=0;
					continue;
				}
				ob.data[prop]=[];
			}
			for (let prop in ob.data_in_view){
				if (prop==="count"){
					ob.data_in_view[prop]=0;
					continue;
				}
				ob.data_in_view[prop]=[];
			}
		}

		this.objects=[];
		this.keys={};
		
	}


	_setUpDocument(width,height){
		if (this.config.brush){
			this.div_container.style.cursor="crosshair";
		}
		if (!height){
            const box= this.div_container.getBoundingClientRect();
			this.height=height=box.height;
			this.width =width=box.width;

		}else{

			this.height=height;
			this.width=width;
			this.div_container.style.height=height+"px";
			this.div_container.style.width = width+"px";
		}

        const attr = {
			height:height,
			width:width,
            styles:{
                //position:"absolute",
                left:"0px",
                top:"0px"
            }
        }
		this.canvas =createEl("canvas",attr);
	    this.label_canvas=createEl("canvas",attr);
		this.label_context=this.label_canvas.getContext("2d");
		this.div_container.append(this.canvas);
        this.div_container.append(this.label_canvas);
	}

	resizeImage(i){
		this.imagetiles.image_scale=i;
		this.maintainImageDimensions();
	}


	setImageSize(height,width){
		this.images.default_size=[height,width];
		this.images.scale_factor=1;
	}




	maintainImageDimensions(){
		if (this.imagetiles.count===0){
			return;
		}
		let width = this.imagetiles.width*this.imagetiles.image_scale;
		let height = this.imagetiles.height*this.imagetiles.image_scale;
		let t_height= height*this.y_scale;
		let t_width = width*this.x_scale;
		if (this.x_scale>this.y_scale){
			width= (this.imagetiles.ratio*height*this.y_scale)/this.x_scale;
		}
		if (this.y_scale>this.x_scale){
			height = (this.imagetiles.ratio*width*this.x_scale)/this.y_scale;
		}
		let pos = this.imagetiles.position;
		let max= pos.length-6
		for (let n=0;n<=max;n+=6){
			let x= this.imagetiles.x_y[n][0];
			let y = this.imagetiles.x_y[n][1];
			pos[n+1][1]=y+height;
			pos[n+2][0]=x+width;
			pos[n+2][1]=y+height;
			pos[n+3][0]=x+width;
			pos[n+3][1]=y+height;
			pos[n+4][0]=x+width;
		}
		
		this.imagetiles.i_w_h = [width,height];
		this.imagetiles.w_h=[width*this.imagetiles.sprite_dim[0],height*this.imagetiles.sprite_dim[1]];
	}



	remove(){
		this.pickbuffer.destroy();
		this.regl.destroy();
		this.canvas.attr({height:1,width:1});
		this.canvas.remove();

	}
	

	setSize(width,height){
		if (this.height===height && this.width===width){
			return
		}
		let self=this;
		width=Math.round(width);
		height=Math.round(height);
		let x_ratio=width/this.width;
	    let y_ratio=height/this.height;
	  
		this.x_scale= this.x_scale*x_ratio;
		this.y_scale= this.y_scale*y_ratio;
		this.maintainImageDimensions();
	
		this.height=height;
		this.width=width;
		this.div_container.style.height=height+"px";
		this.div_container.style.width =width+"px"
		this.canvas.setAttribute("height",height);
		this.canvas.setAttribute("width",width);
		this.label_canvas.setAttribute("height",height);
		this.label_canvas.setAttribute("width",width);	
		
		this.pickbuffer.destroy()
		this.pickbuffer = this.regl.framebuffer({ colorFormat: 'rgba',height:height,width:width});
		//this is necessary, but I  don't know why?
		let loop =this.regl.frame(function(){
			
		});
		setTimeout(()=>{

			self.refresh();
			loop.cancel();

			},100);
		


	}

	_getMousePosition(e){
    	var rect = this.canvas.getBoundingClientRect();
    	return [e.clientX-rect.left,e.clientY-rect.top];
	}

	_getActualPosition(position){
    	var x = (position[0]/this.x_scale) - this.offset[0];
    	var y = (position[1]/this.y_scale) - this.offset[1];
    	return [x,y];
	}

	getRange(){
		let tl = this._getActualPosition([0,0]);
		let br = this._getActualPosition([this.width,this.height]);
		return {
			x_range:[tl[0],br[0]],
			y_range:[tl[1],br[1]],
			offset:[this.offset[0],this.offset[1]],
			scale:[this.x_scale,this.y_scale]

		};
	}

	highlightPoint(key){
		clearInterval(this.an_func);
		let obj = this.keys[key];
		if (!obj){
			return;
		}
		let pt=this.circles.position[this.objects[obj][0]];
		this.label_context.clearRect(0, 0, this.width, this.height);
	    this.highlight_point=pt
		pt = this._getCanvasCoords(pt);
		let ctx= this.label_context;
		let self = this;
		self.highlight_pt_radius=0.1;
		ctx.fillStyle="orange";
		self.an_func=setInterval(function(){
			ctx.beginPath();
            ctx.arc(pt[0], pt[1], self.highlight_pt_radius, 0, 2 * Math.PI);
            ctx.fill();
            self.highlight_pt_radius+=0.1;
            if (self.highlight_pt_radius>(self.universal_circle_radius*self.x_scale)){

            	clearInterval(self.an_func);
            }
		},5);
		

		
	}
	removeHighlightPoint(){
		if (!this.keys["__highlight__"]){
			return;
		}
		let end=this.circles.count-1;
		for (let n in this.circle_properties){
			this.circles[n].splice(end,1)
		}
		this.objects.splice(this.objects.length-1,1);
		delete this.keys["__highlight__"];
		this.ciecles.count--;
		

	}

	_getCanvasCoords(pos){
		let x=(pos[0]+this.offset[0])*this.x_scale;
        let y=(pos[1]+this.offset[1])*this.y_scale;
        return [x,y]

	}

	_drawLabels(){
    	var time =Date.now();
    	if (this.objects_in_view>5000){
        	return;
    	}
    	this.label_context.font = "14px Arial";
   
    	for(var i in this.object_types){
        	var pos =this.object_types[i].data_in_view.position;
        
        	for (var ii=0;ii<pos.length;ii++){
       
          		var x=(pos[ii][0]+this.offset[0])*this.x_scale;
          		var y=(pos[ii][1]+this.offset[1])*this.y_scale;
                
             	this.label_context.fillText("H1",x,y);
           
        	}
    	}    
	}

	setObjectColor(key,color){
		var obj = this.objects[this.keys[key]];
		if (!obj){
			return;
		}
		var obj_type= this.object_types[obj[1]];
		for (var x=obj[0];x<obj[0]+obj_type.vertices;x++){
			obj_type.data.color[x][0]=color[0]/255;
			obj_type.data.color[x][1]=color[1]/255;
			obj_type.data.color[x][2]=color[2]/255;
		}
	}

	setObjectOpacity(key,opacity){
		var obj = this.objects[this.keys[key]];
		if (!obj){
			return;
		}
		var obj_type= this.object_types[obj[1]];
		let arr = obj_type.data.props[obj[0]].opacity
		for (var x=0;x<6;x++){
		
			arr[x]=opacity;
			arr[x]=opacity;
			arr[x]=opacity;
		}
	}

	setObjectPosition(key,x,y){
		var obj = this.objects[this.keys[key]];
		if (!obj){
			return;
		}
		var obj_type= this.object_types[obj[1]];
		for (var i=obj[0];i<obj[0]+obj_type.vertices;i++){
			obj_type.data.position[i][0]=x		
			obj_type.data.position[i][1]=y;
		}
			
	}		



	setUniversalCircleRadius(value){
		if (!value || isNaN(value)){
			value=0;
		}
		this.universal_circle_radius=value;
	}

	getUniversalCircleRadius(){
		return this.universal_circle_radius;
	}


	getObjectColor(key){
		var obj = this.objects[this.keys[key]];
		var obj_type= this.object_types[obj[1]];
		var col= obj_type.data.color[obj[0]];
		return [col[0]*255,col[1]*255,col[2]*255];

	}

	setLogPosition(axis,bool,type){
		let neg_vals=false;
		if (!type){
			neg_vals=true;
		}
		let i=0
		if (axis ==="y"){
			i=1;
		}
		if (bool){
			this.temp_vals={}	
		}
	
		for (let obj of this.objects){
			let obj_type=this.object_types[obj[1]];
			let arr =  obj_type.data.position[obj[0]];
			if (arr[i]!==0){
				let abs = Math.abs(arr[i]);
			    let m=arr[i]<0?-1:1;
				if (bool){
					if ((neg_vals && abs<1) || (!(neg_vals) && m===-1)){
						arr.temp=arr[i];
						arr[i]=0;
					}
                    else{
					    arr[i]= Math.log10(abs)*m;
					}
				}
				else{
					if ((neg_vals && abs<1) || (!(neg_vals) && m===-1)){
						arr[i]=arr.temp;
						delete arr.temp
					}else{
					    arr[i]=Math.pow(10,abs)*m;
				    }
			     }
		    }
		}
	
	}




	addLine(positionTo,positionFrom,color){
    	var index = this.objects.length;
    
     
	
		this.lines.position.push(positionFrom);
		this.lines.position.push(positionTo);
		this.lines.color.push([color[0]/255,color[1]/255,color[2]/255]);
		this.lines.color.push([color[0]/255,color[1]/255,color[2]/255]);
		this.lines.pick_color.push(this._getRGBFromIndex(index+1));
		this.lines.pick_color.push(this._getRGBFromIndex(index+1));
		
		this.lines.count++;
	}

	addThickLine(positionTo,positionFrom,width,color,key){
		var index = this.objects.length;
		if (key && ! this.keys[key]){
			this.keys[key]=index;
		}
		else{
			key=index;
			this.keys[index]=index;
		}
		var rect_index=this.rects.position.length;
		var x_diff= positionTo[0]-positionFrom[0];
		var y_diff = positionTo[1]-positionFrom[1];
		var factor = (0.5*width)/Math.sqrt((x_diff*x_diff)+(y_diff*y_diff));
		var x_offset= factor*y_diff;
		var y_offset= factor*x_diff;
		this.rects.position.push([positionTo[0]+x_offset,positionTo[1]-y_offset]); //TL
		this.rects.position.push([positionFrom[0]+x_offset,positionFrom[1]-y_offset]); //BL
		this.rects.position.push([positionFrom[0]-x_offset,positionFrom[1]+y_offset]); //BR

		this.rects.position.push([positionFrom[0]-x_offset,positionFrom[1]+y_offset]); //BR
		this.rects.position.push([positionTo[0]-x_offset,positionTo[1]+y_offset]); //TR
		this.rects.position.push([positionTo[0]+x_offset,positionTo[1]-y_offset]); //TL


		var c  = [color[0]/255,color[1]/255,color[2]/255];
		var pc = this._getRGBFromIndex(index+1);
		for (var a=0;a<6;a++){
			 this.rects.color.push(c);
			 this.rects.pick_color.push(pc);
		}
		this.objects.push([rect_index,2,key]);
		this.rects.count++;
		return key;
	}


	addImageTile(sprite_sheets,config,positions,callback){
		let height = config.image_height;
		let width = config.image_width;
	    this.imagetiles.width=width;
	    this.imagetiles.height=height;
	    this.imagetiles.original_dim=[width,height];
	    this.imagetiles.ratio = width/height;
	    this.imagetiles.image_scale=1;
	    
		this.imagetiles.w_h=[width*config.sprite_dim[0],height*config.sprite_dim[1]];
		this.imagetiles.i_w_h=[width,height];

		this.imagetiles.sprite_dim=config.sprite_dim
	
		for (let pos of positions){
			var index = this.objects.length;
			
			this.keys[pos.key]=index;
		
			let x = pos.x;
			let y = pos.y;
			
			let imagetile_index=this.imagetiles.position.length;
			this.imagetiles.position.push([x,y]);
			this.imagetiles.position.push([x,y+height]);
			this.imagetiles.position.push([x+width,y+height]),
			this.imagetiles.position.push([x+width,y+height]);
			this.imagetiles.position.push([x+width,y]);
			this.imagetiles.position.push([x,y]);
		
			let pc = this._getRGBFromIndex(index+1);
			let x_y = [x,y];
			let c_r = [pos.sheet,pos.col,pos.row]
			for (var a=0;a<6;a++){
				this.imagetiles.pick_color.push(pc);
			 	this.imagetiles.color.push([1,1,1]);
			 	this.imagetiles.opacity.push(1.0);
			 	this.imagetiles.x_y.push(x_y);
			 	this.imagetiles.c_r.push(c_r);
			}
			this.imagetiles.count++;


			this.objects.push([imagetile_index,5,pos.key]);

		}
	
       
        let promises=[];
        let self =this;
        for (let n=0;n<5;n++){
        	if (sprite_sheets[n]){
				promises.push(this.loadSpriteSheet(sprite_sheets[n],n))
        	}
        	else{
        		self.imagetiles["text"+n]=this.regl.texture({shape: [1, 1]})

        	}
        }
        Promise.all(promises).then(function(){
        	callback();
        })
       

	}


	loadSpriteSheet(url,index){
		let self =this;
		return new Promise(function(fulfill){
			var im = new Image();
			im.onload=function(){

             	self.imagetiles["text"+index]=self.regl.texture({data:im,min:"linear"});
              	fulfill()
        	}

        	im.src = url;
        
		})
	}


	addImage(position,height,width,image,key,opacity){
		let self =this;
		var index = this.objects.length;
		if (key && ! this.keys[key]){
			this.keys[key]=index;
		}
		else{
			key=index;
			this.keys[index]=index;
		}
		let x = position[0];
		let y = position[1];
		var image_index=this.images.position.length;
		this.images.position.push([x,y]);
		this.images.position.push([x,y+height]);
		this.images.position.push([x+width,y+height]),
		this.images.position.push([x+width,y+height]);
		this.images.position.push([x+width,y]);
		this.images.position.push([x,y]);
		this.image_position=[[-1,1],[-1,0],[0,0],[0,0],[0,1],[-1.0]];
		var pc = this._getRGBFromIndex(index+1);
		opacity=opacity==null?1:opacity;
		for (var a=0;a<6;a++){
			 this.images.pick_color.push(pc);
			 this.images.color.push([1,1,1]);
			 this.images.opacity.push(opacity);
		}

		let i_index=this.images.props.length;

		this.images.props.push({count:6,is_buffer:0,x_y:[x,y],w_h:[width,height],opacity:this.images.opacity.slice(image_index,image_index+6),
		text:this.loading_image,position:this.images.position.slice(image_index,image_index+6),globals:this.images.globals,color:this.images.color.slice(image_index,image_index+6)});
		
		var im = new Image();
        im.src = image;
        im.width
        this.images_to_load++;
        im.onload=function(){
              self.images.props[i_index].text=self.regl.texture({data:im,min:"linear"});
              self.images_to_load--;
              self.refresh();
        }
        this.images.count++;


		this.objects.push([image_index,4,key,i_index]);


	}

	addRectangle(position,height,width,color,key){
		var index = this.objects.length;
		if (key && ! this.keys[key]){
			this.keys[key]=index;
		}
		else{
			key=index;
			this.keys[index]=index;
		}
		var rect_index=this.rects.position.length;
		this.rects.position.push(position);
		this.rects.position.push([position[0],position[1]+height]);
		this.rects.position.push([position[0]+width,position[1]+height]);

		this.rects.position.push([position[0]+width,position[1]+height]);
		this.rects.position.push([position[0]+width,position[1]]);
		this.rects.position.push([position[0],position[1]]);
		for (var a=0;a<6;a++){
			 this.rects.opacity.push(1);
		}
       

		var c  = [color[0]/255,color[1]/255,color[2]/255];
		var pc = this._getRGBFromIndex(index+1);
		for (var a=0;a<6;a++){
			 this.rects.color.push(c);
			 this.rects.pick_color.push(pc);
		}
		this.objects.push([rect_index,2,key]);
		this.rects.count++;
		return key;
	}	

	addArc(position,radius,color,start_angle,end_angle,key){
		var index = this.objects.length;
		 if (key && ! this.keys[key]){
			this.keys[key]=index;
		}
		else{
			key=index;
			this.keys[index]=index;
		}
		var circ_index=this.circles.position.length;
		this.circles.position.push(position);
		this.circles.radius.push(radius);
		this.circles.color.push([color[0]/255,color[1]/255,color[2]/255]);
		this.circles.pick_color.push(this._getRGBFromIndex(index+1));
		this.circles.start_angle.push(start_angle);
		this.circles.end_angle.push(end_angle);
		this.objects.push([circ_index,0,key]);
		this.circles.count++;
	}

    colorPoints(colorFunc){
        const len =this.circles.count;
        const colors = this.circles.color;
        for (let n=0;n<len;n++){
            const col = colorFunc(n);
            const cp =n*3;
            colors[cp]=col[0];
            colors[cp+1]=col[1];
            colors[cp+2]=col[2];
        }
    }


    addCircles(config){
        this.circles.x_pos=config.x;
        this.circles.y_pos=config.y;
        this.circles.z_pos=config.z;
        const len = config.x.length;
        this.circles.localFilter= config.localFilter;
		this.circles.globalFilter= config.globalFilter;
        this.circles.pick_color= new Uint8Array(len*3);
        this.circles.color=new Uint8Array(len*3);
        for (let n=0;n<len;n++){
            let p = n*3;
            this.circles.color[p]=40;
            this.circles.color[p+1]=220;
            this.circles.color[p+2]=80;
            let pb= this._getRGBFromIndex(n+1);
            this.circles.pick_color[p]=pb[0];
            this.circles.pick_color[p+1]=pb[1];
            this.circles.pick_color[p+2]=pb[2]
        }
        this.circles.count=len;
    }

	setGlobalFilter(filter){
		this.circles.globalFilter=filter;
	}


	addCircle(position,radius,color,key,opacity){
		var index = this.objects.length;
		 if (key && ! this.keys[key]){
			this.keys[key]=index;
		}
		else{
			key=index;
			this.keys[index]=index;
		}
		var circ_index=this.circles.position.length;
		this.circles.position.push(position);
		this.circles.radius.push(radius);
		this.circles.opacity.push(opacity?opacity:1);
		this.circles.color.push([color[0]/255,color[1]/255,color[2]/255]);
		this.circles.pick_color.push(this._getRGBFromIndex(index+1));
		this.circles.start_angle.push(10);
		this.circles.end_angle.push(0);
		this.objects.push([circ_index,0,key]);
		this.circles.count++;
	}

	removeCircle(key_list){
		let temp_li={};
		for (let n in this.circle_properties){
			temp_li[n]=[]
		}
        
        let keys={};
        for (let key of key_list){
        	keys[key]=true
        }
        let key_to_index={};
		for (let obj of this.objects){
			if (!keys[obj[2]]){
				key_to_index[obj[0]]=obj
			}

		}	

		let count=0;
		let temp_objs=[];
		for (let i =0;i<this.circles.count;i++){
			let obj= key_to_index[i];
			if (obj){
			    for (let n in this.circle_properties){
				    temp_li[n].push(this.circles[n][i]);
			    }
			    obj[0]=count;
			    count++
			    temp_objs.push(obj)

		    }
			   
		 }
		  for (let n in this.circle_properties){
		        this.circles[n]=temp_li[n];
		  }
		  this.circles.count=count;
		  setTimeout(()=>this.refresh(),0);
		  



	}

	addSquare(position,color,size,key){
		
		var index = this.objects.length;
		if (key && ! this.keys[key]){
			this.keys[key]=index;
		}
		else{
			key=index;
			this.keys[index]=index;
		}

		var square_index=this.squares.position.length;
		this.squares.position.push(position);
		this.squares.color.push([color[0]/255,color[1]/255,color[2]/255]);
		this.squares.size.push(size);
		this.squares.opacity.push(1);
		this.squares.pick_color.push(this._getRGBFromIndex(index+1));
		this.objects.push([square_index,3,key]);
		this.squares.count++;
	}




	_getRGBFromIndex(index){
		var b = Math.floor(index/65536);
		var temp = index%65536;
		var g= Math.floor(temp/256);
		var r = temp%256;
		return [r/255,g/255,b/255];
            
	}
	_getIndexFromRGB(rgb){
    	return (rgb[2]*65536)+(rgb[1]*256)+rgb[0];    
	}

	_drawPickBuffer(in_view){
       	this.regl.clear({
        	color: [0, 0, 0, 0],
			depth: 1,
        	framebuffer:this.pickbuffer
    	});
     	this._drawObjects(true,in_view);
    
	}
	//refesh all 
	//in_view only those in view
	refresh(in_view){
    	//this.label_context.clearRect(0, 0, this.width, this.height);
    	this._drawObjects(false,in_view);
    	this._drawPickBuffer(in_view);
    	this.label_context.font = "30px Arial";
   
	}

	zoom(amount){
    	this.x_scale*=amount;
    	this.y_scale*=amount;
   		this._drawObjects(false);
	}

	_drawObject(object,color){
    
		var type =this.object_types[object[1]];
		var obj={
			stage_width:this.width,
			stage_height:this.height,
			x_scale:this.x_scale,
			y_scale:this.y_scale,
			buffer:null,
			offset:this.offset,
			count:type.vertices,
			primitive:type.primitive,
			is_buffer:0,
			universal_radius:this.universal_circle_radius,
			circle_borders:this.circle_borders?1:0

		};
	

		for (var prop in type.properties){   
			if (prop==='pick_color'){
				  continue;
			}
			obj[prop]=[];
			var st= object[0];
			var en =st+type.vertices;
			for (var pos=st;pos<en;pos++){
				if (prop!=='color'){
					obj[prop].push(type.data[prop][pos]);
				}
				else{
					obj[prop].push(color);
				}
			}

		}

		type.method(obj);
	}


	_drawObjects(buffer,in_view){
		//dr
	
		this.images.globals[0]=this.offset[0];
		this.images.globals[1]=this.offset[1]
		this.images.globals[2]=this.x_scale;
		this.images.globals[3]=this.y_scale;
	
		let data_source = "data";
        const cProj = this.camera.getProjection();
		
		for (let i of this.config.draw_order){
			var type =this.object_types[i];
			if (!type[data_source].count){
				continue;
			}

			if (buffer){
				buffer=this.pickbuffer;
			}
			else{
				buffer=null;
			}
			//don't draw images unless buffer
			/*if (!(buffer) && i==4 && this.images.display_as_image){
				continue
			}*/

			var obj={
			   stage_width:this.width,
			   stage_height:this.height,
               cameraProjection:cProj.projection,
               cameraView:cProj.view,
			   x_scale:this.x_scale,
			   y_scale:this.y_scale,
			   globals:this.images.globals,
			   buffer:buffer,
			   offset:this.offset,
			   count:type[data_source].count * type.vertices,
			   primitive:type.primitive,
			   is_buffer:buffer?1:0,
			   universal_radius:this.universal_circle_radius,
			   circle_borders:this.circle_borders?1:0,
               cameraDistance:cProj.distance,
			   cameraEye:cProj.eye


			};
			//dummy values
			if (i==4){
				if (buffer){
					obj.w_h=[300,300];
					obj.x_y=[200,200];
					obj.text=this.loading_image;
					obj.is_buffer=1;
				}
				else{
					obj= this.images.props;
				}
				
				
			

			}
			if (i==5){
				obj.w_h=this.imagetiles.w_h;
				obj.text0=this.imagetiles.text0;
				obj.text1=this.imagetiles.text1;
				obj.text2=this.imagetiles.text2;
				obj.text3=this.imagetiles.text3;
				obj.text4=this.imagetiles.text4;
				obj.i_w_h=this.imagetiles.i_w_h;
			}
			for (var prop in type.properties){

				if (buffer){
					if (prop==='pick_color'){
						obj['color']=type[data_source][prop];
						continue;
					}
					if (prop==='color'){
						continue;
					}
					obj[prop]=type[data_source][prop];
				}
				else{
					if (prop==='pick_color'){
						continue;
					}
					obj[prop]=type[data_source][prop];
				}
			}
			type.method(obj);
		}
	}


	orderColumns(id_to_order,pos_to_id,rows,cols){
		  let max=rows*cols;
            let obj_type= this.object_types[2];
            let rp = this.rects.position;
            let opa = this.rects.opacity;
            let x_pos=0;
			for (let n=0;n<max;n+=rows){
			
                let id = pos_to_id[n];
                let x_pos = id_to_order[id]*2
			
				let st = n*6;
				let en = st+(rows*6);
			
				for (let i=st;i<en;i+=6){

					rp[i][0]=x_pos;
					rp[i+1][0]=x_pos;
					rp[i+2][0]=x_pos+2;
					rp[i+3][0]=x_pos+2;
					rp[i+4][0]=x_pos+2;
					rp[i+5][0]=x_pos;
					opa[i]=op;
					opa[i+1]=op;
					opa[i+2]=op;
					opa[i+4]=op;
					opa[i+5]=op;
		
				
				}
			}
			
	}



	moveAndHide(ids, id_to_pos, field_order){
          
            let obj_type= this.object_types[2];
            let rp = this.rects.position;
            let opa = this.rects.opacity;
            let x_pos=0;
			for (let info of id_to_pos){

				let op =0;
				if (ids[info[0]] !== undefined){
                    op=1;
				}
                let st = info[1]*6;
                let en = info[2]*6;
				
			
				for (let i=st;i<en;i+=6){

					rp[i][0]=x_pos;
					rp[i+1][0]=x_pos;
					rp[i+2][0]=x_pos+2;
					rp[i+3][0]=x_pos+2;
					rp[i+4][0]=x_pos+2;
					rp[i+5][0]=x_pos;
					opa[i]=op;
					opa[i+1]=op;
					opa[i+2]=op;
					opa[i+3]=op;
					opa[i+4]=op;
					opa[i+5]=op;
		
				
				}
				if (op===1){
					x_pos+=2;
				}
			}
			
		}


	makeOpaque(index){
		let op = this.circles.opacity;
		for (let n=0;n<op.length;n++){
			if (index[n]){
				op[n]=1;
			}
			else{
				op[n]=0.5;
			}
		}
			
		
	}


	filterImages(keys){
		let obj_type= this.object_types[4];
		for(let obj of this.objects){	
			let key =obj[2];
			if (keys[key]){
				this.images.props[obj[3]].opacity=[1.0,1.0,1.0,1.0,1.0,1.0]
			}
			else{
				//this.is_filtered[key]=true;
				this.images.props[obj[3]].opacity=[0.6,0.6,0.6,0.6,0.6,0.6];
			}
			if (this.is_hidden[key]){
				this.images.props[obj[3]].opacity=[0,0,0,0,0,0];
			}
			
		}
	}

	hideImages(keys){
		this.is_hidden={};
		let obj_type= this.object_types[4];

		for(let obj of this.objects){	
			let key =obj[2];
			if (keys[key]){
					this.images.props[obj[3]].opacity=[1.0,1.0,1.0,1.0,1.0,1.0]
				
			}
			else{
				this.is_hidden[key]=true;
				this.images.props[obj[3]].opacity=[0,0,0,0,0,0]
				
			}
		
			
		}

	}

	hideObjects(keys,object_type,def_op){
		if (! object_type){
			object_type=0;
		}
		this.is_hidden={};
		if (!def_op){
			def_op=1;
		}
		
		
		let obj_type= this.object_types[object_type];
		let vert=obj_type.vertices;
		for(let obj of this.objects){	
			let key =obj[2];
	        let op=def_op;
			let st = obj[0];
			if (!keys[key]){
				this.is_hidden[key]=true;
				op=0;
				
			}
			for (let n=0;n<vert;n++){
				obj_type.data.opacity[st+n]=op;
			}
			
		
			
		}
	}



	initialise(func,url){
		let self = this;

		var im = new Image()
       
        im.onload=function(){
              self.loading_image=self.regl.texture({data:im,min:"linear"});
              func();
        }
        im.src =url?url:"loading.png";
       
    
	}

	_getObjectAtPosition(position){
		if(position[1]===0){
			return;
		}
		var pixel = this.regl.read({
			x: position[0],
			y: this.height - position[1],
			width: 1,
			height: 1,
			data: new Uint8Array(6),
			framebuffer: this.pickbuffer
		});
		var index = this._getIndexFromRGB(pixel);
		if (index>0){
			return this.objects[index-1];
		}
		return null;
	}

	checkImagesLoaded(){
		let self = this;
		
		if (this.images_to_load>0){
			setTimeout(function(){
				self.refresh();
				self.checkImagesLoaded();
			},2000)
		}
	}

	_getObjectsInView(){
		if (!this.config.in_view_only){
			return;
		}
		var time = Date.now();
		var max = this.width*this.height*4;
		var pixels = this.regl.read({
			x: 0,
			y: 0,
			width:this.width,
			height: this.height,
			data: new Uint8Array(max),
			framebuffer: this.pickbuffer
		});
		var obj={};
		this._clearObjectsInView();
		for (var i=0;i<max-4;i+=4){
			var  index = pixels[i+2]*65536+pixels[i+1]*256+pixels[i];
			if (index>0){
				if(!obj[index-1]){
					obj[index-1]=true;
					this.objects_in_view++;
					if (this.objects_in_view>100000){
						for (var t in this.object_types){
							var type = this.object_types[t];
							for (var prop in type.properties){     
								type.data_in_view[prop]=(type.data[prop]);       
							}
							type.data_in_view.count=type.data.count;
						}

						this.objects_in_view=this.objects.length;
						return;

					}
				}      
			}       
		}
		 //console.log("objects in view old way "+(Date.now()-time));
		var l =  -this.offset[0];
		var r = l+(this.width/this.x_scale);
		var t =  -this.offset[1];
		var b = t+(this.height/this.y_scale);
		var old_count=0;
		var new_count=0;
		for (var i=0;i<this.objects.length;i++){
			if (obj[i]){
				var item= this.objects[i];
				var type =this.object_types[item[1]];
				new_count++;

				var st= item[0];
				var en =st+type.vertices;
				for (var prop in type.properties){    
					for (var pos=st;pos<en;pos++){
						type.data_in_view[prop].push(type.data[prop][pos]);       
					}

				}
				type.data_in_view.count++;

			}
			else{

				var item= this.objects[i];
				var type =this.object_types[item[1]];
				var act_pos =this.object_types[item[1]].data.position[item[0]];
				if (act_pos[0]>l && act_pos[0]<r && act_pos[1] >t && act_pos[1]<b){
					 old_count++;
					var st= item[0];
					var en =st+type.vertices;
					for (var prop in type.properties){    
						for (var pos=st;pos<en;pos++){
							type.data_in_view[prop].push(type.data[prop][pos]);       
						}

					}

					type.data_in_view.count++;
				}
			}





		}
		
		
	}

	_clearObjectsInView(){
		this.objects_in_view=0;
		for (var i in this.object_types){
			var obj = this.object_types[i];

				for (var prop in obj.properties){
					obj.data_in_view[prop]=[];
				}     

			obj.data_in_view.count=0;
		}  
	}

	addHandler(handler_type,handler,name){
		var handler_dict = this.handlers[handler_type];
		if (!handler_dict){
			throw "Handler Not Supported";
		}
		if (!name){
			name = Object.keys(handler_dict).length;
		}
		handler_dict[name]=handler;
		return name;
	}

	removeHandler(handler_type,name){
		var handler_dict = this.handlers[handler_type];
		if (!handler_dict){
			throw "Handler Not Supported";
		}
		delete handler_dict['name'];


	}

	_setUpBrush(origin){
		let self = this;
		let div =createEl("div",{
			classes:["wgl2di-brush"],
			styles:{
				top:origin[1]+"px",
				left:origin[0]+"px"
			}
		},this.div_container);

		
								  
		makeDraggable(div,{
			contain:true,
			onstart:()=>this.brush_moving=true,
			onstop:()=>this._brushingStopped()

		});
	
		/*.resizable({
			handles:"all",
			start:function(ev,ui){
				self.brush.moving=true;
			},
			stop:function (ev,ui){
				self._brushingStopped();

			}

		});*/
		this.brush={origin:origin,div:div,resizing:true};
	}

	_setUpPolyBrush(pos){
		this.poly_brush={
			points:[pos],
			active:true,
		}
		let ctx= this.label_context;
		
		ctx.beginPath()
		ctx.moveTo(pos[0],pos[1]);
	


	}

	_extendPolyBrush(pos,end){
		let ctx= this.label_context;
	    //let prev = this.poly_brush.points[this.poly_brush.points.length-1]
		
		
		ctx.lineTo(pos[0],pos[1]);
		ctx.stroke()
		this.poly_brush.points.push(pos);
		if (end){
			ctx.closePath();
			ctx.fillStyle="lightgray";
			ctx.globalAlpha=0.4;
			ctx.fill();
			let poly = []
			for (let pt of this.poly_brush.points){
				poly.push(this._getActualPosition(pt));
			}
			for (var i in this.handlers.brush_stopped){
			    this.handlers.brush_stopped[i](poly,true);
		    }

		}
	}

	_finishPolyBrush(pos){
		if (this.poly_brush.points.length<4){
			this.clearBrush();
			return;
		}
		let ctx= this.label_context;	
		ctx.closePath();
		ctx.fillStyle="lightgray";
		ctx.globalAlpha=0.2;
		ctx.fill();
		ctx.globalAlpha=1;
		let poly = [];
		this.poly_brush.active=false;
		for (let pt of this.poly_brush.points){
			poly.push(this._getActualPosition(pt));
		}
		for (var i in this.handlers.brush_stopped){
			this.handlers.brush_stopped[i](poly,true);
		}
	}



	clearBrush(){
		if (this.brush){
			this.brush.div.remove();
			this.brush=null;
		}
		if (this.poly_brush){
		    this.label_context.clearRect(0, 0, this.width, this.height);
		    this.poly_brush=null;
		}
	}

	_brushingStopped(){
		this.brush.resizing=false;
		const d= this.brush.div;
		const y= d.offsetTop;
		const x = d.offsetLeft;
		const w= d.offsetWidth;
		const h = d.offsetHeight;
		if (w<2){
			this.clearBrush();
			return;
		}
		
		let lt =this._getActualPosition([x,y]);
		let br = this._getActualPosition([x+w,y+h]);
		let info = {x_min:lt[0],x_max:br[0],y_max:-lt[1],y_min:-br[1]};
		for (var i in this.handlers.brush_stopped){
			this.handlers.brush_stopped[i](info);
		}

	}


	_addHandlers(){
		var self=this;
		this.div_container.addEventListener("mousemove",function(e){
            if (self.panning){
                var pos =self._getMousePosition(e);
                var dx = (pos[0] - self.panning[0]) / 300;
                var dy = (pos[1] - self.panning[1]) / 300;

                self.camera.mouseChange(dx,dy)
                self.panning[0]=pos[0];
                self.panning[1]=pos[1];
            }
			//no drag event going on call any listners if mouse over/out an object
			else if (false){
				var position =self._getMousePosition(e);
				var obj = self._getObjectAtPosition(position);
				if (obj && !self.object_mouse_over){
					for (var i in self.handlers['object_over']){
						self.handlers.object_over[i](obj[2],e);                  
					}
					self.object_mouse_over=obj;
					if (self.mouse_over_color){
						self.object_temp_color=self.getObjectColor(obj[2]);
						self.setObjectColor(obj[2],self.mouse_over_color);
						self.refresh(true);
					}
				}
				else if (!obj && self.object_mouse_over){
					for (var i in self.handlers['object_out']){
						self.handlers.object_out[i](self.object_mouse_over[2]);
					}
					if (self.mouse_over_color){
						self.setObjectColor(self.object_mouse_over[2],self.object_temp_color);
						self.refresh(true);
					}
					self.object_mouse_over=null;

				}
				//move directly from one object to another
				else if(obj && (obj[2]!==self.object_mouse_over[2])){
					for (var i in self.handlers['object_over']){    
						self.handlers.object_over[i](obj[2],e);  
					}


					/*for (var i in self.handlers['object_out']){
						self.handlers.object_out[i](self.object_mouse_over[2]);
					}*/
					if (self.mouse_over_color){
						self.setObjectColor(self.object_mouse_over[2],self.object_temp_color);
						self.object_temp_color=self.getObjectColor(obj[2]);
						self.setObjectColor(obj[2],self.mouse_over_color);
						self.refresh(true);
					}
					self.object_mouse_over=obj;





				}         
			}
		});

		this.div_container.addEventListener("mouseleave",function(evt){
            self.panning=null;
			if (self.object_mouse_over){
				for (var i in self.handlers['object_out']){
							self.handlers.object_out[i](self.object_mouse_over[2]);
						}
						self.object_mouse_over=null;
			}
			if (self.brush && self.brush.resizing){
			    self._brushingStopped();
			}
		})

	
		this.div_container.addEventListener("mouseup",function(evt){
            self.panning = null;
            self.loop.cancel();
            self.loop =null;
			//just a click event - inform handlers
			if (self.config.brush){
				self.div_container.style.cursor="crosshair";
			}
			if (self.brush && self.brush.resizing){
				self._brushingStopped();
				return;
			}
			if (self.poly_brush && self.poly_brush.active){
				self._finishPolyBrush();
			}
			if (!self.dragging){
				//if (self.object_clicked){
					var position =self._getMousePosition(evt);
					var obj = self._getObjectAtPosition(position);
					if (obj){
						for (var i in self.handlers.object_clicked){
							self.handlers.object_clicked[i](obj[2]);
						}  
					}
				//}

			}        
			else{
				//an object has finshed its drag
				if (self.object_clicked){
					self.object_clicked=null;
					self.refresh(true);              
				}
				//update which objects are now in view
				else{
					let ret = self.getRange();
					for (var i in self.handlers.panning_stopped){
							self.handlers.panning_stopped[i](ret);
					}
					self.loop.cancel();
					self.loop=null;
					  

					self._drawPickBuffer();
					self._getObjectsInView();
					if (self.brush){
						self._brushingStopped();
					}
				   self.refresh(true);
				}
				self.dragging=false;
			}
			self.object_clicked=null;
			self.mouse_position=null;   
		});  

		this.div_container.addEventListener('mousewheel', function(event){
			let z= 0;
		if (event.wheelDelta > 0 || event.detail < 0) {
				z=20;

		}
		else {
				z=-20;

			}
          self.camera.mouseWheel(0,z)
          self.refresh();

		});
		this.div_container.addEventListener("mousedown",function (evt){
			if (evt.which===3){
				//add right click behaviour
			}
            const pos=  self._getMousePosition(evt);
            self.panning= pos;
          
            self.loop = self.regl.frame(function(){
				self.regl.clear({
					color: [0, 0, 0, 0],
					depth: 1,
				});
				self._drawObjects(false);
			});
		

		});

	}


	_initDrawMethods(){
		var self = this;


		//loading images
		


		this.__drawCircles = this.regl({
            frag: `
            precision mediump float;
            varying vec3 fragColor;
            void main () {
              if (length(gl_PointCoord.xy - 0.4) > 0.4) {
                discard;
              }
              gl_FragColor = vec4(fragColor,1.0);
            }
            `,
          
            vert: `
            precision mediump float;

            attribute float x_pos;
            attribute float y_pos;
            attribute float z_pos;
            attribute vec3 color;
            attribute float lFilter;
            attribute float gFilter;

			uniform float time;
            uniform float cdistance;
            uniform mat4 view, projection;
			uniform vec3 eye;

            varying vec3 fragColor;

            void main () {
                if (gFilter>0.0 && gFilter != lFilter){
					return;
				}
                vec3 position = vec3(x_pos,y_pos,z_pos);
                fragColor=color/255.0;
                //float ps = 30.0/(pow(distance,2.0));
                //gl_PointSize =ps<30.0?ps:30.0;
				gl_PointSize =5.0;//- (distance(eye, position.xyz) / 10.0);
                gl_Position = projection * view * vec4(2.0 * position, 1.0);
            }
            `,
          
              attributes: {
              
                x_pos: self.regl.prop('x_pos'),
                y_pos:self.regl.prop("y_pos"),
                z_pos:self.regl.prop("z_pos"),
				color: self.regl.prop('color'),
				lFilter:self.regl.prop("localFilter"),
				gFilter:self.regl.prop("globalFilter")
              },
          
              uniforms: {
                view: self.regl.prop("cameraView"),
                cdistance:self.regl.prop("cameraDistance"),
				eye:self.regl.prop("cameraEye"),
                projection: (context,props) =>{
                 return glMatrix.mat4.perspective(
                    props.cameraProjection,
                    Math.PI / 4,
                    context.viewportWidth / context.viewportHeight,
                    0.01,
                    10000
                  )
                 },
                
              },
          
              count:  self.regl.prop('count'),
          
              primitive: 'points',
		
		});
		this.object_types[0]['method']=this.__drawCircles;

        this.__drawLines = this.regl({

            // fragment shader
            frag: ' precision highp float;\n\
                    varying vec3 fragColor;\n\
                    void main () {\n\
                         gl_FragColor = vec4(fragColor,1);\n\
                    }\n',


            vert: `
                    attribute vec3 position;
                    attribute vec3 color;
                    uniform mat4 view, projection;
    
                    varying vec3 fragColor;
                  
                    
                    void main () {
                     
                        fragColor=color;
                        gl_Position = projection * view * vec4(2.0 * position, 1.0);
                       
                       
                    }`
            ,
            attributes: {
                position: self.regl.prop("position"),
                color:self.regl.prop("color")
            


            },
            uniforms: {
                view: self.regl.prop("cameraView"),
                projection: (context,props) =>{
                 return glMatrix.mat4.perspective(
                    props.cameraProjection,
                    Math.PI / 4,
                    context.viewportWidth / context.viewportHeight,
                    0.01,
                    10000
                  )
                 },
                
              },

          
            primitive:self.regl.prop("primitive"),
            framebuffer:self.regl.prop("buffer"),
            count:self.regl.prop("count")



        });
    this.object_types[1]['method']=this.__drawLines;
    this.object_types[2]['method']=this.__drawLines;


    }	
	
}


export {WGL3DI};
