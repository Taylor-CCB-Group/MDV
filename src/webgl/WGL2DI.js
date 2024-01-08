import regl from  "regl";
import * as glMatrix from "gl-matrix";
import {createEl, makeDraggable} from "../utilities/Elements.js";
import {Camera} from "./Camera.js";



class WGL2DI{
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
	constructor(div,config){
	
		var self = this;
		if (!config){
			config={};
		}
		this.config=config;
		if (!config.offset){
			config.offset=[0,0,0,0];
		}
		this.mode = config.mode || "2d";
		this.movingImage=0;
		this.__doc__=document;

		this.config.brush=this.config.brush || "default";

	/*	this.draw_options=config.draw_options?config.draw_options:{depth:{enable:true},
		blend:{}
	}
	*/

		this.draw_options = config.draw_options ? config.draw_options : 
		{
			depth:{enable:false},
			blend:{
				enable:true,
				func: {
					srcRGB: 'src alpha',
					srcAlpha: 'src alpha',
					dstRGB: 'one minus src alpha',
					dstAlpha: 'one minus src alpha'
				}
			}
		}

    
    	this.regl = null;
		this.pickbuffer = null;

		this.hideOnFilter=true;


		this.pointOpacity=0.8;
		this.pointRadius=10;
		this.isFiltered=false;

		
		this.x_scale=1.0;
		this.y_scale=1.0;
		this.offset=[0,0];

	

		//html elements
		if (typeof div === "string"){
			this.div_container=document.getElementById(div);
		}
		else{
			this.div_container=div;
		}
		

		this._setUpDocument(config.width,config.height);
	
	
		//handlers
		this.handlers={
			object_clicked:{},
			object_over:{},
			object_out:{},
			brush_stopped:{},
			zoom_stopped:{},
			panning_stopped:{},
			pan_or_zoom:{}
		};

		//switches
		this.draw_labels=false;

		this.object_types=[];
		this.pointScale=0;

		// circle shapes
		this.circle_properties={x_pos:1,y_pos:1,color:1,pick_color:1,localFilter:1,globalFilter:1};
	
		this.circles={};
		for (var prop in this.circle_properties){
			this.circles[prop]=[];
		}
		this.circles.count=0;
		this.object_types.push({data:this.circles,
								properties:this.circle_properties,
								vertices:1,
								primitive:"points"
							});

	
		// line shapes
		this.line_properties={position:2,color:2,opacity:2};
	
		this.lines={};
		for (var prop in this.line_properties){
			this.lines[prop]=[];   ;
		}
		this.lines.count=0;
		this.object_types.push({data:this.lines,
								properties:this.line_properties,
								vertices:2,
								primitive:"lines"
							});

		//images
		this.image_properties={'position':2,"color":2,"opacity":2};
		this.images={};
		for (var prop in this.image_properties){     
			this.images[prop]=[];   
		
		}	
		this.images.count=0;


		//special propertis for images
		this.images.props={x_y:[],w_h:[],text:[]};
		
		this.object_types.push({data:this.images,
								properties:this.image_properties,
								vertices:6,
								primitive:"triangles"});

        
	   
		if (this.mode==="3d"){
			this.circle_properties.z_pos=1;
			this.camera = new Camera({
				center:config.cameraCenter,
				distance:config.cameraDistance || 500,
				theta:0.75,
				phi:0.5});
			this.axisScales=[1,1,1];
		}


		//squares
		this.square_properties={x_pos:1,y_pos:1,color:1,pick_color:1,localFilter:1,globalFilter:1};

		this.squares={};
		for (var prop in this.square_properties){
			this.squares[prop]=[];
		}
		this.squares.count=0;
		this.object_types.push({data:this.squares,
								properties:this.square_properties,
								vertices:1,
								primitive:"points"
							});
        
        this.rectangle_properties={position:6,color:6,opacity:6};
	
		this.rects={};
		for (var prop in this.rectangle_properties){
			this.rects[prop]=[];   ;
		}
		this.rects.count=0;
		this.object_types.push({data:this.rects,
								properties:this.rectangle_properties,
								vertices:6,
								primitive:"triangles"
							});

		//The last mouse position recorded
		this.mouse_position=null;
		//Was an object clicked
		this.object_clicked=null;
		//an object was clicked
		this.dragging=false;
		//object which mouse is over
		this.object_mouse_over=null;


		this.zoom_amount=0;
		this.draw_order= this.mode==="3d"?[2,0,1,3,4]:[2,0,1,3,4];
	

		regl({
			onDone: function(err,regl){
				self.regl=regl;
				regl._refresh()
				self.pickbuffer = regl.framebuffer({ colorFormat: 'rgba',height:self.height,width:self.width});
				self._initDrawMethods();
				if (self.mode==="3d"){
					self.object_types[0]['method']=self.__draw3DCircles;
					self.object_types[1]['method']=self.__draw3DLines;
				}
		
				self._addHandlers();
				 
			},
			canvas:self.canvas,
			attributes:{
				antialias:false,
			}
		});

	}

	setBackGroundColor(color){
		color = color || "";
		this.div_container.style.background=color;
	}



	_setUpDocument(width,height){
		if (this.config.brush){
			this.div_container.style.cursor="crosshair";
		}
		if (!height){
            const box= this.div_container.getBoundingClientRect();
			this.height=height=Math.round(box.height);
			this.width =width=Math.round(box.width);

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
                position:"absolute",
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
		//color highlights
		if (this.highlightPoints){
			const ps = this.highlightPoints;
			for (let i=0;i<ps.length;i++){
				const ci = ps.indexes[i]*3;
				const hi= i*3;
				ps.color[hi]=this.circles.color[ci];
				ps.color[hi+1]=this.circles.color[ci+1];
				ps.color[hi+2]=this.circles.color[ci+2];	
			}
		}
    }


	remove(){
		this.pickbuffer.destroy();
		this.regl.destroy();
		//this.canvas.attr({height:1,width:1});
		this.canvas.remove();
		this.label_canvas.remove();

	}
	

	setSize(width,height,rescale=null){
		if (this.height===height && this.width===width){
			return
		}
		let self=this;
		width=Math.round(width);
		height=Math.round(height);
		if (rescale){
			const scale= rescale(width,height);
			this.x_scale=scale[0] || this.x_scale;
			this.y_scale=scale[1] ||  thisy_scale
			
		}
		else{

			let x_ratio=width/this.width;
			let y_ratio=height/this.height;
		  
			this.x_scale= this.x_scale*x_ratio;
			this.y_scale= this.y_scale*y_ratio;
		
		}

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
			},200);
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

	
	_getCanvasCoords(pos){
		let x=(pos[0]+this.offset[0])*this.x_scale;
        let y=(pos[1]+this.offset[1])*this.y_scale;
        return [x,y]

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

	setHighlightPoints(indexes){
		if (indexes==null){
			this.highlightPoints=null;
			return;
		}
		const hlp={
			x_pos:[],
			y_pos:[],
			color:[],
			length:indexes.length,
			indexes:indexes
		}
		if (this.mode==="3d"){
			hlp.z_pos=[];
		}
		for (let index of indexes){
			hlp.x_pos.push(this.circles.x_pos[index]);
			hlp.y_pos.push(this.circles.y_pos[index]);
			if (this.mode==="3d"){
				hlp.z_pos.push(this.circles.z_pos[index]);
			}
		
			const st = index*3;
			hlp.color.push(this.circles.color[st]);
			hlp.color.push(this.circles.color[st+1]);
			hlp.color.push(this.circles.color[st+2]);
		}
		this.highlightPoints=hlp;		
	}


	setPointRadius(value){
		if (!value || isNaN(value)){
			value=0;
		}
		this.pointRadius=value;
	}
	setFilter(val){
		this.isFiltered=val;
	}

	setHideOnFilter(value){
		this.hideOnFilter=value;
	}

	setPointOpacity(value){
		value=value>1.0?1.0:value;
		value = value<0.0?0.0:value;
		this.pointOpacity=value
	}

	getPointOpacity(){
		return this.pointOpacity;
	}

	getPointRadius(){
		return this.pointRadius;
	}

	addLine(positionTo,positionFrom,color=[0,0,0],opacity=1){
		this.lines.position= this.lines.position.concat(positionTo,positionFrom);
		this.lines.color = this.lines.color.concat(color,color);
		this.lines.opacity = this.lines.opacity.concat([opacity,opacity]);
		this.lines.count++;
		return this.lines.count-1
	}

    removeAllLines(){
        this.lines.position= [];
		this.lines.color = [];
		this.lines.opacity = [];
		this.lines.count=0;
    }
    

    addRectangle(position,width,height,color=[0,0,0],opacity=1){	
		this.rects.position.push(position);
		this.rects.position.push([position[0],position[1]+height]);
		this.rects.position.push([position[0]+width,position[1]+height]);
		this.rects.position.push([position[0]+width,position[1]+height]);
		this.rects.position.push([position[0]+width,position[1]]);
		this.rects.position.push([position[0],position[1]]);
		for (let a=0;a<6;a++){
			 this.rects.opacity.push(opacity);
		}
		const c  = [color[0],color[1],color[2]];
		
		for (let a=0;a<6;a++){
			 this.rects.color.push(c);
		}
		this.rects.count++;
		return this.rects.count-1;
	}	
    


	changeImage(image,config,index){
		this.images.props.text[index]=this.regl.texture({data:image,min:"linear"});
		this.resizeImage([config.width,config.height],index);
		this.setImagePosition(config.position[0],config.position[1],index);
	}

	changeImageOpacity(index,opacity){
		const st = index*6;
		const end = st+6;
		for (let a=st;a<end;a++){
			this.images.opacity[a]=opacity;
	   	}

	}

	moveImage(x,y){
		const st = this.movingImage*6;
		const en = st+6;
		for (let n= st;n<en;n++){
			const pos = this.images.position[n];
			pos[0]+=x;
			pos[1]+=y;
		}
		const xy =this.images.props.x_y[this.movingImage];
		xy[0]+=x;
		xy[1]+=y;
	}

	setImagePosition(x,y,index){
		const p = this.images.position;
		const i = index*6;
		const width= this.images.props.w_h[index][0];
		const height = this.images.props.w_h[index][1];
		y=-y-height;
		this.images.props.x_y[index]=[x,y];
		p[i]= [x,y]
		p[i+1]=[x,y+height];
		p[i+2]=[x+width,y+height];
		p[i+3]=[x+width,y+height];
		p[i+4]= [x+width,y];
		p[i+5]=[x,y];
	}

	getImageDetails(index){
		const wh = this.images.props.w_h[index];
		const xy =this.images.props.x_y[index];
		return {
			width:wh[0],
			height:wh[1],
			position:[xy[0],xy[1]],
			index:index
		}
	}

	resizeImage(factor,index){
		const st = index*6;
	    const wh = this.images.props.w_h[index];
		const xy =this.images.props.x_y[index];
		let width,height = null;
		if (Array.isArray(factor)){
			width = factor[0];
			height= factor[1];
			wh[0]=width;
			wh[1]=height;
		}
		else{
			width = wh[0]*=factor;
			height= wh[1]*=factor;
		}
	
		const x= xy[0];
		const y =xy[1];
		this.images.position[1]=[x,y+height];
		this.images.position[2]=[x+width,y+height];
		this.images.position[3]=[x+width,y+height];
		this.images.position[4]=[x+width,y];		
	}

	addImage(image,config){
		const c = config;
		let self =this;
	
		let x = c.position[0];
		let y = -c.position[1];
		const height=c.height;
		const width =c.width;
		y=y-height;
		var image_index=this.images.position.length;
		this.images.position.push([x,y]);
		this.images.position.push([x,y+height]);
		this.images.position.push([x+width,y+height]),
		this.images.position.push([x+width,y+height]);
		this.images.position.push([x+width,y]);
		this.images.position.push([x,y]);
		this.image_position=[[-1,1],[-1,0],[0,0],[0,0],[0,1],[-1.0]];
		
		const opacity=c.opacity==null?1:c.opacity;
		for (var a=0;a<6;a++){
			 this.images.color.push([1,1,1]);
			 this.images.opacity.push(opacity);
		}
		this.images.count++;
		this.images.props.x_y.push([x,y]);
		this.images.props.w_h.push([width,height]);
		this.images.props.text.push(self.regl.texture({data:image,min:"linear"}));
		return this.images.count-1;
	}

	removeImages(){
		this.images.position=[];
		this.images.color = [];
		this.images.opacity = [];
		this.images.props.x_y=[];
		this.images.props.w_h=[];
		this.images.props.text=[];
		this.images.count=0;

	}

    addSquares(config){
        this.squares.x_pos=config.x;
        this.squares.y_pos=config.y;	
        const len = config.x.length;
        this.squares.localFilter= config.localFilter;
		this.squares.globalFilter= config.globalFilter;
        this.squares.pick_color= new Uint8Array(len*3);
        this.squares.color=config.colors;
        for (let n=0;n<len;n++){
            let p = n*3;
            let pb= this._getRGBFromIndex(n+1);
            this.squares.pick_color[p]=pb[0];
            this.squares.pick_color[p+1]=pb[1];
            this.squares.pick_color[p+2]=pb[2]
        }
        this.squares.count=len;
    }

	addCircles(config){
        this.circles.x_pos=config.x;
        this.circles.y_pos=config.y;
		if (this.mode==="3d"){
			this.circles.z_pos=config.z;
		}
        const len = config.x.length;
        this.circles.localFilter= config.localFilter;
		this.circles.globalFilter= config.globalFilter;
        this.circles.pick_color= new Uint8Array(len*3);
        this.circles.color=new Uint8Array(len*3);
        for (let n=0;n<len;n++){
            let p = n*3;
			const col = config.colorFunc(n);
            this.circles.color[p]=col[0];
            this.circles.color[p+1]=col[1];
            this.circles.color[p+2]=col[2];
            let pb= this._getRGBFromIndex(n+1);
            this.circles.pick_color[p]=pb[0];
            this.circles.pick_color[p+1]=pb[1];
            this.circles.pick_color[p+2]=pb[2]
        }
        this.circles.count=len;
    }

	updateSize(newSize,config){

		this.circles.x_pos=config.x;
        this.circles.y_pos=config.y;
		if (this.mode==="3d"){
			this.circles.z_pos=config.z;
		}

        this.circles.localFilter= config.localFilter;
		this.circles.globalFilter= config.globalFilter;
		let newPickColor= new Uint8Array(newSize*3);
		let newColor=new Uint8Array(newSize*3);
		newColor.set(this.circles.color);
		newPickColor.set(this.circles.pick_color);
		this.circles.color= newColor;
		this.circles.pick_color=newPickColor;
		for (let n = this.circles.count;n<newSize;n++){
			let p = n*3;
			const col = config.colorFunc(n);
            this.circles.color[p]=col[0];
            this.circles.color[p+1]=col[1];
            this.circles.color[p+2]=col[2];
            let pb= this._getRGBFromIndex(n+1);
            this.circles.pick_color[p]=pb[0];
            this.circles.pick_color[p+1]=pb[1];
            this.circles.pick_color[p+2]=pb[2]
		}
		this.circles.count=newSize;
	}


	_getRGBFromIndex(index){
	
		var b = Math.floor(index/65536);
		var temp = index%65536;
		var g= Math.floor(temp/256);
		var r = temp%256;
		return [r,g,b];
            
	}

	_getIndexFromRGB(rgb){
    	return (rgb[2]*65536)+(rgb[1]*256)+rgb[0];    
	}

	_drawPickBuffer(){
       	this.regl.clear({
        	color: [0, 0, 0, 0],
			depth: 1,
        	framebuffer:this.pickbuffer
    	});
     	this._drawObjects(true);
    
	}
	//refesh all 
	//in_view only those in view
	refresh(){
    	//this.label_context.clearRect(0, 0, this.width, this.height);
		
		this.regl.clear({
			color: [0, 0, 0, 0],
			depth:1
			});
    	this._drawObjects(false);
    	this._drawPickBuffer();
    	this.label_context.font = "30px Arial";
   
	}

	setCamera(distance,theta,phi){
		if (this.mode!=="3d"){
			return;
		}
		const s = this.camera.cameraState;
		s.distance=Math.log(distance);
		s.theta=theta;
		s.phi=phi;
		this.camera.updateCamera();
	}

	getCameraSettings(){
		if (this.mode!=="3d"){
			return;
		}
		const s = this.camera.cameraState;
		return{
			distance:Math.pow(Math.E,s.distance),
			theta:s.theta,
			phi:s.phi

		}
	
	}

	zoom(amount){
    	this.x_scale*=amount;
    	this.y_scale*=amount;
   		this._drawObjects(false);
	}



	setFilterAction(type){
		if (type==="grey"){
			this.hideOnFilter=false;

		}
		else{
			this.hideOnFilter=true;
		}
	}


	_drawObjects(buffer){
	
		let obj=null;
	
		if (this.mode==="3d"){
			const cProj = this.camera.getProjection();
			obj={
				cameraProjection:cProj.projection,
				cameraView:cProj.view,
				cameraDistance:cProj.distance,
				axisScales:this.axisScales,
				cameraEye:cProj.eye
 
			 };
		}
		else{
			obj = {
				x_scale:this.x_scale,
				y_scale:this.y_scale,
				offset:this.offset
			}
		}
		obj.point_scale =this.pointScale?(this.x_scale+this.y_scale)/2/this.pointScale:1;
		obj.point_radius=this.pointRadius;
		obj.is_filtered=this.isFiltered?1:0;
		obj.point_opacity=this.pointOpacity;
		obj.hide_on_filter=this.hideOnFilter?1:0;
	
				
		for (let i of this.draw_order){
			var type =this.object_types[i];
			//no objects of this type
	
			if (!type.data.count){
				continue;
			}
	

			if (buffer){
				if (!type.properties.pick_color){
					continue;
				}
				buffer=this.pickbuffer;
			}
			else{
				buffer=null;
			}

			obj.buffer=buffer;
			obj.count=type.data.count * type.vertices;
			obj.primitive=type.primitive;
			obj.is_buffer=buffer?1:0;
			  
		
			//images special case - a draw commnad for each image
			if (i===2){
				for (let i=0;i<type.data.count;i++){
					for (let prop in type.properties){
						obj[prop]= type.data[prop].slice(i*6,(i*6)+6)
					}
					for (let prop in this.images.props){
						obj[prop]=this.images.props[prop][i];
					}
					type.method(obj)

				}
				continue;
			}
			if(i===3){
				if (buffer){
					//continue;
				}
				
			
				if (this.x_scale>=this.y_scale){
					obj.side_length=20*(this.x_scale)
					obj.bottom_clip=(this.y_scale/this.x_scale);
					obj.right_clip=1.0;
				}
				else{

					obj.side_length=20 *(this.y_scale);
					obj.right_clip=(this.x_scale/this.y_scale);
					obj.bottom_clip=1.0;


				}
			
			}
				
			
			for (var prop in type.properties){

				if (buffer){
					//swap color for pock buffer
					if (prop==='pick_color'){
						obj['color']=type.data[prop];
						continue;
					}
					if (prop==='color'){
						continue;
					}
					obj[prop]=type.data[prop];
				}
				else{
					if (prop==='pick_color'){
						continue;
					}
					obj[prop]=type.data[prop];
				}
			}
			type.method(obj);
			
		}
		if (this.highlightPoints && !(buffer)){
			obj.x_pos=this.highlightPoints.x_pos;
			obj.y_pos=this.highlightPoints.y_pos;
			obj.color=this.highlightPoints.color;
			obj.primitive="points";
			obj.count=this.highlightPoints.length;
			if (this.mode==="3d"){
				obj.z_pos=this.highlightPoints.z_pos;
				this.__draw3DHighlights(obj);
			}
			else{
				this.__drawHighlights(obj);
			}
			
		}
	}


	getObjectsInRange(x,y,w,h){
		console.log(`${x},${y},${w},${h}`);
		const  max = w*h*4;
    	var pixels = this.regl.read({
			x: x,
			y: this.height -y-h,
			width:w,
			height:h ,
		
			framebuffer: this.pickbuffer
    	});
		const s = new Set()
		for (let i=0;i<max;i+=4){
			const index =this._getIndexFromRGB([pixels[i],pixels[i+1],pixels[i+2]]);
			if (index!==null){
				s.add(index-1)
			}
		}
		return s;

	}


	_getObjectAtPosition(position){
		if(position[1]<=0 || position[0]<=0 || position[0]>=this.width || position[1]>=this.height){
			return;
		}
		try{
			var pixel = this.regl.read({
				x: position[0],
				y: this.height - position[1],
				width: 1,
				height: 1,
				framebuffer: this.pickbuffer
			});
		}catch(e){
			console.log(position);
			return;
		}
		
		if (pixel[0]===0 && pixel[1]===0 && pixel[2]===0){
			return null;
		}
		//console.log(pixel);
		var index = this._getIndexFromRGB(pixel);
		if (index>0){
			return index-1;
		}
		return null;
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
			ondragstart:()=>this.brush_moving=true,
			ondragend:()=>this._brushingStopped(),
			doc:this.__doc__

		});
	
		this.brush={origin:origin,div:div,resizing:true};
	}

	_setUpPolyBrush(pos){
		
		this.poly_brush={
			points:[pos],
			active:true,
		}
		let ctx= this.label_context;
		ctx.strokeStyle = getComputedStyle(this.canvas).getPropertyValue('color');
		ctx.beginPath()
		ctx.moveTo(pos[0],pos[1]);
	


	}

	_extendPolyBrush(pos, threshold = 15){
		const ctx= this.label_context;
		if (this.poly_brush.points.length) {
			const prev = this.poly_brush.points[this.poly_brush.points.length-1];
			const dx = prev[0]-pos[0];
			const dy = prev[1]-pos[1];
			const d = Math.sqrt(dx**2 + dy**2);
			if (d < threshold) return;
		}
		
		ctx.lineTo(pos[0],pos[1]);
		ctx.stroke();
		this.poly_brush?.points.push(pos);
	}

	_finishPolyBrush(){
		if (this.poly_brush.points.length<3){
			this.clearBrush();
			return;
		}
		let ctx= this.label_context;	
		const start = this.poly_brush.points[0];
		ctx.lineTo(start[0], start[1]);
		ctx.stroke();
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
			setTimeout(()=>this.handlers.brush_stopped[i](poly,true), 0);
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
		if (!this.brush) return; //'cannot set properties of null' was possible here
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
		if (this.mode==="3d"){
			const s = this.getObjectsInRange(x,y,w,h);
			
			for (var i in this.handlers.brush_stopped){
				this.handlers.brush_stopped[i](s);
			}


		}
		else{
			let lt =this._getActualPosition([x,y]);
			let br = this._getActualPosition([x+w,y+h]);
			let info = {x_min:lt[0],x_max:br[0],y_max:-lt[1],y_min:-br[1]};
			for (var i in this.handlers.brush_stopped){
				this.handlers.brush_stopped[i](info);
			}

		}
		
	

	}

	setGlobalFilter(filter){
		this.circles.globalFilter=filter;
	}


    _finish(evt){
        if (this.config.brush){
            this.div_container.style.cursor="crosshair";
        }
        if (this.brush && this.brush.resizing){
            this._brushingStopped();
            //return;
        }
        if (this.poly_brush && this.poly_brush.active){
            this._finishPolyBrush();
        }
        
        
            //an object has finshed its drag
            if (this.object_clicked){
                this.object_clicked=null;
                this.refresh(true);              
            }
            
            //update which objects are now in view
            else{
                let ret = this.getRange();
                if (this.loop){
                    this.loop.cancel();
                    this.loop=null;
                }
            
                
            
                this.refresh();
                ret.imageMoved= this.imageMoved;
                this.imageMoved=null;
                for (var i in this.handlers.panning_stopped){
                        this.handlers.panning_stopped[i](ret);
                }
            }
            this.dragging=false;
    
        this.object_clicked=null;
        this.mouse_position=null;   

    }


	_addHandlers(){
		var self=this;
		//some listeners are on window while brush is active, and removed when brush is finished
		//allows dragging outside the canvas
		const mousemove = (e) => {
			if (self.brush){
				if (self.brush.resizing){
					let origin =self.brush.origin;
					let now = self._getMousePosition(e);
					let left = Math.round((origin[0]<now[0]?origin[0]:now[0]))+"px";
					let top =Math.round((origin[1]<now[1]?origin[1]:now[1]))+"px";
					let width = (Math.abs(origin[0]-now[0]))+"px";
					let height= (Math.abs(origin[1]-now[1]))+"px";
					self.brush.div.style.top=top;
					self.brush.div.style.left=left;
					self.brush.div.style.width=width;
					self.brush.div.style.height=height;
					
					return;
				}
				else if (self.brush.moving){
					self.dragging=false;
					return;

				}
				
			}
			if (self.poly_brush && self.poly_brush.active){
				let pt =self._getMousePosition(e);
				self._extendPolyBrush(pt);
			}
			//is this a drag or just a click without the mouse moving
			if (self.mouse_position &&  ! self.dragging){
				var x_amount= (e.pageX-self.mouse_position[0]);
				var y_amount = (e.pageY-self.mouse_position[1]);
				if (Math.abs(x_amount) > 3 || Math.abs(y_amount)>3){
					self.dragging = true;
				}
			}

			if (self.dragging){
				const dx= e.pageX-self.mouse_position[0];
				const dy =e.pageY-self.mouse_position[1]

				if (self.mode==="3d"){
					self.camera.mouseChange(dx/self.width,dy/self.height);
				}
				else{		
					const x_amount= dx/self.x_scale;
					const y_amount = dy/self.y_scale;
					if (self.movingImage != null && e.ctrlKey && self.images.count>0){
						self.moveImage(x_amount,y_amount);
						self.imageMoved=self.getImageDetails(self.movingImage);
					}
					else if(self.offsets && e.which ==3 		){
						self.alterOffsets(x_amount,y_amount);
					
					}
					else{
						if (!self.config.lock_x_axis){
							self.offset[0]+=x_amount;
						}		
						self.offset[1]+=y_amount;
						for (var i in self.handlers.pan_or_zoom){
							self.handlers.pan_or_zoom[i](self.offset,self.x_scale,self.y_scale);
						}

					}
				}
			
				if (!self.loop){
					//self.label_context.clearRect(0, 0, self.width, self.height);
					self.loop = self.regl.frame(()=>{
						self.regl.clear({
							color: [0, 0, 0, 0],
							depth:1
					  	});
						self._drawObjects(false);
					});

				}
				self.mouse_position[1]=e.pageY;
				self.mouse_position[0]=e.pageX;      
			}
			//no drag event going on call any listners if mouse over/out an object
			else{
				var position =self._getMousePosition(e);
				var obj = self._getObjectAtPosition(position);
				if (obj!=null && self.object_mouse_over==null){
					for (var i in self.handlers['object_over']){
						self.handlers.object_over[i](e,obj);		           
					}
				
					self.object_mouse_over=obj;
					
				}
				else if (obj==null && self.object_mouse_over){
					for (var i in self.handlers['object_out']){
						self.handlers.object_out[i](e,self.object_mouse_over);
					}
				
					
					self.object_mouse_over=null;

				}
				//move directly from one object to another
				else if(obj!=null && (obj!==self.object_mouse_over)){
					for (var i in self.handlers['object_over']){    
						self.handlers.object_over[i](e,obj);  
					}
					self.object_mouse_over=obj;
				}         
			}
		};

		const mouseup = (evt) => {
			window.removeEventListener("mousemove",mousemove);
			window.removeEventListener("mouseup",mouseup);

            this._finish(evt);
            if (!this.dragging && !(this.brush || this.poly_brush) && evt.button===0){
                    var position =this._getMousePosition(evt);
                    var obj = this._getObjectAtPosition(position);
                    if (obj){
                        for (var i in this.handlers.object_clicked){
                            this.handlers.object_clicked[i](obj);
                        }  
                    }
            }
        
            
			/*//just a click event - inform handlers
			if (self.config.brush){
				self.div_container.style.cursor="crosshair";
			}
			if (self.brush && self.brush.resizing){
				self._brushingStopped();
				//return;
			}
			if (self.poly_brush && self.poly_brush.active){
				self._finishPolyBrush();
			}
			if (!self.dragging && !(self.brush || self.poly_brush)){
				//if (self.object_clicked){
					var position =self._getMousePosition(evt);
					var obj = self._getObjectAtPosition(position);
					if (obj){
						for (var i in self.handlers.object_clicked){
							self.handlers.object_clicked[i](obj);
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
					if (self.loop){
						self.loop.cancel();
						self.loop=null;
					}
				
					
					
					  

					//self._drawPickBuffer(false);
					//self._getObjectsInView();
				
					self.refresh();
					ret.imageMoved= self.imageMoved;
					self.imageMoved=null;
					for (var i in self.handlers.panning_stopped){
							self.handlers.panning_stopped[i](ret);
					}
				}
				self.dragging=false;
			}
			self.object_clicked=null;
			self.mouse_position=null;   */
		};

		this.div_container.addEventListener('wheel', function(event){
			event.preventDefault();
			var position =self._getActualPosition(self._getMousePosition(event));
			let wdelta=0;
			if (event.wheelDelta > 0 || event.detail < 0) {
					self.zoom_amount+=0.01;
					wdelta=-30;
			}
			else {
					self.zoom_amount-=0.01;
					wdelta=30;
			}
			if (self.mode==="3d"){
				self.camera.mouseWheel(0,wdelta)
			}
			else {
				
				if (self.movingImage != null && event.ctrlKey){
					self.resizeImage(1+self.zoom_amount,self.movingImage);
					const new_position=self._getActualPosition(self._getMousePosition(event));
					self.moveImage(new_position[0]-position[0],new_position[1]-position[1]);
					self.imageMoved =self.getImageDetails(self.movingImage);
					
				}
				else{

					if (!self.config.lock_x_axis){
						self.x_scale*=(1+self.zoom_amount);
					}
					if (!self.config.lock_y_axis){
						self.y_scale*=(1+self.zoom_amount);
					}			
					const new_position=self._getActualPosition(self._getMousePosition(event));
					self.offset[0]+=new_position[0]-position[0];
					self.offset[1]+=new_position[1]-position[1];
					for (var i in self.handlers.pan_or_zoom){
						self.handlers.pan_or_zoom[i](self.offset,self.x_scale,self.y_scale);					
					}
				}
			}
           
			if (!self.loop){
		        self.label_context.clearRect(0, 0, self.width, self.height);
				self.loop = self.regl.frame(function(){
					self.regl.clear({
						color: [0 ,0, 0, 0],
						depth:1,
					  })
					self._drawObjects(false);
				});
			}


			
			//clear the timeout user has not finished zooming
			clearTimeout(self._timer);
			//when user finishes call the esxpensive methods;
			self._timer= setTimeout(function() {
				self.zoom_amount=0;
				self.loop.cancel();
				self.loop=null;
				self._drawPickBuffer(false);
				if (self.brush){

				}
				let ret = self.getRange();
				ret.imageMoved = self.imageMoved;
				self.imageMoved =null;
				for (let name in self.handlers.zoom_stopped){	
					self.handlers.zoom_stopped[name](ret);
				}
                		
				self.refresh();
			}, 350);

		});
		this.div_container.addEventListener("mousedown",function (evt){
			window.addEventListener("mousemove",mousemove);
			window.addEventListener("mouseup",mouseup);

			if (evt.which===3){
				//add right click behaviour
			}
			let otherKey = evt.shiftKey || evt.ctrlKey || evt.which==2 || evt.which==3;
			//create brush
			if ((self.config.brush && !(otherKey)) || (!(self.config.brush)&&otherKey)){
				let origin = self._getMousePosition(evt);
				if (self.config.brush=="default"){
					let t = evt.target;
					if(t.classList.contains("wgl2di-brush")){
						return;
					}
					if (self.brush){
						self.clearBrush();
					}
					let origin =self._getMousePosition(evt);
					self._setUpBrush(origin);
					return;
				}
				else if (self.config.brush==="poly"){
				    if(self.poly_brush){
			    	    self.clearBrush();
			        }
			        self._setUpPolyBrush(origin)		    
			        return;
				}
			}
		
			if (evt.otherKey){
				self.div_container.style.cursor="move";
			}
			
			var position =self._getMousePosition(evt);
			self.mouse_position= [evt.pageX, evt.pageY];
			self.clearBrush();
			evt.preventDefault();
		});
	}


	_initDrawMethods(){
		var self = this;

		//loading images
		this.__drawCircles = this.regl({
			depth:self.draw_options.depth,
	        blend:self.draw_options.blend,
	
			frag: 
		    `precision mediump float;
			varying vec3 fragColor;
			varying float op;
			varying float has_border;
			varying float grey_out;
			uniform float is_buffer;
			
			void main(){
			
				vec2 cxy = 2.0 * gl_PointCoord - 1.0;
                
				float r = dot(cxy, cxy);
				if (r > 1.0) {
					discard;
				}
				else{					
					if(r>0.60 && has_border==-1.0 && is_buffer==0.0 && grey_out==0.0){
					    gl_FragColor=vec4(0.1,0.1,0.1,1.0);
					}
					else{
						if (grey_out==1.0 && is_buffer==0.0){
							gl_FragColor = vec4(0.2,0.2,0.2,0.4);
						}
						else{
							gl_FragColor = vec4(fragColor,is_buffer==1.0?1.0:op);
						}
						
					}
				}
				
			}`,

			vert: 
	        `attribute float x_pos;
            attribute float y_pos;
			attribute vec3 color;
			attribute float lFilter;
			attribute float gFilter;
			varying vec3 fragColor;
			varying float op;
			varying float has_border;
			varying float grey_out;
			uniform float x_scale;
			uniform float y_scale;
			uniform vec2 offset;
			uniform float stage_height;
			uniform float stage_width;
			uniform float point_radius;
			uniform float is_filtered;
			uniform float point_opacity;
			uniform float hide_on_filter;
			uniform float point_scale;
			
			vec2 normalizeCoords(float posX, float posY){
				float x = (posX+offset[0])*x_scale;
				float y = (posY+offset[1])*y_scale;
				return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
			}
	
			void main() {
				float dd = 0.0;
				grey_out=0.0;
				if (lFilter==2.0){
					return;
				}
			    if (gFilter>0.0 && gFilter != lFilter) {
					if(hide_on_filter==1.0){
						return;
					}
					else{
						grey_out=1.0;
					}
				}
				
			
				float r=point_radius;
			
				gl_PointSize = r*point_scale;
				vec3 c = vec3(255.0,255.0,255.0);
				if (color==c){
					op=0.1;
					fragColor=vec3(0.1,0.1,0.1);
				}
				else{
					op=point_opacity;
					fragColor = color/255.0;
				}
			
				
				
				has_border=lFilter-is_filtered;
			
				
				vec2 real_position = normalizeCoords(x_pos,-y_pos);
				gl_Position = vec4(real_position, 0.0, 1.0);
			}
			`,

			attributes: {
				x_pos: self.regl.prop('x_pos'),
                y_pos:self.regl.prop("y_pos"),
				color: self.regl.prop('color'),
				lFilter:self.regl.prop("localFilter"),
				gFilter:self.regl.prop("globalFilter"),

			},
		
			uniforms: {
				x_scale:self.regl.prop('x_scale'),
				y_scale:self.regl.prop('y_scale'),
				stage_width: self.regl.context("viewportWidth"),
				stage_height: self.regl.context('viewportHeight'),
				offset:self.regl.prop("offset"),
				point_radius:self.regl.prop("point_radius"),
				point_opacity:self.regl.prop("point_opacity"),
				is_filtered:self.regl.prop("is_filtered"),
				is_buffer:self.regl.prop("is_buffer"),
				hide_on_filter:self.regl.prop("hide_on_filter"),
				point_scale:self.regl.prop('point_scale')
			},

			count:  self.regl.prop('count'),
			primitive: self.regl.prop('primitive'),
			framebuffer:self.regl.prop("buffer")
		});

		this.__drawHighlights = this.regl({
			depth:self.draw_options.depth,
	        blend:self.draw_options.blend,
			depth:self.draw_options.depth,
	        blend:self.draw_options.blend,
	
			frag: 
		    `precision highp float;
			varying vec3 fragColor;
			
			void main(){

				vec2 cxy = 2.0 * gl_PointCoord - 1.0;
                
				float r = dot(cxy, cxy);
				if (r > 1.0) {
					discard;
				}
				else{					
					if(r>0.50) {
					    gl_FragColor=vec4(0.1,0.1,0.1,1.0);
					}
					else{
						gl_FragColor = vec4(fragColor,1);
					}
				}	
			}`,

			vert: 
	        `attribute float x_pos;
            attribute float y_pos;
			attribute vec3 color;
			uniform float point_radius;
			
			varying vec3 fragColor;
		
			uniform float x_scale;
			uniform float y_scale;
			uniform vec2 offset;
			uniform float stage_height;
			uniform float stage_width;
			uniform float point_scale;
			
			vec2 normalizeCoords(float posX, float posY){
				float x = (posX+offset[0])*x_scale;
				float y = (posY+offset[1])*y_scale;
				return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
			}
	
			void main() {
			
				float pr = point_radius*point_scale;
				pr = pr<10.0?10.0:pr;
				gl_PointSize = pr;
				fragColor = color/255.0;
					
				vec2 real_position = normalizeCoords(x_pos,-y_pos);
				gl_Position = vec4(real_position, 0.0, 1.0);
			}
			`,

			attributes: {
				x_pos: self.regl.prop('x_pos'),
                y_pos:self.regl.prop("y_pos"),
				color: self.regl.prop('color'),

			},
		
			uniforms: {
				x_scale:self.regl.prop('x_scale'),
				y_scale:self.regl.prop('y_scale'),
				stage_width: self.regl.context("viewportWidth"),
				stage_height: self.regl.context('viewportHeight'),
				offset:self.regl.prop("offset"),
				point_radius:self.regl.prop("point_radius"),
				point_scale:self.regl.prop('point_scale')
			
			},

			count:  self.regl.prop('count'),
			primitive: self.regl.prop('primitive')
			
		});


		this.object_types[0]['method']=this.__drawCircles;


	

		this.__drawLines = this.regl({

				// fragment shader
				frag: `precision highp float;
						varying vec3 fragColor;
						void main () {
							 gl_FragColor = vec4(fragColor,1);
						}`,


				vert: `
						attribute vec2 position;
						attribute vec3 color;
						attribute float opacity;
						uniform float x_scale;
						uniform float y_scale;
						uniform vec2 offset;
						uniform float stage_height;
						uniform float stage_width;
						varying vec3 fragColor;
						//varying float op;
						vec2 normalizeCoords(vec2 position){	    
							float x = (position[0]+offset[0])*x_scale;
							float y = (position[1]+offset[1])*y_scale;
				            return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
						}
						void main () {
							if (opacity==0.0){
								return;
							}
							fragColor=color/255.0;
							vec2 norm_pos =normalizeCoords(position);
							gl_Position = vec4(norm_pos, 0.0, 1.0);
						}`
				,
				attributes: {
					position: self.regl.prop("position"),
					color:self.regl.prop("color"),
					opacity:self.regl.prop("opacity")
				},

				uniforms: {
					  x_scale:self.regl.prop('x_scale'),
					  y_scale:self.regl.prop('y_scale'),
					  stage_height:self.regl.context("viewportHeight"),
					  stage_width:self.regl.context("viewportWidth"),
					  offset:self.regl.prop("offset")
				},
				primitive:self.regl.prop("primitive"),
				count:self.regl.prop("count")



			});
		this.object_types[1]['method']=this.__drawLines;
        this.object_types[4]['method']=this.__drawLines;
		

		this.__drawImages = this.regl({
			frag: `
				precision mediump float;
				uniform sampler2D text;
				varying vec2 uv;
				varying vec3 fragColor;
				varying float op;
				void main () {
					gl_FragColor = vec4(fragColor,op)*texture2D(text, uv);							
				}`,

 			vert: `
				precision mediump float;

				attribute vec2 position;
				attribute vec3 color;
				attribute float opacity;

				uniform vec2 x_y;
				uniform vec2 w_h;
				uniform float stage_height;
				uniform float stage_width;
				uniform float x_scale;
				uniform float y_scale;
				uniform vec2 offset;

				varying vec2 uv;
				varying vec3 fragColor;
				varying float op;
		
				
				vec2 normalizeCoords(vec2 pos){
					float x = (pos[0]+offset[0])*x_scale;
					float y = (pos[1]+offset[1])*y_scale;
					return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
				}

				void main () {
					if (opacity==0.0){
						return;
					}
					op=opacity;
					
					vec2 new_pos=normalizeCoords(position);
					fragColor = color;
				
					float x_factor = 1.0/(((w_h[0]*x_scale)/stage_width)*2.0);
					float y_factor = 1.0/(((w_h[1]*y_scale)/stage_height)*2.0);
					

					float x_offset=(((x_y[0]+offset[0])*x_scale)/stage_width)*2.0*x_factor;
					float y_offset=(((x_y[1]+offset[1])*y_scale)/stage_height)*2.0*y_factor;
					uv = vec2((new_pos[0]*x_factor)+x_factor-x_offset,-(new_pos[1]*y_factor)+y_factor-y_offset);
					gl_Position = vec4(new_pos, 0.0, 1.0);

				}`,

			attributes: {
    			position:self.regl.prop("position"),
    			color:self.regl.prop("color"),
    			opacity:self.regl.prop("opacity")
   			},

  			uniforms: {
    			stage_height:self.regl.context("viewportHeight"),
    			stage_width:self.regl.context("viewportWidth"),
    			w_h:self.regl.prop("w_h"),
    			x_y:self.regl.prop("x_y"),
    			text:self.regl.prop("text"),
				
				x_scale:self.regl.prop('x_scale'),
				y_scale:self.regl.prop('y_scale'),
				offset:self.regl.prop("offset")
    			

  			},

  			count: self.regl.prop("count"),
			depth:self.draw_options.depth,
			blend:self.draw_options.blend
		});

		this.object_types[2]['method']=this.__drawImages;


		this.__drawSquares = this.regl({
			frag: `
				precision highp float;
				varying vec3 fragColor;
				varying float r_clip;
				varying float b_clip;
				uniform int is_buffer;
				void main(){
					if (gl_PointCoord[1]<1.0-b_clip){
						discard;
					}
					if (gl_PointCoord[0]> r_clip){
						discard;
					}				
					gl_FragColor = vec4(fragColor,1);					
				}
				`,

			vert:` 
				attribute float x_pos;
				attribute float y_pos;			
				attribute vec3 color;
			
				varying float r_clip;
				varying float b_clip;
				varying vec3 fragColor;

				uniform float x_scale;
				uniform float y_scale;
				uniform vec2 offset;
				uniform float stage_height;
				uniform float stage_width;
				uniform float right_clip;
				uniform float bottom_clip;
				uniform float side_length;
				
				vec2 normalizeCoords(float posX,float posY){
					float x = (posX+offset[0])*x_scale;
					float y = (posY+offset[1])*y_scale;
					return vec2(2.0 * ((x / stage_width) - 0.5),-(2.0 * ((y / stage_height) - 0.5)));
				}
				void main() {
					gl_PointSize = side_length;
					fragColor = color/255.0;
					r_clip=right_clip;
					b_clip= bottom_clip;
					float y = y_pos;
					float x = x_pos;
					if (bottom_clip!=1.0){
						y-=10.0*(1.0/b_clip)-10.0;
					}
					if (r_clip!=1.0){
						x+=10.0*(1.0/r_clip)-10.0;
					}

					vec2 real_position = normalizeCoords(x,y);
				
					gl_Position = vec4(real_position, 0.0, 1.0);
				}
			`
			,
				
			attributes: {
				x_pos: self.regl.prop('x_pos'),
                y_pos:self.regl.prop("y_pos"),
				color: self.regl.prop('color')	
			},
	
			uniforms: {
				x_scale:self.regl.prop('x_scale'),
				y_scale:self.regl.prop('y_scale'),
				stage_height:self.regl.context("viewportHeight"),
				stage_width:self.regl.context("viewportWidth"),
				offset:self.regl.prop("offset"),
				is_buffer:self.regl.prop("is_buffer"),
				right_clip:self.regl.prop("right_clip"),
				bottom_clip:self.regl.prop("bottom_clip"),
				side_length:self.regl.prop('side_length'),
			},
	
			count:  self.regl.prop('count'),
			primitive: self.regl.prop('primitive'),
			framebuffer:self.regl.prop("buffer")
		});
		this.object_types[3]['method']=this.__drawSquares;





		this.__draw3DCircles = this.regl({
            frag: `
            precision highp float;
            varying vec3 fragColor;
			varying float op;
			varying float has_border;
			uniform float is_buffer;
			varying float grey_out;
            void main () {
			  float r = length(gl_PointCoord.xy - 0.4) ;
              if (r> 0.4) {
                discard;
              }
			  else{
				if(r>0.3 && has_border==-1.0 && is_buffer==0.0 && grey_out==0.0){
					gl_FragColor=vec4(0.1,0.1,0.1,1.0);
				}
				else{
					if (grey_out==1.0 && is_buffer==0.0){
						gl_FragColor = vec4(0.2,0.2,0.2,0.4);
					}
					else{
						gl_FragColor = vec4(fragColor,is_buffer==1.0?1.0:op);
					}
					
				}

			  }
             
            }
            `,
          
            vert: `
            precision highp float;

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
			uniform float point_radius;
			uniform float point_opacity;
			uniform float is_buffer;
			uniform float is_filtered;
			uniform float hide_on_filter;
			uniform vec3 axis_scales;


            varying vec3 fragColor;
			varying float op;
			varying float has_border;
			varying float grey_out;

            void main () {
				grey_out=0.0; 
			    if (gFilter>0.0 && gFilter != lFilter) {
					if(hide_on_filter==1.0){
						return;
					}
					else{
						grey_out=1.0;
					}
				}
                vec3 position = vec3(-x_pos *axis_scales.x,y_pos*axis_scales.y,-z_pos*axis_scales.z);
                fragColor=color/255.0;
				op = is_buffer==1.0?1.0:point_opacity;
				has_border=lFilter-is_filtered;
               
				gl_PointSize = point_radius; //(distance(eye, position.xyz) /30.0);
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
				point_radius:self.regl.prop("point_radius"),
				point_opacity:self.regl.prop("point_opacity"),
				is_buffer:self.regl.prop("is_buffer"),
                view: self.regl.prop("cameraView"),
                cdistance:self.regl.prop("cameraDistance"),
				eye:self.regl.prop("cameraEye"),
				is_filtered:self.regl.prop("is_filtered"),
				hide_on_filter:self.regl.prop("hide_on_filter"),
				axis_scales:self.regl.prop("axisScales"),

                projection: (context,props) =>{
                 return glMatrix.mat4.perspective(
                    props.cameraProjection,
                    Math.PI / 4,
                    context.viewportWidth / context.viewportHeight,
                   1,
                    100000
                  )
                 },
                
              },
          
              count:  self.regl.prop('count'),
			  framebuffer:self.regl.prop("buffer"),
			 
			  blend:self.draw_options.blend,
              primitive: 'points',
		
		});

		this.__draw3DHighlights = this.regl({
            frag: `
            precision mediump float;
            varying vec3 fragColor;
		
            void main () {
			  float r = length(gl_PointCoord.xy - 0.4) ;
              if (r> 0.4) {
                discard;
              }
			  else{
				if(r>0.3){
					gl_FragColor=vec4(0.1,0.1,0.1,1.0);
				}
				else{
					
						gl_FragColor = vec4(fragColor,1.0);
				}
					
			  }
             
            }
            `,
          
            vert: `
            precision mediump float;

            attribute float x_pos;
            attribute float y_pos;
            attribute float z_pos;
            attribute vec3 color;
			
			uniform float time;
    
            uniform mat4 view, projection;
			uniform vec3 eye;
			uniform float point_radius;
			uniform vec3 axis_scales;
		


            varying vec3 fragColor;
		

            void main () {
			
                vec3 position = vec3(-x_pos*axis_scales.x,y_pos*axis_scales.y,-z_pos*axis_scales.z);
                fragColor=color/255.0;
               
				gl_PointSize = point_radius*2.0; //(distance(eye, position.xyz) /30.0);
                gl_Position = projection * view * vec4(2.0 * position, 1.0);
            }
            `,
          
              attributes: {
              
                x_pos: self.regl.prop('x_pos'),
                y_pos:self.regl.prop("y_pos"),
                z_pos:self.regl.prop("z_pos"),
				color: self.regl.prop('color'),
				
              },
          
              uniforms: {
				point_radius:self.regl.prop("point_radius"),
		
                view: self.regl.prop("cameraView"),
                cdistance:self.regl.prop("cameraDistance"),
				eye:self.regl.prop("cameraEye"),
				axis_scales:self.regl.prop("axisScales"),
				
                projection: (context,props) =>{
                 return glMatrix.mat4.perspective(
                    props.cameraProjection,
                    Math.PI / 4,
                    context.viewportWidth / context.viewportHeight,
                    1,
                    100000
                  )
                 },
                
              },
          
              count:  self.regl.prop('count'),
			  blend:self.draw_options.blend,
              primitive: 'points',
		
		});

		this.__draw3DLines = this.regl({

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
					uniform vec3 axis_scales;
    
                    varying vec3 fragColor;
                  
                    
                    void main () {
                     
                        fragColor=color;
                        gl_Position = projection * view * vec4(2.0 * position*axis_scales*vec3(-1.0,1.0,-1.0), 1.0);
                       
                       
                    }`
            ,
            attributes: {
                position: self.regl.prop("position"),
                color:self.regl.prop("color")
            


            },
            uniforms: {
                view: self.regl.prop("cameraView"),
				axis_scales:self.regl.prop("axisScales"),
                projection: (context,props) =>{
                 return glMatrix.mat4.perspective(
                    props.cameraProjection,
                    Math.PI / 4,
                    context.viewportWidth / context.viewportHeight,
                    0.01,
                    100000
                  )
                 },
                
              },

          
            primitive:self.regl.prop("primitive"),
            framebuffer:self.regl.prop("buffer"),
            count:self.regl.prop("count"),
			depth:{enable:false}
        });	
	}
	
}


export {WGL2DI};



