import BaseChart from "./BaseChart.js";
import {createEl} from "../utilities/Elements.js";
import { timeMilliseconds } from "d3";


class ImageBox extends BaseChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config);
        const c = this.config
        this.contentDiv.style.overflowY="auto";
        c.images_per_row= c.images_per_row || 2;
        c.images = c.images || [];
        this.refreshImages();
     
        
	}

    refreshImages(){
        this.contentDiv.innerHTML="";
        this.divs= [];
        this.images=[];
        for (let n=0;n<this.config.images.length;n++){
            const i = this.config.images[n];
            const d = createEl("div",{
                style:{
                    display:"inline-block"
                }
              
            },this.contentDiv);
            createEl("div",{
                style:{
                    height:"20px",
                    fontSize:"14px"
                   
                },
                text: i[1] || ""
            },d);
            const im = createEl("img",{
                src:i[0]

            },d);
            this.divs.push(d);
            this.images.push(im);
            if (n===this.config.images.length-1){
                im.onload= ()=>{
                    this.ratio= im.width/im.height;
                    this.setImageSizes();
                }
            }
        }
    }

    setSingleImage(url,name){
        this.divs=[];
        this.images=[];
        this.contentDiv.innerHTML="";
        this.config.images_per_row=1;
        this.config.images=[[url,name]];
        this._createImages();
        this.images[0].onload= ()=>{
            const im = this.images[0]
            this.ratio= im.width/im.height;
            this.setImageSizes();
        }  
    }

    setImageSizes(){
        const imwidth = Math.round(this.width/this.config.images_per_row-20);
        const imheight= Math.round(imwidth/this.ratio);
        for (let d of this.divs){
            d.style.height=imheight+"px";
            d.style.width=imwidth+"px";
        }
        for (let i of this.images){
            i.style.height=(imheight-10)+"px";
            i.style.width=imwidth+"px";
        }
    }

    setSize(x,y){
        super.setSize(x,y);
        this.setImageSizes();
    }
}

export default ImageBox;

BaseChart.types["image_box"]={
    "class":ImageBox,
    name:"Image Box",
    params:[],
    allow_user_add:false

}