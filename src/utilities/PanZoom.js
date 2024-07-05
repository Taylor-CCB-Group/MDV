
import { makeDraggable } from "./Elements";
class ImagePanZoom{
    constructor(container,image,config={}){
        this.container=container;
     
        this.img = new Image();
       
        this.img.style.position="absolute";
       
        this.img.onload= e=>{
            this.orig_dim= [this.img.naturalWidth,this.img.naturalHeight,this.img.naturalWidth/this.img.naturalHeight];
            this.container.append(this.img);
            this.fit();
        }
        makeDraggable(this.img);
        this.container.style.overflow="hidden";
        this.container.addEventListener("wheel", (event) => {
            event.preventDefault();
            this.zoom(Math.sign(event.deltaY) > 0 ? -1 : 1,event);
               
        });
        this.setImage(image)
    }

    setImage(url){
        this.img.src= url;
    }

    fit(){
        const box = this.container.getBoundingClientRect();
        this.img.height= box.height;
        this.img.width= this.img.height*this.orig_dim[2];
        if (this.img.width > box.width){
            this.img.width= box.width;
            this.img.height= this.img.width/this.orig_dim[2];
            this.img.style.top= `${(box.height-this.img.height)/2}px`;
            this.img.style.left="0px";
        }
        else{
            this.img.style.top="0px";
            this.img.style.left= `${(box.width-this.img.width)/2}px`;
        }
    }

   

    zoom(amount,event){
        amount = amount>0?1.1:.9;
        const cbox = this.container.getBoundingClientRect();
        const ibox =  this.img.getBoundingClientRect();
        const  dx = ((event.clientX-cbox.left)  - this.img.offsetLeft) * (amount-1);
        const dy = ((event.clientY-cbox.top) - this.img.offsetTop) * (amount-1);
      
        this.img.width=ibox.width*amount;
        this.img.height=ibox.height*amount;
        this.img.style.left =`${this.img.offsetLeft-dx}px`;
        this.img.style.top= `${this.img.offsetTop -dy}px`;
    }
}



export default ImagePanZoom;