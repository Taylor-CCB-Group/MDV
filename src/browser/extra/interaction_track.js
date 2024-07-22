import {TabixBedFeatureSource} from "../feature.js";
import {MLVBedTrack,MLVTrack} from "../tracks.js";

let decode_function= function(tokens,feature){
    //extract interactions with the first token
    feature.data= [parseInt(tokens[0])];
    feature.end1= parseInt(tokens[1]);
    feature.start2= parseInt(tokens[2]);
    feature.val= parseFloat(tokens[3]);
}

class InteractionTrack extends MLVBedTrack{
    _setFeatureSource(){
        this.feature_source = new TabixBedFeatureSource(this.config,decode_function);
    }



    getFeatureAt(gl, chr, coord, bpPerPixel,ctx){
        const d=  ctx.getImageData(coord.x, coord.y, 1, 1).data;
        if (d[0]!==0 || d[1]!==0 || d[2]!==255){
            return null;
        }
        
        const tol = bpPerPixel;
        const inters=[];

        const fts = this.feature_source.featureCache.queryFeatures(chr,gl-tol,gl+tol);
        for (let f of fts){
            if (gl>=f.start-tol && gl <= f.end1+tol){
                inters.push(f);
            }
            else if (gl>=f.start2-tol && gl <=f.end+tol){
                inters.push(f)
            }
        }
        return inters.length>0?inters:null;
    }

    onFeatureOver(container,{track,feature,event}){
        this.drawCurves({
            ctx:container.ctx,
            features:feature,
            bpStart:container.start,
            yOffset:container.yOffset,
            bpEnd:container.end,
            bpPerPixel:(container.end-container.start)/container.canvas.width,
            color:"red",
            region:false
        })
    }

    onFeatureOut(container,{track,feature,event}){	
        this.drawCurves({
            ctx:container.ctx,
            features:feature,
            yOffset:container.yOffset,
            bpStart:container.start,
            bpEnd:container.end,
            bpPerPixel:(container.end-container.start)/container.canvas.width,
            region:false
        });
    }

    getSettings(panel){
        const s= super.getSettings(panel);
        const v = this.config.value;
        return s.concat([
        /*{
            label:v.name,
            type:"doubleslider",
            min:v.extent[0],
            max:v.extent[1],
            doc:panel.__doc__ || document,
        
            current_value:v.cutoff, 
            func:(x,y)=>{
                v.cutoff=[x,y];
                panel.update();
            }
        },
        {
            label:"Length (kb",
            type:"doubleslider",
            min:0,
            max:1000,
            doc:panel.__doc__ || document,
        
            current_value:this.config.range, 
            func:(x,y)=>{
                this.config.range=[x,y];
                panel.update();
            }
        },*/
    {
        type:"check",
        label:"Only show local",
        current_value:this.config.onlyShowLocal,
        func:(x)=>{
            this.config.onlyShowLocal=x;
            panel.update();
        }

    }]);
    }


    drawCurves({ctx,features,bpStart,bpPerPixel,color,region,pixelWidth,yOffset=0}){

        //filtering based on config
      
        const bpEnd= pixelWidth*bpPerPixel+bpStart;
        const bottom = this.bottom+yOffset;
        const c_bottom= bottom-6;
        const local = this.config.onlyShowLocal;
        if (this.config.value){
            const mm= this.config.value.cutoff;
            const minr= this.config.range[0]*1000;
            const maxr= this.config.range[1]*1000;
            const local = this.config.onlyShowLocal;
            features = features.filter(x=>x.val>mm[0] && x.val<mm[1] && x.end-x.start>minr && x.end-x.start<maxr && (!local || (x.end1>bpStart && x.start2<bpEnd)));
        }
       
        let mod=1;
        if (features.length>2000 && !this.filter_function){
            //mod += Math.floor(features.length/2000);
        }
        ctx.beginPath();
        
        let drawn=0;
        for (let n=0;n<features.length;n+=mod){
            const feature = features[n];
            if ((this.filter_function && !this.filter_function(feature)) || (local && (feature.end>bpEnd || feature.start<bpStart))){
                feature.display=false;
                continue;
            }
            let st1 = (feature.start-bpStart)/bpPerPixel;
            const en1 = (feature.end1-bpStart)/bpPerPixel;
            let st2 = (feature.start2-bpStart)/bpPerPixel;
            const en2 = (feature.end-bpStart)/bpPerPixel;
            const p1 = st1+(en1-st1)/2;
            let p2 = st2+(en2-st2)/2;
            let h = ((feature.end-feature.start)/(bpPerPixel*200))*this.config.height;
            h=h>this.config.height*1.2?this.config.height*1.2:h;
            h = h<10?10:h;
            const c = color?color:this.color_function?this.color_function(feature):"gray";
            ctx.strokeStyle=c;
            ctx.beginPath();
            ctx.moveTo(p1,c_bottom);
            if(p2-p1 <3){
                p2=p1+3;
            }
            
            const mp = p1+(p2-p1)/2;
            ctx.bezierCurveTo(mp,c_bottom-h,mp,c_bottom-h,p2,c_bottom);
            ctx.stroke();
            drawn++;	
            
            if (region){
                ctx.fillStyle="blue";
                let w = en1-st1;
                if (w<5){
                    w=5;
                    st1-=2;
                }
                

                ctx.fillRect(st1,this.bottom-6,w,6);
                w=en2-st2;
                if (w<5){
                    w=5;
                    st2-=2;
                }
                ctx.fillRect(st2,this.bottom-6,w,6);
            }
        
        }
        
    }

    drawFeatures(options){
        this.drawCurves({
            ctx:options.context,
            features:options.features,
            bpStart:options.bpStart,
            bpPerPixel:options.bpPerPixel,
            region:true,
            pixelWidth:options.pixelWidth
        });
        
    }
}

MLVTrack.track_types["interaction_track"]={
    "class":InteractionTrack
}
export default InteractionTrack;

