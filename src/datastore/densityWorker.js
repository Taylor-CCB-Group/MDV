import {contourDensity} from "d3-contour";
import {scaleLinear} from "d3-scale";


// biome-ignore lint/suspicious/noGlobalAssign: relatively innocuous in simple web worker
onmessage= (e)=> {
    const config = e.data[5];
    const cat = new Uint8Array(e.data[4]);
    const y =e.data[3][1]==="int32"? new Int32Array(e.data[3][0]):new Float32Array(e.data[3][0]);
    const x= e.data[2][1]==="int32"? new Int32Array(e.data[2][0]):new Float32Array(e.data[2][0]);
    const lFilter=new  Uint8Array(e.data[0]);
    const gFilter = new  Uint8Array(e.data[1]);
    const height= config.yscale[1][1];
    const yscale = scaleLinear().domain(config.yscale[0]).range(config.yscale[1]);
    const xscale = scaleLinear().domain(config.xscale[0]).range(config.xscale[1]);
    function calcDensity(category){
        return contourDensity()
        .x((d,i) => xscale(x[i]))
        .y((d,i) => yscale(y[i]))
        .weight((d,i)=>{
            if (lFilter[i]===2){
                return 0;
            }
            if (gFilter[i]!==0){
                if  (gFilter[i] !==lFilter[i] ){
                    return 0;
                }           
            }
            if (cat[i] !== category){
                return 0;
            }
            return 1;
    
        })
        .size([config.xscale[1][1],height ])
        .bandwidth(config.bandwidth)
        (gFilter);

    }
    const data=[];
    for (const i of config.categories){
        if (i===-1){
            data.push(null)
        }
        else{
            data.push(calcDensity(i))
        }
    }

    postMessage(data);
}