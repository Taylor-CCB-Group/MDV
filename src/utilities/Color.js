import {createEl,createSVGEl, makeResizable,makeDraggable} from "./Elements.js";
import {select} from "d3-selection";
import {scaleLinear} from "d3-scale";
import {axisBottom} from "d3-axis";
import {getRandomString} from "./Utilities.js";

function getColorLegend(colors,names,config={}){
    const len = colors.length;
    const height=(len*12)+ (len+2)*2;
    const container = createEl("div",{
        styles:{
            width:"120px",
            height:((height>215?215:height)+35)+"px",
            position:"absolute",
            border:"0.5px solid black" 
        }
    });

    createEl("div",{
        styles:{
            height:"20px",
        },
        classes:["legend-title"],
        text:config.label
    },container);
    const body = createEl("div",{
        styles:{
            overflowY:"auto",
            overflowX:"hidden",
            "height":"calc(100% - 25px)",
            "width":"100%",
        },
        classes:["legend-body"]
    },container);

    const legend= createSVGEl("svg",{
        height:height,
        width:180,
        styles:{
            position:"relative"
        }
    },body);

    const legendg=createSVGEl("g",{},legend);


    for (let i =0;i<len;i++){
        createSVGEl("rect",{
            y: ((i+1)*2)+(i*12),
        	x:2,
        	height:"10",
        	width:"10",
            styles:{
                fill:colors[i]
            }
        },legendg);

        createSVGEl("text",{
            y:((i+1)*2)+(i*12)+6,
            x:14,
            "alignment-baseline":"middle",
            styles:{
                "font-size":"12px"
            },
            text:names[i]
        },legendg)
    }

    makeDraggable(container,{handle:".legend-body"});
    makeResizable(container);
    return container;
}

function getColorBar(colors,config={}){
    const c = config;
    const len =colors.length;
    const colorPct= colors.map((c,i)=>[Math.floor((i/len)*100)+"%",c]);
    const width =c.width || 120;
    const height = c.height || 45;
    const svg = createSVGEl("svg",{
        height:height,
        width:width,
        styles:{
            position:"relative"
        }
    });
    const g= select(svg).append("g");
    const id = getRandomString();
    if (c.label){
        g.append("text")
        .text(c.label)
        .attr("x",10)
        .attr("alignment-baseline","hanging")
        .style("font-size","12px");

    }

    const grad = g.append('defs')
        .append('linearGradient')
        .attr('id', id)
        .attr('x1', '0%') // bottom
        .attr('y1', '0%')
        .attr('x2', '100%') // to top
        .attr('y2', '0%')
        .attr('spreadMethod', 'pad');
    let  bH=height;
    if (c.label){
        bH-=10
    }
    if (c.range){
        bH-=25;
    }
    
    const bar =g.append('rect')
        .attr('x', c.range?10:0)
        .attr('y', c.label?10:0)
        .attr("width",c.range?width-20:20)
        .attr("height",bH)
        .style('fill', `url(#${id})`);

    colorPct.forEach((d)=> {
        grad.append('stop')
        .attr('offset', d[0])
        .attr('stop-color', d[1])
        .attr('stop-opacity', 1);
    });

    //add bottom axis
    if (c.range){
        const scale = scaleLinear().domain([c.range[0],c.range[1]]).range([0,width-20]);
        const axis = axisBottom(scale)
            .tickFormat((v,i)=>v>=10000?Number.parseFloat(v).toPrecision(2):v)
            .ticks(Math.ceil(width/20));
               
        const axisg=g.append("g")
            .call(axis)
            .attr("transform",`translate(10,${height-25})`)
            .selectAll("text")
            .style("text-anchor", "end")
            .attr("dx", "-.4em")
            .attr("dy", ".4em")
            .attr("transform","rotate(-45)")
    }
    const container = createEl("div",{
        height:height,
        width:width,
        styles:{
            position:"absolute"
        }
    });
    container.append(svg);
    makeDraggable(container);
    return container;
   
}


export {getColorLegend,getColorBar}
