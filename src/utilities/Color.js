import {createEl,createSVGEl, makeResizable,makeDraggable} from "./Elements.js";
import {select} from "d3-selection";
import {scaleLinear} from "d3-scale";
import {axisBottom} from "d3-axis";
import {getRandomString} from "./Utilities.js";



function getColorLegendCustom(scale,config={}){
    const ticks = scale.ticks(config.ticks || 4);
    const widths= ticks.map(x=>scale(x));
    return getColorLegend(widths,ticks,config);
}

function getColorLegend(colors,names,config={}){
    const len = colors.length;
    const type = config.type || "color";
    const h_fac= type==="circle"?20:12;
    const height=(len*h_fac)+ (len+2)*2;
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
            whiteSpace:"nowrap"

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
    const t_offset= type==="color"?14:type==="line"?20:30;

    
    for (let i =0;i<len;i++){
        if (type === "color"){
            createSVGEl("rect",{
                y: ((i+1)*2)+(i*h_fac),
                x:2,
                height:"10",
                width:"10",
                styles:{
                    fill:colors[i]
                }
            },legendg);
        }
        else if (type==="line"){
            createSVGEl("rect",{
                y: ((i+1)*2)+(i*h_fac),
                x:2,
                height:"10",
                width:colors[i],
                styles:{
                    fill:"gray"
                }
            },legendg);

        }
        else if (type==="circle"){
            createSVGEl("circle",{
                cy: ((i+1)*2)+(i*h_fac)+4,
                cx:15,
                r:colors[i],
                styles:{
                    fill:"gray"
                }
            },legendg);

        }

        const t= createSVGEl("text",{
            y:((i+1)*2)+(i*h_fac)+6,
            x:t_offset,
            "alignment-baseline":"middle",
            //"fill":"currentColor",
            styles:{
                "font-size":"12px",
                "fill":"currentcolor"
            },
            text:names[i]===""?"none":names[i]
        },legendg);
        select(t).style("fill","currentcolor");
        
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
        .style("fill","currentColor")
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


export {getColorLegend,getColorBar,getColorLegendCustom}
