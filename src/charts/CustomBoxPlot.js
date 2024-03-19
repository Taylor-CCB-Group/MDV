import SVGChart from "./SVGChart.js";
import BaseChart from "./BaseChart.js";
import {select} from "d3-selection";
import {easeLinear} from "d3-ease";

function getYLabel(c){
    return c.denominators?`${c.category} per ${c.denom_unit}`: `% abundance ${c.category}`
}

class CustomBoxPlot extends SVGChart{
    constructor(dataStore,div,config){
        //default values
        if (!config.category){
            config.category = dataStore.getColumnValues(config.param[2])[0];
        }
        if(!config.axis){
            config.axis={
                y:{label:getYLabel(config)},
                x:{size:60,rotate_labels:true}
            }
        }
        config.title = config.title || config.category;     
		super(dataStore,div,config,{x:{type:"band"},y:{}});
        const p= this.config.param;
        const vals= this.dataStore.getColumnValues(p[0]);
        this.x_scale.domain(vals.slice(0));       
        this.addToolTip();
        this.dim = this.dataStore.getDimension("category_dimension");
        this.addDenominators();      
        this.onDataFiltered();      
    }

   

    onDataFiltered(){
        const p =this.config.param;
        let cat =null;
        if (p.length===3){
            cat  = this.dataStore.getColumnValues(this.config.param[2]).indexOf(this.config.category);
        }
        this.dim.getSankeyData(data=>{
            this.data=data;
            this.drawChart();
        },p,{method:"proportion",category:cat,denominators:this.denom});
    }

    remove(notify=true){
        this.dim.destroy(notify);
        super.remove();
    }

    setSize(x,y){
        super.setSize(x,y);
        this.drawChart();
    }

    addDenominators(){
        const d= this.config.denominators;
        this.denom=null;
        if (d){
            this.denom=[];
            const vals= this.dataStore.getColumnValues(this.config.param[1])
            for (let n=0;n<vals.length;n++){
                this.denom.push(d[vals[n]] || 100)
            }
        }
    }

    getSettings(){
        const c = this.config;
        const values= this.dataStore.getColumnValues(c.param[2]).map(x=>({t:x}));
        //allow user to change category
        const s =[
            {
                values:[values,"t","t"],
                current_value:c.category,
                label:"Category",
                type:"dropdown",
                func: x =>{
                    c.category=x;
                    c.axis.y.label = getYLabel(c);
                    this.setTitle(x);
                    this.onDataFiltered();
                }
            }
        ];
        return super.getSettings().concat(s);
    }

    drawChart(){
        const trans =  select(this.contentDiv).transition()
        .duration(200).ease(easeLinear);
        const p = this.config.param;
        this.y_scale.domain([this.data.max+this.data.max/10,0])
        this.updateAxis();
        const colors = this.dataStore.getColumnColors(p[0]);
        const vals= this.dataStore.getColumnValues(p[0]);
        const vals2 =  this.dataStore.getColumnValues(p[1]); 
        const dim  = this._getContentDimensions();
        const recWidth= dim.width/vals.length;
        const self = this;
        const xpos = recWidth/5;
        const margin = recWidth/10;
        const radius = 5;
        const g= recWidth/4;
        const g1= recWidth/5;
        this.graph_area.selectAll(".cbp-col")
        .data(this.data)
        .join("g")
        .attr("transform",(d,i)=>{

            return `translate(${i*recWidth},0)`
        })
        .attr("class","cbp-col")
        .selectAll(".cbp-circle")
        .data(d=>d)
        .join("circle")
        .attr("class","cbp-circle")
        .on("mouseover pointermove",(e,d)=>{    
            const h= `${vals[d[1]]}<br>${vals2[d[2]]}`;             
            self.showToolTip(e,h)
        }).
        on("mouseleave",()=>{
            self.hideToolTip();
        })   
        .transition(trans)
        .attr("cy",(d,i)=>{
            return self.y_scale(d[0]);
        })
        .attr("fill",(d,i)=>colors[d[1]])
        .attr("r",radius)
        .attr("cx",d=>xpos+(margin*d[3]));

        this.graph_area.selectAll(".cbp-rect")
        .data(this.data)
        .join("path")
        .attr("d",(d,i)=>{
            const x  = (i*recWidth) + (recWidth/2);
           
            const av  = self.y_scale(d.av);
            let sd1 = d.av+d.std;
            sd1= self.y_scale(sd1>d.max?d.max:sd1);
            let sd2 = d.av-d.std;
            sd2= self.y_scale(sd2<d.min?d.min:sd2)
            return `M${x} ${sd1} L${x} ${sd2}
                    M${x-g} ${av} L${x+g} ${av} 
                    M${x-g1} ${sd1} L${x+g1} ${sd1}     
                    M${x-g1} ${sd2} L${x+g1} ${sd2}      `;
        })
        .style("fill","none")
        .style("stroke","currentColor")
        .attr("stroke-width", 1)
        .attr("class","cbp-rect")
    }
}

BaseChart.types["custom_box_plot"]={
    name:"Abundance Box Plot",
    params:[
        {type:"text",name:"Sample Group"},
        {type:"text",name:"Sample"},
        {type:"text",name:"Categories"}
    ],
    class:CustomBoxPlot 
}

//special chart where instead of abundance, number per area of region is displayed
//the sample column is already known - the areas (and units) are added
//depending on the regions config in the datastore
BaseChart.types["density_box_plot"]={
    name:"Density Box Plot",
    class:CustomBoxPlot,
    required:["regions"],
    extra_controls:(ds)=>{
        const vs = ds.getColumnList("text").map(x=>({name:x.name,value:x.field}));
        return [
            {
                type:"dropdown",
                name:"group",
                label:"Group",
                values:vs
            },
            {
                type:"dropdown",
                name:"cat",
                label:"Category",
                values:vs
            }
        ]
    },
    init:(config,ds,ec)=>{   
        config.denominators={};
        config.denom_unit= ds.regions.scale_unit+"²"; //should be power 2
        config.param=[ec.group,ds.regions.region_field,ec.cat];
        //work out area of each region
        for (let reg_name in ds.regions.all_regions){
            const reg = ds.regions.all_regions[reg_name];
            const area = (reg.roi.max_x - reg.roi.min_x) * (reg.roi.max_x - reg.roi.min_x);
            config.denominators[reg_name]=area*(ds.regions.scale*ds.regions.scale);
        }
    }

}

export default CustomBoxPlot;