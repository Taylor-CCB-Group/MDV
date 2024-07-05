import {select} from "d3-selection";
import {easeLinear} from "d3-ease";
import BaseChart  from "./BaseChart.js";
import SVGChart from "./SVGChart.js";

class SingleSeriesChart extends SVGChart{
    constructor(dataStore,div,config){
        if (!config.axis){
            config.axis={
                x:{label:dataStore.getColumnName(config.param[4])},
                y:{label:dataStore.getColumnName(config.param[0])}
            }
        }    
	    super(dataStore,div,config,{x:{type:"band"},y:{}});
        if (!config.title){
            this._setTitle();
        }     
        const p= this.config.param;
        const c= this.config;
        const vals= this.dataStore.getColumnValues(p[0]);
        this.x_scale.domain(vals.slice(0));
        this.config.scale= this.config.scale || this.dataStore.getColumnQuantile(p[4],"0.001");
        this.initialSetUp();
        this.drawChart();
        this.addToolTip();       
    }

    remove(notify=true){
        super.remove();
    }

    onDataHighlighted(data){ 
        const i = data.indexes[0];
        const s =  this.graph_area.selectAll("rect");
        s.attr("stroke","white").style("stroke-width", "0px");
        s.filter((d)=>{
            if (d && d[3]===i){
                return true
            }
            return false;
        }).attr("stroke","black").style("stroke-width", "4px");
        this.highlight=i;    
    }

    _setTitle(){
        const a = this.config.categories;
        const v = this.dataStore.getColumnName(this.config.param[4]);
        this.setTitle(`${v} ${a[0]} x ${a[1]}`);
    }

    initialSetUp(){
        const p = this.config.param;
        //interactions
        const c1 =this.dataStore.columnIndex[p[1]].data;
        const c2= this.dataStore.columnIndex[p[2]].data;
        //groups
        const fc = this.dataStore.columnIndex[p[0]].data;
        //samples
        const sc = this.dataStore.columnIndex[p[3]].data;
        const fcVals= this.dataStore.columnIndex[p[0]].values;
        const errorBars= fcVals.map(x=>{
            return {
                values:[],
                totals:0
            };
        })
        const values= this.dataStore.columnIndex[p[4]].data;
        const len = this.dataStore.size;
        const data= [];     
        const cat1 = this.dataStore.columnIndex[p[1]].values.indexOf(this.config.categories[0]);
        const cat2 = this.dataStore.columnIndex[p[2]].values.indexOf(this.config.categories[1]);
        let min = null
        let max=  null;
        for (let i =0;i<len;i++){
            if (c1[i]===cat1 && c2[i]===cat2){
                if (isNaN(values[i])){
                    continue
                }
                max = max===null?values[i]:Math.max(values[i],max);
                min = min==null?values[i]:Math.min(values[i],min);
                
                data.push([values[i],fc[i],sc[i],i,Math.floor(Math.random()*6)+1]);
                const eb = errorBars[fc[i]];
                eb.values.push(values[i]);  
                eb.totals+=values[i];               
            }
        }
        for (const eb of errorBars){
            eb.mean = eb.totals/eb.values.length;
            let n = eb.values.length-1;
            n= n===0?1:n;
            eb.std= eb.values.reduce((prev,cur)=>prev+(Math.pow(cur-eb.mean,2)),0);
            eb.std= Math.sqrt(eb.std/n);
            delete eb.values
        }
        this.errorBars=errorBars;
        const margin = (max-min)/10;
        this.config.scale = [min-margin,max+margin]
        this.data=data;   
    }

    getSettings(){
        const c = this.config;
        const mm= this.dataStore.getColumnQuantile(c.param[4],"0.001")
        const s =[
            //user can alter scale
            {        
                type:"doubleslider",
                max:mm[1],
                min:mm[0],
                doc:this.__doc__,
                current_value:c.scale,
                label:"Scale",
                func:(x,y)=>{     
                    c.scale=[x,y];
                    this.drawChart();
                }
            },
            //user can change interaction pair
            {
                label:this.dataStore.getColumnName(c.param[1]),
                type:"dropdown",
                values:[this.dataStore.getColumnValues(c.param[1],"name_value"),"name","value"],
                current_value:c.categories[0],
                func:(v)=>{      
                   c.categories[0]=v;
                   this.initialSetUp();
                   this.drawChart();        
                }
            },
            {
                label:this.dataStore.getColumnName(c.param[2]),
                type:"dropdown",
                values:[this.dataStore.getColumnValues(c.param[2],"name_value"),"name","value"],
                current_value:c.categories[1],
                func:(v)=>{      
                   c.categories[1]=v;
                   this.initialSetUp();
                   this.drawChart();        
                }
            }
        ];
        return super.getSettings().concat(s);
    }


    drawChart(tTime=300){
        const p = this.config.param;
        this.y_scale.domain([this.config.scale[1],this.config.scale[0]]).clamp(true);   
        this.updateAxis();
        const trans =  select(this.contentDiv).transition()
        .duration(tTime).ease(easeLinear);
        const xvals= this.dataStore.getColumnValues(p[0]);
        const svalues = this.dataStore.getColumnValues(p[3])
        const dim = this._getContentDimensions();
        const recWidth= dim.width/xvals.length;
        const colors = this.dataStore.getColumnColors(p[0]);
        const g= recWidth/4;
        const g1= recWidth/5;
        this.graph_area.selectAll("circle")
            .data(this.data)
            .join("circle")
            .on("mouseover mousemove",(e,d)=>{          
                const c1 = xvals[d[1]];
                const sn = svalues[d[2]];
                this.showToolTip(e,`${c1}<br>${sn}<br>${d[0]}`);
            }).
            on("mouseleave",()=>{
                this.hideToolTip();
            })
            .on("click",(e,d)=>{
                this.dataStore.dataHighlighted([d[3]],this);
            })      
  
            .transition(trans)
            .attr("cx",d=>(recWidth/8)+(recWidth*(d[4]/8))+(d[1]*recWidth))
            .attr("cy",d=>this.y_scale(d[0]))
            .attr("r",5)
            .attr("fill",d=>colors[d[1]]);

            this.graph_area.selectAll(".cbp-rect")
            .data(this.errorBars)
            .join("path")
            .attr("d",(d,i)=>{
                const x  = (i*recWidth) + (recWidth/2);
                
                const av  = this.y_scale(d.mean);
                const sd1 =this.y_scale(d.mean+d.std);
                const sd2  =this.y_scale(d.mean-d.std)
                return `M${x} ${sd1} L${x} ${sd2}
                        M${x-g} ${av} L${x+g} ${av} 
                        M${x-g1} ${sd1} L${x+g1} ${sd1}     
                        M${x-g1} ${sd2} L${x+g1} ${sd2}`;
            })
            .style("fill","none")
            .style("stroke","black")
            .attr("stroke-width", 1)
            .attr("class","cbp-rect")   
    }
  
    setSize(x,y){
        super.setSize(x,y);
        this.drawChart();
    }
}

BaseChart.types["single_series_chart"]={
    name:"Single Interaction Chart",
    required:ds=>ds.interactions?.is_single_region,
    class:SingleSeriesChart,
    extra_controls:(dataSource)=>{
        const i= dataSource.interactions;
        const c1 = i.interaction_columns[0];
        const c2 =i.interaction_columns[1];
        const ex =new Set([c1,c2,i.pivot_column]);
        const groups =dataSource.getColumnList("text").filter(x=>!ex.has(x.field));
        return [
            {
                type:"dropdown",
                name:"c1",
                label:dataSource.getColumnName(c1),
                values:dataSource.getColumnValues(c1).map(x=>({name:x,value:x}))
            },
            {
                type:"dropdown",
                name:"c2",
                label:dataSource.getColumnName(c2),
                values:dataSource.getColumnValues(c2).map(x=>({name:x,value:x}))
            },
            {
                type:"dropdown",
                name:"data",
                label:"Data to display",
                values:dataSource.getColumnList("number").map(x=>({name:x.name,value:x.field}))
            },
            {
                type:"dropdown",
                name:"group",
                label:"Grouping",
                values:groups.map(x=>({name:x.name,value:x.field}))
            }
            
        ];
    },
    init:(config, ds, ec)=>{
        const i =  ds.interactions;
        const param = ds.interactions.cell_radial_chart;
        config.param=[
            ec.group,
            i.interaction_columns[0],
            i.interaction_columns[1],
            i.pivot_column, 
            ec.data         
        ];
        config.categories=[ec.c1,ec.c2];
    }
  
}

export default SingleSeriesChart;