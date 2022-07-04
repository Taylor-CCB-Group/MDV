import {select} from "d3-selection";
import {easeLinear} from "d3-ease";
import SVGChart from "./SVGChart";
import BaseChart from "./BaseChart.js";
import {pie,arc} from 'd3-shape'


class MultiPieChart extends SVGChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config,{});

        this.config.type = "multi_pie_chart";        
        this.pie= pie().sort(null).value(d=>d[2]);
        this.arc= arc();
        //this will draw the chart
        //this.onDataFiltered(null);
        this.initialSetUp();
        this.drawChart();
        this.addToolTip();
	}

    initialSetUp(){
        const p = this.config.param;
        const so = this.config.specific_only;
        const c1 =this.dataStore.columnIndex[p[0]]
        const c1d = c1.data;
        const c2 = this.dataStore.columnIndex[p[1]]
        const c2d = c2.data;
        const c1v= c1.values;
        const fc = this.dataStore.columnIndex[p[3]].data;
        const vc  = this.dataStore.columnIndex[p[2]].data;
        const len = this.dataStore.size;
        let data= new Array(c1v.length);
        for (let n=0;n<c1.values.length;n++){
            data[n]={
                data:new Array(c2.values.length),
                name:c1v[n],
                index:n
            }
        }
        const cat = this.dataStore.columnIndex[p[3]].values.indexOf(this.config.category);
        for (let i =0;i<len;i++){
            if (fc[i]===cat){
                const d = data[c1d[i]].data;
                let val = d[c2d[i]]
                if (val == null){
                    d[c2d[i]]=[c1d[i],c2d[i],vc[i],i,1]
                }
                else{
                    val[2]+=vc[i];
                    val[4]++;
                }
               
            }
        }
        for (let n=0;n<c1.values.length;n++){
            let d= data[n].data;
            if (this.config.hide_self_interaction){
                d.splice(n,1);
            }
            let sum =0;
            for (let p of d){
                if (isNaN(p[2])){
                    p[2]=0;
                }
                else{
                    p[2]=p[2]/p[4];
                }
                
                sum+=p[2]
            }
            data[n].sum=sum;
        }
        if (so){
          data = data.filter(x=>so[0].indexOf(x.name)!==-1);
          const set2 = new Set(so[1].map(x=>c2.values.indexOf(x)))
          for (let n=0;n<data.length;n++){
              data[n].data= data[n].data.filter(x=>set2.has(x[1]));
              data[n].sum=data[n].data.reduce((a,b)=>a+b[2],0)

          }
        }
        
        this.data=data;
      
    }

    setSize(x,y){
        super.setSize(x,y);
        this.drawChart();
    }

    onDataHighlighted(data){ 
            const i = data.indexes[0];
            const s =  this.graph_area.selectAll("path");
           s.attr("stroke","white").style("stroke-width", "2px");
            s.filter((d)=>{
                return d.data[3]===i
            }).attr("stroke","black").style("stroke-width", "4px");
            this.highlight=i;
        
    }

    drawChart(tTime=1000){
        //can't get transitions to work properly
        const trans =  select(this.contentDiv).transition()
        .duration(tTime).ease(easeLinear);
        const dim = this._getContentDimensions();
        const ratio = dim.width/dim.height;
        const rc = Math.ceil(Math.sqrt(this.data.length));
        const pw = dim.width/rc;
       
        const ph = dim.height/rc
        let radius =pw;
        if (ph/1.2<pw){
            radius = ph/1.2;
        }

         radius/=2;

     
        //calculate the radius and center of the plot based on the
        //container dimensions
        
        
      
        this.arc.innerRadius(radius*0.5).outerRadius(radius*0.9)
        let fontSize= Math.floor(pw/10);
        fontSize=fontSize>12?12:fontSize;
        const v1 = this.dataStore.getColumnValues(this.config.param[0]);
        const c1  = this.dataStore.getColumnColors(this.config.param[0]); 
        const v2 = this.dataStore.getColumnValues(this.config.param[1]);
        const c2  = this.dataStore.getColumnColors(this.config.param[1]);  
        
    
        //update the pie chart
        const self = this;
        const pas = this.graph_area.selectAll(".pie-area")
        .data(this.data)
        .join("g")
        .attr("transform",(d,i)=>{
            const y= Math.floor(i/rc) *(ph);
            const x=  i%rc *pw;
            return `translate(${x+pw/2},${y+ph/2})`
            }
        )
        .attr("class","pie-area");

        pas.selectAll(".pie-inner")
        .data((d,i)=>{
            return [d.index];
        })
        .join("circle")
        .attr("class","pie-inner")
        .attr("r",radius*0.4)
        .style("fill",(d,i)=>{    
            return c1[d]    
         });


         pas.selectAll(".pie-text-header")
         .data((d,i)=>{
             return [d.name]
         })
         .join("text")
         .attr("text-anchor", "middle")
         .attr("class","pie-text-header")
         .attr("y",-(ph/2)+(ph/10))
         .attr("font-size",fontSize+"px")
         .attr("dominant-baseline", "bottom") 
        .text(d=>d);


        pas.selectAll(".pie-outer")
        .data(d=>{
            const r =self.pie(d.data.filter(x=>x[2]>0));
            for (let i of r){
                i.sum=d.sum;
            }
            return r;
        })
        .join("path")
        .attr("class","pie-outer")
        .attr("stroke",d=>d.data[3]===self.highlight?"black":"white")
        .attr("d",this.arc)
        .on("click",(e,d)=>{
            self.dataStore.dataHighlighted([d.data[3]],self);
        })
    .on("mouseover mousemove",(e,d)=>{          
        const c1 = v1[d.data[0]];
        const c2 = v2[d.data[1]];
        const p = ((d.data[2]/d.sum)*100).toFixed(2)
        self.showToolTip(e,`${c1}<br>${c2}<br>${d.data[2].toFixed(2)}<br>${p}%`);
    }).
    on("mouseleave",()=>{
        self.hideToolTip();
    })   
        .style("stroke-width",d=>d.data[3]===self.highlight?"4px":"2px")
        .style("fill",(d,i)=>{    
                return c2[d.data[1]]    
        });
        



    


        //update the text
       /* this.graph_area.selectAll(".pie-text") 
        .data(data,(d,i)=>{
            d.data[1]})
        .join("text")
        .on("click",(e,d)=>{
            self.filterCategories(vals[d.data[1]],e.shiftKey);
        })
        .attr("class","pie-text")
        .text((d,i)=>{
            //only show text if segment is big enough
            const angle=(d.endAngle-d.startAngle);
            return angle > 0.4?vals[d.data[1]]:"";
        })
        .attr("x",(d)=>{    
           return  self.arc.centroid(d)[0]
        })
        .attr("font-size",fontSize+"px")
        .attr("y",(d)=>self.arc.centroid(d)[1])
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "central")
        */
    }

    getSettings(){
        const s = super.getSettings();
        const c = this.config;
        return s.concat(
            {
                type:"check",
                label:"Hide self intearction",
                current_value:c.hide_self_interaction,
                func:(x)=>{
                    c.hide_self_interaction=x
                    this.initialSetUp();
                    this.drawChart();
                }
            });
     
    }



}

BaseChart.types["multi_pie_chart"]={
    name:"Multi Pie Chart",
    allow_user_add:false,
    class:MultiPieChart,
    params:[{
        type:"text",
        name:"first cat"
    },
    {
        type:"text",
        name:"second cat"
    },
    {
        type:"number",
        name:"Value to display in each chart"
    },
    {
        type:"text",
        name:"column to filter"
    }



    ]

}

export default MultiPieChart;