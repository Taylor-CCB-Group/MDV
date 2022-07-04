import {select} from "d3-selection";
import {easeLinear} from "d3-ease";
import BaseChart  from "./BaseChart.js";
import SVGChart from "./SVGChart.js";
import {line} from "d3-shape";



class SingleSeriesChart extends SVGChart{
    constructor(dataStore,div,config){    
		super(dataStore,div,config,{x:{type:"band"},y:{}});
        if (!config.title){
            this._setTitle();
        }     
        const p= this.config.param;
        const c= this.config;
        let vals= this.dataStore.getColumnValues(p[0]);
        this.x_scale.domain(vals.slice(0)); 
       /*  c.color_scale = c.color_scale ||  {log:false};
        if (!c.color_legend){
            c.color_legend={display:true}
        }
       this.setColorFunction()
       //no filtering one time set up
       */
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
       
        this.setTitle(`${a[0]} x ${a[1]}`);
    }

    changeCategory(cat){
        this.config.category=cat;
        this.initialSetUp();
        this.drawChart();
        this._setTitle();
     
    }

    changeValues(col){
        this.config.param[2]=col;
        this.setColorFunction()
        this.initialSetUp();
        this.drawChart();
        this._setTitle();
    }

    initialSetUp(){
        //samples p[0]
        //cat1/cat2 1 and 2
        //data to disply 3....
       
        const p = this.config.param;
        const c1 =this.dataStore.columnIndex[p[1]].data;
        const c2= this.dataStore.columnIndex[p[2]].data;
       
  
        const fc = this.dataStore.columnIndex[p[0]].data;
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

        const max=this.config.scale[1];
       
        const cat1 = this.dataStore.columnIndex[p[1]].values.indexOf(this.config.categories[0]);
        const cat2 = this.dataStore.columnIndex[p[2]].values.indexOf(this.config.categories[1]);

        for (let i =0;i<len;i++){
            if (c1[i]===cat1 && c2[i]===cat2){
                data.push([values[i],fc[i],sc[i],i,Math.floor(Math.random()*6)+1]);
                const eb = errorBars[fc[i]];
                const m = values[i]>max?max:values[i];
                eb.values.push(values[i]);
                eb.totals+=values[i];

                      
            }
        }

        for (let eb of errorBars){
            eb.mean = eb.totals/eb.values.length;
            let n = eb.values.length-1;
            n= n===0?1:n;
            eb.std= eb.values.reduce((prev,cur)=>prev+(Math.pow(cur-eb.mean,2)),0);
            eb.std= Math.sqrt(eb.std/n);
            delete eb.values
        }
        this.errorBars=errorBars;
        this.data=data;
     
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
        const self = this;
        //const recHeight = dim.height/yvals.length;


        const colors = this.dataStore.getColumnColors(p[0]);
        const g= recWidth/4;
        const g1= recWidth/5;
       
   

      this.graph_area.selectAll("circle")
      .data(this.data)
      .join("circle")
      .on("mouseover mousemove",(e,d)=>{          
        const c1 = xvals[d[1]];
        const sn = svalues[d[2]];
        self.showToolTip(e,`${c1}<br>${sn}<br>${d[0]}`);
    }).
    on("mouseleave",()=>{
        self.hideToolTip();
    })
    .on("click",(e,d)=>{
        self.dataStore.dataHighlighted([d[3]],self);
    })      
  
      .transition(trans)
      .attr("cx",d=>(recWidth/8)+(recWidth*(d[4]/8))+(d[1]*recWidth))
      .attr("cy",d=>self.y_scale(d[0]))
      .attr("r",5)
      .attr("fill",d=>colors[d[1]]);

      this.graph_area.selectAll(".cbp-rect")
      .data(this.errorBars)
      .join("path")
      .attr("d",(d,i)=>{
          const x  = (i*recWidth) + (recWidth/2);
         
          const av  = self.y_scale(d.mean);
          let sd1 =self.y_scale(d.mean+d.std);
          let sd2  =self.y_scale(d.mean-d.std)
          return `M${x} ${sd1} L${x} ${sd2}
                  M${x-g} ${av} L${x+g} ${av} 
                  M${x-g1} ${sd1} L${x+g1} ${sd1}     
                  M${x-g1} ${sd2} L${x+g1} ${sd2}      `;
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

    getColorLegend(){
        const cs = this.config.color_scale;
        const conf={
            overideValues:{
                colorLogScale:cs.log,
                colors:["#313695","#4575B4","#74ADD1","#ABD9E9", "#E0F3F8", "#E0F3F8","#FFFFBF", "#FEE090","#FDAE61","#F46D43" ,"#D73027" ,"#A50026" ]
            },
            name:"Scale"
        };
        const mm =this.config.color_scale.min_max;
        if (mm){
            conf.overideValues.min=mm[0];
            conf.overideValues.max=mm[1];
        }
        return this.dataStore.getColorLegend(this.config.param[2],conf);
    }


 
    



    setColorFunction(){
        const p =this.config.param; 
        const conf = {
            useValue:true,
            overideValues:{
                colorLogScale:this.config.color_scale.log,
                colors:["#313695","#4575B4","#74ADD1","#ABD9E9", "#E0F3F8", "#E0F3F8","#FFFFBF", "#FEE090","#FDAE61","#F46D43" ,"#D73027" ,"#A50026" ],
        
            }
        };
        const mm =this.config.color_scale.min_max;
        if (mm){
            conf.overideValues.min=mm[0];
            conf.overideValues.max=mm[1];
        }
        this.colorFunction=(this.dataStore.getColorFunction(p[2],conf));
        this.setColorLegend(); 
    }
   
}

BaseChart.types["single_series_chart"]={
    name:"Single Series Chart",
    allow_user_add:false,
    class:SingleSeriesChart,
    params:[{
        type:"text",
        name:"Categories on y-axis"
    },
    {
        type:"text",
        name:"categories on x-axis"
    },
    {
        type:"number",
        name:"Value to display"
    },
    {
        type:"text",
        name:"column to filter"
    }



    ]
}

export default SingleSeriesChart;