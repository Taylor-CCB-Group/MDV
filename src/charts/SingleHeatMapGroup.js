import {select} from "d3-selection";
import {easeLinear} from "d3-ease";
import BaseChart  from "./BaseChart.js";
import SVGChart from "./SVGChart.js";



class SingleHeatMapGroup extends SVGChart{
    constructor(dataStore,div,config){
          
		super(dataStore,div,config,{x:{type:"band"},y:{type:"band"},tx:{size:50}});
        if (!config.title){
            this._setTitle();
        }
        
        const c= this.config;
      
        c.color_scale = c.color_scale ||  {log:false};
        if (!c.color_legend){
            c.color_legend={display:true}
        }
       this.setColorFunction()
       //no filtering one time set up
       this.initialSetUp();
       this.drawChart();
       
         
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
        const a = this.config.category;
        const b = this.dataStore.getColumnName(this.config.param[2]);
        this.setTitle(`${a}-${b}`);
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
        const c = this.config;
        const p = c.param;
      
        const y_col =this.dataStore.columnIndex[p[0]];
        const y = y_col.data;
        const x_col= this.dataStore.columnIndex[p[1]];
        const x = x_col.data;
        const z_col = this.dataStore.columnIndex[p[3]];
        let oset1= null;
        let oset2 = null;
        let y_len = y_col.values.length;
        let x_len = x_col.values.length;
        const z_len = z_col.values.length;
        

        const so = c.specific_only;
        if (so){
            oset1 = so[0].map(x=>y_col.values.indexOf(x));
            oset2 = so[1].map(x=>x_col.values.indexOf(x));
            y_len = so[0].length;
            x_len =  so[1].length;
          
        }
        const sc= [];

        this.y_scale.domain(so?so[0].slice(0):y_col.values.slice(0)); 
        this.x_scale.domain(so?so[1].slice(0):x_col.values.slice(0));   
        
        


        const fc = z_col.data;
        const vc  = this.dataStore.columnIndex[p[2]].data;
        const len = this.dataStore.size;
        const data= new Array(y_len);
        for (let n=0;n<y_len;n++){
            data[n]=new Array(x_len);
            for (let i=0;i<x_len;i++){
                data[n][i]= new Array(z_len);
            }
        }
       
        for (let i =0;i<len;i++){
            
                let py = y[i];
                let px = x[i];
                if (so){
                    px= oset1.indexOf(px);
                    py= oset2.indexOf(py);
                    if (px===-1 || py===-1){
                        continue;
                    }
                }
                const d = data[px][py][fc[i]];
                if (!d){
                    data[px][py][fc[i]]=[y[i],x[i],vc[i],i,1,fc[i]];
                }
                else{
                    d[2]+=vc[i];
                    d[4]++;
                }
               
            
        }
        for (let y1=0;y1<y_len;y1++){
            const new_arr=[]
            for (let x1=0;x1<x_len;x1++){
                for (let z1=0;z1<z_len;z1++){
                    let d= data[y1][x1][z1];
                    if (!d){
                        d= [0,0,0,0,0,0];
                    }
                    else{
                        d[2]/=d[4];
                    }
                    
                    new_arr.push(d);
                }

            }
            data[y1]=new_arr;
        }
        this.data=data;
    }


    drawChart(tTime=300){
        const trans =  select(this.contentDiv).transition()
        .duration(tTime).ease(easeLinear);
        const yvals= this.dataStore.getColumnValues(this.config.param[0]);
        const xvals= this.dataStore.getColumnValues(this.config.param[1]);
        const zvals = this.dataStore.getColumnValues(this.config.param[3]);
        const so = this.config.specific_only
        const y_len = so?so[0].length:yvals.length;
        const z_len =zvals.length;
        const x_len_i = so?so[1].length:xvals.length;
        const x_len= x_len_i*zvals.length;
        const dim = this._getContentDimensions();
        const recWidth= dim.width/x_len;
        
        const grpWidth= dim.width/x_len_i;
        const arecWidth= (grpWidth-4)/z_len;
        const recHeight = dim.height/y_len;
        this.graph_area.selectAll(".heatmap-row")
        .data(this.data)
        .join("g")
        .attr("transform",(d,i)=>`translate(0,${i*recHeight})`)
        .attr("class","heatmap-row")
        .selectAll(".heatmap-rect")
        .data(d =>{
            return d;
        })
        .join("rect")
        .attr("class","heatmap-rect")
        .attr("x",(d,i)=>{
            let x = (Math.floor(i/z_len)*grpWidth)+2
            x+=(i%z_len)*arecWidth;
            return x;
        })
        .attr("height",recHeight)
        .attr("width",arecWidth)
        .on("click",(e,d)=>{
           this.dataStore.dataHighlighted([d[3]],this);
        })
        .transition(trans)
        .attr("fill",(d,i)=>{
            if (!d){
                return "#f4f0ec";
            }
            if (isNaN(d[2])){
                return "#f4f0ec";
            }
            return this.colorFunction(d[2]);
        });

        this.updateAxis();
        this.tx_axis_svg.selectAll(".tx-text")
        .data(this.data[0])
        .join("text")
        .attr("class","tx-text")
        .attr("font-size","10px")
        .text((d,i)=>{
           return zvals[i%zvals.length]
        })
        .style("text-anchor", "end")
       .attr("transform",(d,i)=> ` translate(${(i*recWidth)+(recWidth/2)},${dim.top}) rotate(60)`)


    
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


    getSettings(){
        const settings= super.getSettings();
        const c = this.config
        const vals = this.dataStore.getColumnValues(c.param[3]);
        const li  = vals.map((x)=>{
            return {t:x,v:x}
        })
       
        settings.push({
            label:"Change Category",
            type:"dropdown",
            values:[li,"t","v"],
            current_value:c.category,
            func:(x)=>{
               this.changeCategory(x);
            }
        });


        const params= this.dataStore.getColumnList("number");
        settings.push({
            label:"Change Values",
            type:"dropdown",
            values:[params,"name","field"],
            current_value:c.param[2],
            func:(x)=>{
               this.changeValues(x);
            }
        });
        return settings;
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

BaseChart.types["single_heat_map_group"]={
    name:"Single Heat Map",
    allow_user_add:false,
    class:SingleHeatMapGroup,
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

export default SingleHeatMapGroup;