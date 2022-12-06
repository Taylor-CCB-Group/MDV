import {select} from "d3-selection";
import {easeLinear} from "d3-ease";
import BaseChart  from "./BaseChart.js";
import SVGChart from "./SVGChart.js";



class SingleHeatMap extends SVGChart{
    constructor(dataStore,div,config){
          
		super(dataStore,div,config,{x:{type:"band"},y:{type:"band"}});
        if (!config.title){
            this._setTitle();
        }      
        const p= this.config.param;
        let vals= this.dataStore.getColumnValues(p[0]);
        this.y_scale.domain(vals.slice(0)); 
        vals=this.dataStore.getColumnValues(p[1]);
        this.x_scale.domain(vals.slice(0));   
        const c= this.config;   
        c.color_scale = c.color_scale ||  {log:false};
        if (!c.color_legend){
            c.color_legend={display:true}
        }
        this.setColorFunction()
        //no filtering one time set up
        this.initialSetUp();
        this.addToolTip();
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
    onDataFiltered(dim){
        this.drawChart();  
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
        const y_col =this.dataStore.columnIndex[p[0]]
        const y = y_col.data;
        const x_col= this.dataStore.columnIndex[p[1]]
        const x = x_col.data;
        let oset1= null;
        let oset2 = null;
        let y_len = y_col.values.length;
        let x_len = x_col.values.length;
        const so = c.specific_only;
        if (so){
            oset1 = so[0].map(x=>y_col.values.indexOf(x));
            oset2 = so[1].map(x=>x_col.values.indexOf(x));
            y_len = so[0].length;
            x_len =  so[1].length;
            this.y_scale.domain(so[0].slice(0)); 
            this.x_scale.domain(so[1].slice(0));   
        }

        const ci = this.dataStore.columnIndex;
        const fc = ci[p[3]].data;
        const vc  = ci[p[2]].data;
     
        
        const len = this.dataStore.size;
        const data= new Array();
        for (let n=0;n<y_len;n++){
            data[n]=new Array(x_len);
        }
        const cat = this.dataStore.columnIndex[p[3]].values.indexOf(c.category);
        let ct1 = null;
        let ct2 = null;
        if (c.show_group_interactions){
            ct1= this.dataStore.getRawColumn(p[4]);
            ct2= this.dataStore.getRawColumn(p[5]);     
        }
      
        for (let i =0;i<len;i++){
            if (fc[i]===cat){
                let py = y[i];
                let px = x[i];
                if (so){
                    px= oset1.indexOf(px);
                    py= oset2.indexOf(py);
                    if (px===-1 || py===-1){
                        continue;
                    }
                }
                let d = data[px][py];
                if (!d){
                    data[px][py]=[y[i],x[i],vc[i],i,1];
                    if (ct1){
                        data[px][py].push(`${ct1[i]}-${ct2[i]}`);
                    }               
                }
                else{
                    d[2]+=vc[i];
                    d[4]++;
                }      
            }
        }
        for (let y1=0;y1<y_len;y1++){
            for (let x1=0;x1<x_len;x1++){
                const d= data[y1][x1];
                d[2]/=d[4];

            }
        }
        this.data=data;
    }


    drawChart(tTime=300){
        const trans =  select(this.contentDiv).transition()
        .duration(tTime).ease(easeLinear);
        const c  =this.config;
        const p = c.param
        const yvals= this.dataStore.getColumnValues(this.config.param[0]);
        const xvals= this.dataStore.getColumnValues(this.config.param[1]);
        const so = this.config.specific_only
        const y_len = so?so[0].length:yvals.length;
        const x_len = so?so[1].length:xvals.length;
        const dim = this._getContentDimensions();
        const recWidth= dim.width/x_len;
        const self = this;
        const recHeight = dim.height/y_len;
        const f  = this.dataStore.filterArray;
        const mcol = "#E0E0E0";//"#f4f0ec"
        let interactions= null;
        if (this.config.show_group_interactions){
            interactions={};
            const gv1 = this.dataStore.getColumnValues(p[4]);
            const gv2= this.dataStore.getColumnValues(p[5]);
            for (let v1 in gv1 ){
                for (let v2 in gv2){
                    interactions[`${v1}-${v2}`]=0;
                }
            }
        }

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
        
        .attr("x",(d,i)=>(i*recWidth)+1)
        .attr("height",recHeight-2)
        .attr("width",recWidth-2)
        .on("click",(e,d)=>{
           self.dataStore.dataHighlighted([d[3]],self);
        })
        .on("mouseover mousemove",(e,d)=>{ 
            const col = xvals[d[0]];
            const row  =yvals[d[1]];
            self.showToolTip(e,`${col}<br>${row}<br>${d[2].toPrecision(3)}`);
        }).
        on("mouseleave",()=>{
            self.hideToolTip();
        })   
        .transition(trans)
        .attr("fill",(d,i)=>{
            if (!d){
                return mcol;
            }
            if (isNaN(d[2])){
                return mcol;
            }
            if (f[d[3]]>0){
                return mcol;
            }
           
            if (interactions){
                interactions[d[5]]++;
            }
            return self.colorFunction(d[2]);
        });
        if (interactions){
            const v1 = this.dataStore.getColumnValues(p[4]);
            const v2= this.dataStore.getColumnValues(p[5]);
            const msg=[];
            for (let i in interactions){
                const arr = i.split("-");
                msg.push(`${v1[arr[0]]}-${v2[arr[1]]}:${interactions[i]}`);
                this.legendIcon.setAttribute("aria-label",msg.join(", "));

            }
        }
    
        this.updateAxis();
        
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

BaseChart.types["single_heat_map"]={
    name:"Single Heat Map",
    allow_user_add:false,
    class:SingleHeatMap,
    params:[
        {
            type:"text",
            name:"Categories on y-axis" //0
        },
        {
            type:"text",
            name:"categories on x-axis" //1
        },
        {
            type:"number",
            name:"Value to display" //2
        },
        {
            type:"text",
            name:"group to display" //3
        },
    
        {
            type:"text",
            name:"group1" //4
        },
        {
            type:"text",
            name:"group2" //5
        }
    ]
};

export default SingleHeatMap;