import SVGChart from "./SVGChart.js";
import BaseChart from "./BaseChart.js";
import { forceSimulation,forceLink,forceManyBody,forceCenter,drag,scaleSqrt,scaleLinear,schemeReds} from "d3";
import {getColorLegendCustom} from "../utilities/Color.js";

const color_schemes= {
    "blue yellow":["#115f9a", "#1984c5", "#22a7f0", "#48b5c4", "#76c68f", "#a6d75b", "#c9e52f", "#d0ee11", "#d0f400"],
    "red":schemeReds[8].slice(1),
    "blue yellow red":["blue","yellow","red"]


}

class CellNetworkChart extends SVGChart{
    constructor(dataStore,div,config){ 
           
        if (!config.title){
            config.title = config.category;
        }
        config.link_strength = config.link_strength || 1;
        config.node_repulsion = config.node_repulsion || -500;
		super(dataStore,div,config,{});
        const c = this.config;
        const self = this;

        //thickness
        this.linkThicknessScale= scaleLinear();
        c.scales={};
        if (!c.link_thickness){
            c.link_thickness={
                show_legend:true,
                legend_position:[20,80]
            };
            const tq= this.dataStore.getColumnQuantile(c.param[4],0.01);
            this._changeLinkThicknessScale(tq,[1,10]);
        }
        else{
            this._changeLinkThicknessScale();
        }


        //length
        this.linkLengthScale= scaleLinear();
        this.forceLink = forceLink().id(function(d) { return d.id; });
        if (!c.link_length){
            c.link_length={};
            const tq= this.dataStore.getColumnQuantile(c.param[3],0.01);
            this._changeLinkLengthScale(tq,[1,100]);
        }
        else{
            this._changeLinkLengthScale();
        }

        //node size
        this.nodeScale= scaleSqrt();
        if (!c.node_size){
            c.node_size={
                show_legend:true,
                legend_position:[20,180]
            };
            const tq= this.dataStore.getMinMaxForColumn(c.param[6]);
            this._changeNodeSizeScale(tq,[1,15]);
        }
        else{
            this._changeNodeSizeScale();
        }

        if (!c.link_color){
            c.link_color={};
            const tq= this.dataStore.getColumnQuantile(c.param[5],0.01);
            this._changeLinkColorScale("red",tq);

        }
        else{
            this._changeLinkColorScale();
        }
        if (!c.color_legend){
            c.color_legend={display:true}
        }

        c.node_color = c.node_color || "cells"






        c.levels= c.levels || 1;
        c.label_size= c.label_size || 13;

        this.extra_legends= ["linkThicknessLegend","nodeColorLegend"];

      
       
        
      
        this.forceManyBody  = forceManyBody().distanceMax(150);
        this.simulation = forceSimulation()
            .force("link", this.forceLink.strength(this.config.link_strength))
            .force("charge", this.forceManyBody.strength(config.node_repulsion))
            .force("center", forceCenter(this.width/2, this.height/2));
        if (this.config.base_cell){
            this.getLinks2();
        }
        else{
            this.getLinks();
        }
     
       this.setColorLegend();
       this.reCalculate();

	}

    remove(){
        this.simulation.on("tick",null);
        super.remove();
    }

    drawChart(){
        const self = this;
        const c = this.config;
        this.svg.selectAll("*").remove();
        
        const colors  = this.dataStore.getColumnColors(c.param[1]);
        const ncolor = this.dataStore.getColorForCategory(c.param[0],c.category);

        const thick = this.dataStore.getRawColumn(c.param[4]);
    

        const lnks = this.svg.append("g")
            .attr("class", "links")
            .selectAll("line")
            .data(this.linkData)
            .enter();
        this.link = lnks.append("line")
            .attr("stroke-width", d=>{
                return self.linkThicknessScale(thick[d.d_index]);
            })
            .style("stroke",d=>{
                if (self.linkColorScale){
                    return self.linkColorScale(d.d_index)
                }
                return "currentcolor"
            })
            .style("opacity",0.6)
            .on("click",(e,d)=>{
                self.dataStore.dataHighlighted([d.d_index],self);
            })

        this.arrows = lnks.append("path")
        .attr("stroke-width",3)
        .attr("fill","none")
        .style("opacity",c.show_directionality?1.0:0.0)
          .attr("stroke",d=>{
            if (self.linkColorScale){
                return self.linkColorScale(d.d_index);
            }
            return "currentcolor"
        });

  
        this.node = this.svg.append("g")
            .attr("class", "nodes")
            .selectAll("g")
            .data(this.nodeData)
            .enter().append("g")
    
        this.circles = this.node.append("circle")
            .attr("r", d=>{
                return self.nodeScale(d.node_size)
            })
            .attr("fill", function(d) {
                if (c.node_color==="cells"){
                    return colors[d.d_index] 
                }
                return ncolor;
                
            })
            .call(this.setUpDrag());

        this.lables = this.node.append("text")
            .text(function(d) {
                return d.id;
            })
            .style("fill","currentcolor")
            .style("font-size",`${self.config.label_size}px`)
            .attr('x', 6)
            .attr('y', 3);
          
          

        this.simulation
            .nodes(this.nodeData)
            .on("tick", ()=>self.ticked());

            this.simulation.force("link")
            .links(this.linkData);
    }


    _changeLabelSize(s){ 
        const c = this.config;      
        c.label_size=s;
        this.lables.style("font-size",`${s}px`);
    }

    _changeNodeSizeScale(domain,range,update){
        const c = this.config;
        if (range){
            c.node_size.range=range;
        }
        if (domain){
            c.node_size.domain=domain;
        }   
        const self = this;
        this.nodeScale= scaleSqrt().domain(c.node_size.domain)
                    .range(c.node_size.range).
                    clamp(true);

        const l = this.nodeColorLegend;
       
        if (l){
            c.node_size.legend_position=[l.offsetLeft,l.offsetTop];
            l.remove();
            
            
        }
        if (!c.node_size.show_legend){
            delete this.nodeColorLegend;
            return;
        }
        const name = this.dataStore.columnIndex[c.param[6]].name
        let pos = c.node_size.legend_position;
        this.nodeColorLegend=getColorLegendCustom(this.nodeScale,{label:name,type:"circle"});
        this.contentDiv.append(this.nodeColorLegend);
        this.nodeColorLegend.style.top= pos[1]+"px";
        this.nodeColorLegend.style.left= pos[0]+"px";
        if (update){
            this.circles.attr("r", d=>{
                return self.nodeScale(d.node_size);
            })
        }
    }

    _changeLinkThicknessScale(domain,range,update){
        const c = this.config;
        if (range){
            c.link_thickness.range=range;
        }
        if (domain){
            c.link_thickness.domain=domain;
        }
        
        this.linkThicknessScale.domain(c.link_thickness.domain)
                            .range(c.link_thickness.range)
                            .clamp(true);
        
        if (update){
            const self = this;
            const thick = this.dataStore.getRawColumn(c.param[4]);
            this.link.attr("stroke-width", d=>{
                return self.linkThicknessScale(thick[d.d_index]);
            })
        }
        const l = this.linkThicknessLegend;
       
        if (l){
            c.link_thickness.legend_position=[l.offsetLeft,l.offsetTop];
            l.remove();
            
            
        }
      
        if (!c.link_thickness.show_legend){
            delete this.linkThicknessLegend;
            return;
        }
        let pos = c.link_thickness.legend_position;
        const name = this.dataStore.columnIndex[c.param[4]].name

        this.linkThicknessLegend=getColorLegendCustom(this.linkThicknessScale,{label:name,type:"line"});
        this.contentDiv.append(this.linkThicknessLegend);
        this.linkThicknessLegend.style.top= pos[1]+"px";
        this.linkThicknessLegend.style.left= pos[0]+"px";
    }

    _changeLinkLengthScale(domain,range,update){
        const c = this.config;
      
        const self = this;
        if (range){
            c.link_length.range=range;
        }
        if (domain){
            c.link_length.domain=domain;
        }
        //const fudge = Math.max(this.height,this.width)/300;
        const r= [c.link_length.range[0],c.link_length.range[1]];
        this.linkLengthScale.domain([...c.link_length.domain].reverse()).range(r).clamp(true);
      
        const len = this.dataStore.getRawColumn(c.param[3]);
        const maxLen = this.linkLengthScale.range()[1];
           
        this.forceLink.distance(d=>{
                const l = len[d.d_index];
                if (isNaN(l)){
                    return maxLen;
                }
                return self.linkLengthScale(l);
             
        });
        if (update){
            this.simulation.alphaTarget(0.2).restart();
            setTimeout(()=>{
                this.simulation.alphaTarget(0);
            },1000)
        }
    }

    _changeNodeColor(v){
        const c = this.config;
        c.node_color=v;
        const colors  = this.dataStore.getColumnColors(c.param[1]);
        const ncolor = this.dataStore.getColorForCategory(c.param[0],c.category)
        this.circles 
        .attr("fill", function(d) {
            if (c.node_color==="cells"){
                return colors[d.d_index] 
            }
            return ncolor;
            
        });
    }
    

    getConfig(){
        const config = super.getConfig();
        const sl = this.linkThicknessLegend;
        if (sl){
            config.link_thickness.legend_position=[sl.offsetLeft,sl.offsetTop];
        }
        const ns = this.nodeColorLegend;
        if (ns){
            config.node_size.legend_position=[ns.offsetLeft,ns.offsetTop];
        }
        return config;    
    }


    _changeLinkColorScale(colors,domain,update){
        if (domain){
            this.config.link_color.domain=domain;
        }
        if (colors){
            this.config.link_color.colors=colors;
        }
        const d =  this.config.link_color.domain;
        this.linkColorScale= this.dataStore.getColorFunction(this.config.param[5],{
            overideValues:{
                min:d[0],
                max:d[1],
                colors:color_schemes[this.config.link_color.colors]
                //schemeReds[8].slice(1)    
            }
        });
        const self = this;
        if (update){
            this.link.style("stroke",d=>{
                if (self.linkColorScale){
                    return self.linkColorScale(d.d_index)
                    }
                    return "black"
                })
                this.setColorLegend();
            this.arrows
                .attr("stroke",d=>{
                  if (self.linkColorScale){
                      return self.linkColorScale(d.d_index);
                  }
                  return "black"
              });
        }

    }

    onDataFiltered(dim){
        this.reCalculate();  
	}

    getColorLegend(){
        const c = this.config;
        const d =  c.link_color.domain;
        return this.dataStore.getColorLegend(c.param[5],
            {overideValues:{
                min:d[0],
                max:d[1],
                colors:color_schemes[c.link_color.colors]
         
            }});
    }


	colorByDefault(){
        delete this.config.color_by;
        if (this.legend){
            this.legend.remove();
        }
        this.linkColorScale = null;
        this._changeLinkColor()
	}

   

    ticked(){
        const self = this;

        this.lables.attr("text-anchor",d=>{
            return d.x>this.width*2/3?"end":"start";
        }).attr("x",d=>{
            return d.x>this.width*2/3?-6:6;
        })
       

        this.node
        .attr("transform", function(d) {
            d.x = Math.min(d.x,self.width-20);       
            d.x = Math.max(d.x,15);
            d.y = Math.min(d.y,self.height-20);       
            d.y = Math.max(d.y,15);
        
      
          return "translate(" + d.x + "," + d.y + ")";
        })
        this.link
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

        this.arrows.attr("d",d=>{
      
            let len = Math.sqrt(Math.pow(d.source.y-d.target.y,2)+Math.pow(d.target.x-d.source.x,2));
            len = len===0?0.1:len;
            const f = ((len/2)-5)/len;
            const x2 = ((d.target.x)- (d.source.x))/2 + (d.source.x);
            const y2 = ((d.target.y)- (d.source.y))/2 + (d.source.y);
            const xmido = ((d.target.x)- (d.source.x))*f + (d.source.x);
            const ymido = ((d.target.y)- (d.source.y))*f + (d.source.y);
            const xo=(5/len)*(d.source.y-d.target.y)
            const x1=  xmido + xo
            const x3 = xmido - xo;
            const yo=(5/len)*(d.target.x-d.source.x)
            const y1=  ymido + yo
            const y3 = ymido - yo;

            return `M ${x1} ${y1} ${x2} ${y2} ${x3} ${y3}`;       
        });
    }

    setUpDrag(){
            const self = this;
            function dragstarted(event) {
              if (!event.active) self.simulation.alphaTarget(0.3).restart();
              event.subject.fx = event.subject.x;
              event.subject.fy = event.subject.y;
            }
            
            function dragged(event) {
              event.subject.fx = event.x;
              event.subject.fy = event.y;
            }
            
            function dragended(event) {
              if (!event.active) self.simulation.alphaTarget(0);
              event.subject.fx = null;
              event.subject.fy = null;
            }
            
            return drag()
              .on("start", dragstarted)
              .on("drag", dragged)
              .on("end", dragended);
         
    }

    reCalculate(){
        if (this.config.base_cell){
            this.getLinks2();
        }
        else{
            this.getLinks();
        }
        
        this.drawChart();
        this.simulation.alphaTarget(0.2).restart();
        setTimeout(()=>{
            this.simulation.alphaTarget(0);
        },1000)

    }

    changeParameter(){
        const c = this.config;
  
        this.forceLink.strength(c.link_strength);
        this.forceManyBody.strength(c.node_repulsion);
        this.simulation.alphaTarget(0.2).restart();
        setTimeout(()=>{
            this.simulation.alphaTarget(0);
        },1000)


    }

    getLinks2(){
        const p = this.config.param
        const cells = this.dataStore.getColumnValues(p[1]);
        const cells2 = this.dataStore.getColumnValues(p[2]);
        const cat_needed= this.dataStore.getColumnValues(p[0]).indexOf(this.config.category);
        const index = this.dataStore.columnIndex;
        const cat= index[p[0]].data;
        const c1= index[p[1]].data;
        const c2= index[p[2]].data;
        const ns = index[p[6]].data;
     
        this.nodeData=[];
        let firstRound=new Set();
        
        const bc =cells.indexOf(this.config.base_cell);
        firstRound.add(bc);
        const linksNeeded= new Set();
        let allowed= new Set();
        const so =this.config.specific_only;
        if (so){
           allowed= new Set(so[0].map(x=>cells.indexOf(x)))
        }
        else{
            for (let n=0;n<cells.length;n++){
                allowed.add(n)
            }
        }
        const f = this.dataStore.filterArray;
        //get initial Links
        //just connections between center 
        for (let n=0;n<this.dataStore.size;n++){
            if (f[n]>0){
                continue;
            }
            if (cat[n] === cat_needed  && c1[n]==bc && allowed.has(c2[n]) && c1[n] !== c2[n]){
                linksNeeded.add(`${c1[n]}|${c2[n]}`);     
                    firstRound.add(c2[n]);          
            }
        }
        let needed = new Set(firstRound);
        if (this.config.levels==2){
            for (let n=0;n<this.dataStore.size;n++){
                if (f[n]>0){
                    continue;
                }
                if (cat[n] === cat_needed  && firstRound.has(c1[n]) && allowed.has(c2[n]) && c1[n] !== c2[n]){
                    if (c2[n]===bc){
                        continue;
                    }
                    if (!linksNeeded.has(`${c2[n]}|${c1[n]}`)){
                        linksNeeded.add(`${c1[n]}|${c2[n]}`); 
                        needed.add(c2[n]);
                    }
                }
            }
        }
      
       
       for (let n1=0;n1<cells.length;n1++){
            for (let n2=0;n2<cells2.length;n2++){
                if (needed && (!(needed.has(n1))||!needed.has(n2))){
                    continue
                }
               // if (!linksNeeded.has(`${n2}|${n1}`)){
                    linksNeeded.add(`${n1}|${n2}`);
              //  }
            }
        }

        this.linkData= [];
        const ns_map = {};
      
        for (let n=0;n<this.dataStore.size;n++){
            if (f[n]>0){
                continue;
            }
            ns_map[cells[c1[n]]]=ns[n];
            const l = `${c1[n]}|${c2[n]}`;
            

            if (cat[n] === cat_needed && linksNeeded.has(l)){
                linksNeeded.delete(`${c2[n]}|${c1[n]}`)
               
                this.linkData.push({
                    source:cells[c1[n]],
                    target:cells2[c2[n]],
                    d_index:n
                });
            }
        }  
        for (let n=0;n<cells.length;n++){
            if (needed && !(needed.has(n))){
                continue;
            }
            this.nodeData.push({id:cells[n],d_index:n,node_size:ns_map[cells[n]]});
        }
    }

  
    getLinks(){
        const p = this.config.param
        const cells = this.dataStore.getColumnValues(p[1]);
        const cells2 = this.dataStore.getColumnValues(p[2]);
        this.nodeData=[];
        let needed=null;
        const so =this.config.specific_only;
        if (so){
            needed= new Set(so[0].map(x=>cells.indexOf(x)))
        }
       
        const linksNeeded= new Set()
        for (let n1=0;n1<cells.length;n1++){
            for (let n2=n1+1;n2<cells2.length;n2++){
                if (needed && (!(needed.has(n1))||!needed.has(n2))){
                    continue
                }
                linksNeeded.add(`${n1}|${n2}`);
            }
        }
        const cat_needed= this.dataStore.getColumnValues(p[0]).indexOf(this.config.category);
        const index = this.dataStore.columnIndex;
        const cat= index[p[0]].data;
        const c1= index[p[1]].data;
        const c2= index[p[2]].data;
        const ns = index[p[6]].data;
        this.linkData= [];
        const ns_map = {};
      
        const f = this.dataStore.filterArray;
        for (let n=0;n<this.dataStore.size;n++){
            if (f[n]>0){
                continue;
            }
            ns_map[cells[c1[n]]]=ns[n];
            const l = `${c1[n]}|${c2[n]}`;
            if (cat[n] === cat_needed && linksNeeded.has(l)){       
                this.linkData.push({
                    source:cells[c1[n]],
                    target:cells2[c2[n]],
                    d_index:n
                });      
            }
        }

     
        for (let n=0;n<cells.length;n++){
            if (needed && !(needed.has(n))){
                continue;
            }
            this.nodeData.push({id:cells[n],d_index:n,node_size:ns_map[cells[n]]});
        }

    }

  

    setSize(x,y){
        super.setSize(x,y);
        this.simulation.force("center", forceCenter(this.width/2, this.height/2));
        this._changeLinkLengthScale(null,null,true);
    }

    getSettings(){

   
        const settings= super.getSettings();
        const c = this.config
        const vals = this.dataStore.getColumnValues(c.param[0]);
        const cells=  this.dataStore.getColumnValues(c.param[1]).map(x=>{
            return {t:x}
        });
        cells.unshift({t:"None"});
        const li  = vals.map((x)=>{
            return {t:x,v:x}
        });

      
        settings.push({
            
                type:"slider",
                max:1,
                min:0,
                step:0.01,
                doc:this.__doc__,
                current_value:c.stat_cutoff,
                label:"significance cut off",
                func:(x)=>{
                    c.stat_cutoff=x
                    this.reCalculate();
            }
        });
        settings.push({
            
            type:"slider",
            max:20,
            min:5,
            doc:this.__doc__,
            current_value:c.label_size,
            label:"Label Size",
            func:(x)=>{
               this._changeLabelSize(x);
        }
    });


        
            settings.push({
            
                type:"dropdown",
                current_value:c.base_cell || "None",
                values:[cells,"t","t"],
             
                label:"Center Cell",
                func:(x)=>{
                   if (x==="None"){
                    delete c.base_cell;
                   }
                   else{
                    c.base_cell=x;
                   } 
                   this.reCalculate();
            }
        });

        settings.push({
            type:"radiobuttons",
            label:"Number of Levels",
            current_value:c.levels,
            choices:[["1",1],["2",2]],             
            func:(v)=>{
                c.levels=parseInt(v);
                if (c.base_cell){
                    this.reCalculate();
                }
               
               
            }
        });


        const tCols = this.dataStore.getColumnList("number");
   
        let  tQuant= this.dataStore.getColumnQuantile(c.param[5],0.01);

        //link color
        settings.push({
            label:"Link Color",
            type:"dropdown",
            values:[tCols,"name","field"],
            current_value:c.param[5],
            func:(x)=>{
                c.param[5]=x;

            },
            onchange:(controls,v)=>{
                c.param[5]=v;
                const sl =controls["Link Color Domain"].noUiSlider;
                const tq= this.dataStore.getColumnQuantile(v,0.01);
                sl.updateOptions({range:{max:tq[1],min:tq[0]}});
                sl.set(tq);
                this._changeLinkColorScale(null,[tq[0],tq[1]],true);
            }
        });
        const colors=[]
        for (let c in color_schemes){
            colors.push([c,c]);
        }
        settings.push({
            type:"radiobuttons",
            label:"Link Color",
            current_value:c.link_color.colors,
            choices:colors,             
            func:(v)=>{
               this._changeLinkColorScale(v,null,true);   
            }
        });

    
        settings.push({
                
            type:"doubleslider",
            max:tQuant[1],
            min:tQuant[0],
            doc:this.__doc__,
            current_value:this.config.link_color.domain,
            label:"Link Color Domain",
            func:(x,y)=>{       
                this._changeLinkColorScale(null,[x,y],true);
            }
        });


        //link thickness
        tQuant= this.dataStore.getColumnQuantile(c.param[4],0.01);
        settings.push({
            label:"Link Thickness",
            type:"dropdown",
            values:[tCols,"name","field"],
            current_value:c.param[3],
            func:(x)=>{
                c.param[5]=x;

            },
            onchange:(controls,v)=>{
                c.param[4]=v;
                const sl =controls["Link Thickness Domain"].noUiSlider;
                const tq= this.dataStore.getColumnQuantile(v,0.01);
                sl.updateOptions({range:{max:tq[1],min:tq[0]}});
                sl.set(tq);
                this._changeLinkThicknessScale(tq,null,true);
            }
        });

        settings.push({
                
            type:"doubleslider",
            max:20,
            min:1,
            doc:this.__doc__,
            current_value:this.config.link_thickness.range,
            label:"Link Thickness Range",
            func:(x,y)=>{       
                this._changeLinkThicknessScale(null,[x,y],true);
            }
        });
        settings.push({
                
            type:"doubleslider",
            max:tQuant[1],
            min:tQuant[0],
            doc:this.__doc__,
            current_value:this.config.link_thickness.domain,
            label:"Link Thickness Domain",
            func:(x,y)=>{       
                this._changeLinkThicknessScale([x,y],null,true);
            }
        });

        settings.push({
                
            type:"check",
            current_value:c.link_thickness.show_legend,
            label:"Show Link Thickness Legend",
            func:(x)=>{ 
                c.link_thickness.show_legend=x;    
                this._changeLinkThicknessScale();
            }
        });


        //link length
        tQuant= this.dataStore.getColumnQuantile(c.param[3],0.01);
        settings.push({
            label:"Link Length",
            type:"dropdown",
            values:[tCols,"name","field"],
            current_value:c.param[3],
            func:(x)=>{
                c.param[3]=x;

            },
            onchange:(controls,v)=>{
                c.param[3]=v;
                const sl =controls["Link Length Domain"].noUiSlider;
                const tq= this.dataStore.getColumnQuantile(v,0.01);
                sl.updateOptions({range:{max:tq[1],min:tq[0]}});
                sl.set(tq);
                this._changeLinkLengthScale(tq,null,true);
            }
        });

        settings.push({
                
            type:"doubleslider",
            max:200,
            min:1,
            doc:this.__doc__,
            current_value:this.config.link_length.range,
            label:"Link Length Range",
            func:(x,y)=>{       
                this._changeLinkLengthScale(null,[x,y],true);
            }
        });
        settings.push({
                
            type:"doubleslider",
            max:tQuant[1],
            min:tQuant[0],
            doc:this.__doc__,
            current_value:this.config.link_length.domain,
            label:"Link Length Domain",
            func:(x,y)=>{       
                this._changeLinkLengthScale([x,y],null,true);
            }
        });
        settings.push({
                
            type:"check",
            current_value:c.show_directionality,
            label:"Show Directionality",
            func:(x)=>{ 
                c.show_directionality=x;    
                this.arrows.style("opacity",x?1.0:0.0);
            }
        });


        //node  color
      
        settings.push({
            type:"radiobuttons",
            label:"Node Color",
            current_value:c.node_color,
            choices:[["Cells","cells"],["State","state"]],             
            func:(v)=>{
               this._changeNodeColor(v)
                     
            }
        });



        //node size
        tQuant= this.dataStore.getMinMaxForColumn(c.param[6]);
        settings.push({
                
            type:"doubleslider",
            max:50,
            min:1,
            doc:this.__doc__,
            current_value:this.config.node_size.range,
            label:"Node Size Range",
            func:(x,y)=>{       
                this._changeNodeSizeScale(null,[x,y],true);
            }
        });
        settings.push({
                
            type:"doubleslider",
            max:tQuant[1],
            min:tQuant[0],
            doc:this.__doc__,
            current_value:this.config.node_size.domain,
            label:"Node Size Domain",
            func:(x,y)=>{       
                this._changeNodeSizeScale([x,y],null,true);
            }
        });
        settings.push({
                
            type:"check",
            current_value:c.node_size.show_legend,
            label:"Show Node Size Legend",
            func:(x)=>{ 
                c.node_size.show_legend=x;    
                this._changeNodeSizeScale();
            }
        });
    
       
        settings.push({
                
                type:"slider",
                max:2,
                min:0,
                doc:this.__doc__,
                current_value:c.link_strength,
                label:"Link Strength",
                func:(x)=>{
                    c.link_strength=x
                    this.changeParameter();
            }
        });
        settings.push({
                
            type:"slider",
            max:0,
            min:-1000,
            doc:this.__doc__,
            current_value:c.node_repulsion,
            label:"Node Repulsion",
            func:(x)=>{
                c.node_repulsion=x
                this.changeParameter();
        }
    })
        return settings;
    }
}

BaseChart.types["cell_network_chart"]={
    name:"Cell Network Chart",
    allow_user_add:false,
    class:CellNetworkChart,
    params:[
        {
            type:"text",
            name:"category"  //0
        },
        {
            type:"text",
            name:"Cells 1"  //1
        },
        {
            type:"text",
            name:"Cells 2"  //2
        },
        {
            type:"number",
            name:"link length" //3
        },
        {
            type:"number",
            name:"link thickness" //4
        },
        {
            type:"number",
            name:"link color"  //5
        },
        {
            type:"number",
            name:"cell Number"  //6
        }
    ]
}

export default CellNetworkChart;