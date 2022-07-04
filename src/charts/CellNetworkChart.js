import SVGChart from "./SVGChart.js";
import BaseChart from "./BaseChart.js";
import { forceSimulation,forceLink,forceManyBody,forceCenter,drag } from "d3";

class CellNetworkChart extends SVGChart{
    constructor(dataStore,div,config){       
        if (!config.title){
            config.title = config.category;
        }
        config.link_strength = config.link_strength || 1;
        config.node_repulsion = config.node_repulsion || -500;
		super(dataStore,div,config,{});
        this.forceLink = forceLink().id(function(d) { return d.id; }).distance(d=>(d.value)*20);
        this.forceManyBody  = forceManyBody().distanceMax(150);
        this.simulation = forceSimulation()
            .force("link", this.forceLink.strength(this.config.link_strength))
            .force("charge", this.forceManyBody.strength(config.node_repulsion))
            .force("center", forceCenter(this.width/2, this.height/2));
        this.getLinks();
        this.drawChart();
	}

    drawChart(){
        const self = this;
        this.svg.selectAll("*").remove();
        const colors  = this.dataStore.getColumnColors(this.config.param[1]);
        this.link = this.svg.append("g")
            .attr("class", "links")
            .selectAll("line")
            .data(this.linkData)
            .enter().append("line")
            .attr("stroke-width", d=>{
                return (5-d.value)*2;
            })
            .style("stroke","black")
            .style("opacity",0.2)
  
        this.node = this.svg.append("g")
            .attr("class", "nodes")
            .selectAll("g")
            .data(this.nodeData)
            .enter().append("g")
    
        this.circles = this.node.append("circle")
            .attr("r", 5)
            .attr("fill", function(d) { return colors[d.index] })
            .call(this.setUpDrag());

        this.lables = this.node.append("text")
            .text(function(d) {
                return d.id;
            })
            .attr('x', 6)
            .attr('y', 3);
            this.node.append("title")
            .text(function(d) { return d.id; });

        this.simulation
            .nodes(this.nodeData)
            .on("tick", ()=>self.ticked());

            this.simulation.force("link")
            .links(this.linkData);
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
        this.getLinks();
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
        for (let n=0;n<cells.length;n++){
            if (needed && !(needed.has(n))){
                continue;
            }
            this.nodeData.push({id:cells[n],index:n});
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
        const v= index[p[3]].data;
        const s= index[p[4]].data;
        this.linkData= [];

        for (let n=0;n<this.dataStore.size;n++){
            const l = `${c1[n]}|${c2[n]}`;
            if (cat[n] === cat_needed && linksNeeded.has(l) && s[n]<this.config.stat_cutoff){
                let val = 5-v[n];
                val= Math.max(val,0);
            
                this.linkData.push({
                    source:cells[c1[n]],
                    target:cells2[c2[n]],
                    value:val
                });
            }
        }

    }

  

    setSize(x,y){
        super.setSize(x,y);
        this.simulation.force("center", forceCenter(this.width/2, this.height/2));
        this.simulation.alphaTarget(0.2).restart();
        setTimeout(()=>{
            this.simulation.alphaTarget(0);
        },1000)
    }

    getSettings(){
        const settings= super.getSettings();
        const c = this.config
        const vals = this.dataStore.getColumnValues(c.param[0]);
        const li  = vals.map((x)=>{
            return {t:x,v:x}
        })

      
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
            })
        
       
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
    params:[{
        type:"text",
        name:"cat1"
    },
    {
        type:"text",
        name:"cat2"
    },
    {
        type:"number",
        name:"interaction"
    },
    {
        type:"number",
        name:"fdr"
    }



    ]
}



export default CellNetworkChart;