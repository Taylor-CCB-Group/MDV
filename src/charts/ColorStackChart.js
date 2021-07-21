import {DCChart} from "./DCChart.js"

class ColorStackChart extends DCChart{
	constructor(ndx,chart_type,div,config,){
		super(ndx,chart_type,div,config);
		this.init();

		this.setColorFromConfig(); 
        if (this.config.size){
        	this.setSize(this.config.size[0],this.config.size[1]);
        }
        else{
        	this.setSize();
        }
	}
	 setFilter(a,b){
        this.chart.filter(null);
       setTimeout(()=>{
            this.chart.filter(dc.filters.RangedFilter(a,b));
             dc.redrawAll();
     },50);
       //dc.redrawAll();
    }

  

    setField(new_field){
        let self = this;
        this.param=new_field;
         this.dim = this.ndx.dimension(
            function(d){return d[self.param];}
        );
        this.max = this.dim.top(1)[0][this.param];
        this.min= this.dim.bottom(1)[0][this.param];
        this.param=new_field;
        this.setParameters({max:this.max,min:this.min,bin_number:this.bin_number});
    }




    getFilter(){
    	let filters = this.chart.filters();
    	if (filters.length===0){
    		return [];
    	}
    	return [{field:this.config.param,operand:"between",value:[filters[0][0],filters[0][1]]}];
    }

   	/**
   	* params
   	* max
   	* min
   	* bin_number
   	* max_y
   	* color_by
   	*/



	getColorOptions(){
		//let d = this.div.parent();
        //let id = this.div.attr("id")+"-parent";
        //d.attr("id",id);
		return {
			datatype:"text",
			div:"sss"
		}
	}

	setColorFromConfig(){
		let cb = this.config.color_by;
		let params =null;
		if (cb){
			let opt = this.getColorOptions();
			params = FilterPanel.getColorScale(cb.column,this.ndx.getOriginalData(),cb.scheme,opt.div);
			
		}
		this.colorByField(params);
	}
	setColorStack(param){
    	let self =this;
    	if (this.config.color_by){
        	this.categories=[];
        	let cols=[]
        	for (let cat in this.config.color_by.value_to_color){
        		this.categories.push(cat);
        		cols.push(this.config.color_by.value_to_color[cat])

        	}
        	this.categories.push("Other")
        	cols.push("rgb(220, 220, 220)")
        	this.chart.ordinalColors(cols);
        	this.chart.group(this.group,this.categories[0],function(d){
            	let val= d.value[self.categories[0]];
            	if (val===undefined){
            		return 0;
            	}
            	return val;
            });
        	let field =self.config.color_by.column.field;
        	let value_to_color= self.config.color_by.value_to_color;
        	this.group.reduce(
				function(p,v,nf){
					let val = v[field];
					if (!value_to_color[val]){
						val="Other"
					}
					let num = p[val];
					let inc = param?v[param]:1
					inc = inc?inc:0
					if (!num){
						p[val]=inc;
					}
					else{
						p[val]+=inc;
					}
					return p;
				},
				function(p,v,nf){
					let val = v[field];
					if (!value_to_color[val]){
						val="Other"
					}

					let num = p[val];
					if (num){
						let inc = param?v[param]:1
						inc = inc?inc:0
						p[val]-=inc;
					}

					return p;
				},
				function(){
					return {};
				});
			for (let i=1;i<this.categories.length;i++){
        		this.chart.stack(this.group,this.categories[i],this._getColorFunction(i));
        	}
              
        }
        else{
        	this.chart.valueAccessor(function(d){return d.value});
        	this.chart.group(this.group);
        }
        if (this.config.max_y){
            this.chart.elasticY(false);
            this.chart.y(d3.scaleLinear().domain([0,this.config.max_y]))
        }
        else{
            this.chart.elasticY(true);
        }              
    }

	 _getColorFunction(i){
    	let self=this;
    	return function(d){
    		let val =d.value[self.categories[i]];
    		return val?val:0
    	}
    }

}

export {ColorStackChart};
