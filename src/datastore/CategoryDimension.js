import Dimension from "./Dimension.js";
class CategoryDimension extends Dimension{
    constructor(column,parent){
        super(column, parent)
        this.worker= new Worker (new URL("./catWorker.js",import.meta.url));
      
    }

    filterCategories(category,index=0){
        const parent = this.parent;
        const col = parent.columnIndex[this.column[index]];
        const data = col.data;
        const filter = parent.filterArray;
        const localFilter= this.filterArray;
        const t = performance.now();
     
        const vals = col.values;
        let ind = 0;
        const cats=new Set();
        const len = this.parent.size;
        
        if (typeof category === "string"){
            const ind = vals.indexOf(category);
            for (let i=0;i<len;i++){
            
                if (data[i]===ind){
                    if (localFilter[i]===1){
                       if (--filter[i] === 0){
                            parent.filterSize++;
                        }
                    }
                    localFilter[i]=0;              
                }
                else{                  
                    if (localFilter[i]===0){
                        if(++filter[i]===1){
                            parent.filterSize--;
                        }
                    }              
                    localFilter[i]=1
                }
            }
           
        }
        else{
            for (let cat of category){
                cats.add(vals.indexOf(cat));
            }
            for (let i=0;i<len;i++){
          
                if (cats.has(data[i])){
                    if (localFilter[i]===1){
                        if (--filter[i] === 0){
                            parent.filterSize++;
                        }
                    }
                    localFilter[i]=0;              
                }
                else{                  
                    if (localFilter[i]===0){
                        if(++filter[i]===1){
                            parent.filterSize--;
                        }
                    }              
                    localFilter[i]=1
                }
               
            }
        }
        console.log(`cat filtered ${category} : ${performance.now()-t}`);

        this.parent._callListeners("filtered",this);
    }

    getCategories(callback,config={},index=0){
        const col = this.parent.columnIndex[this.column[index]];
        config.values=col.values;
        const t = performance.now();
           
    

        const action = (e)=>{
            console.log(`calc categories ${col.name} : ${performance.now()-t}`);
            callback(e.data);
            this.worker.removeEventListener("message",action)    
        }
        
        this.worker.addEventListener("message",action);
        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            col.buffer,
            config  
         ]);
           

    }

    destroy(){
        super.destroy();
        this.worker.terminate();
    }
}


Dimension.types["category_dimension"]=CategoryDimension;

export default CategoryDimension; 
