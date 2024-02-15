import Dimension from "./Dimension.js";
class CategoryDimension extends Dimension{
    constructor(parent){
        super(parent)
        this.worker= new Worker (new URL("./catWorker.js?v=11",import.meta.url));    
    }

    filterCategoriesPredicate(args,columns){
        const category =args;
        const col = this.parent.columnIndex[columns[0]];
        const data = col.data;
        const vals = col.values;
        
        if (typeof category === "string"){
            const ind = vals.indexOf(category);
            
            if (col.datatype==="multitext"){
                const int = col.stringLength;
                const predicate = i => {
                    const st = i*int;
                    for (let n=st; n < st+int; n++) {
                        if (data[n] === ind) return true;
                    }
                    return false;
                }
                return this.filterPredicate({predicate});
            } //end multitext
            else {
                const predicate = i => data[i] === ind;
                return this.filterPredicate({predicate});
            }     
        } //end typeof category === "string"
        else { //category is an array
            const cats = new Set();
            for (let cat of category){
                cats.add(vals.indexOf(cat));
            }
            if (col.datatype==="multitext"){
                const int = col.stringLength;
                let ao = category.operand === "and";
                const catsArr = category.map(c => vals.indexOf(c));
                const singleCat = cats.size === 1;

                const predicate = i => {
                    const st = i*int;
                    //nb, this is based on the premise that we want "and" of ['a', 'a'] to require two 'a's, so using array rather than Set
                    //in any case, we need to remove found categories for a given row as we go,
                    //otherwise we end up with a false positive if the same category is used twice
                    let catsToFind = ao ? catsArr.slice(0) : cats; //nb slice is faster than [...spread]
                    for (let n=st; n < st+int; n++) {
                        if (data[n]===65535) {
                            //... if this row only has one item, we can't possibly have two of the same item,
                            //but I think we still want to match if we're looking for a single category...
                            //e.g. if we have 'Count Periostin' and we click the 'Periostin / Periostin' element of an interaction matrix.
                            //In that case, `catsToFind` will have fewer items than `catsArr`
                            //(testing with cellpairs where I want but should be right for arbitrary multitext??)
                            if (ao && singleCat && catsToFind.length < catsArr.length) {
                                return true;
                            }
                        }
                        if (ao) {
                            const index = catsToFind.indexOf(data[n]);
                            if (index !== -1) {
                                catsToFind.splice(index, 1);
                                if (catsToFind.length === 0) {
                                    return true;
                                }
                            }
                        } else {
                            if (catsToFind.has(data[n])){
                                return true;
                            }
                        }
                    }
                    return false;
                }
                return this.filterPredicate({predicate});
            } //end multitext
            else{
                //this appears to be taking about twice as long as non-predicate version??
                const catsArr = Array.from(cats);
                const predicate = i => cats.has(data[i]);
                // const predicate = i => catsArr.includes(data[i]);
                return this.filterPredicate({predicate});
            }
        }  
    }
    filterCategories(args,columns){
        const category =args;
        const parent = this.parent;
        const col = parent.columnIndex[columns[0]];
        const data = col.data;
        const filter = parent.filterArray;
        const localFilter= this.filterArray;
      
     
        const vals = col.values;
        const cats = new Set();
        const len = this.parent.size;
        
        if (typeof category === "string"){
            const ind = vals.indexOf(category);
           
            if (col.datatype==="multitext"){
                const int = col.stringLength;
                for (let i=0;i<len;i++){
                    const st = i*int;
                    let has =false;
                    for (let n=st;n<st+int;n++){
                        if (data[n]===ind){
                            has=true;
                            break
                        }
                    }
                
                    if (has){
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
        }
        else {
            for (let cat of category){
                cats.add(vals.indexOf(cat));
            }
            if (col.datatype==="multitext"){
                const int = col.stringLength;
                let ao = category.operand === "and";
                const catsArr = category.map(c => vals.indexOf(c));
                const singleCat = cats.size === 1;
                for (let i=0;i<len;i++){
                    const st = i*int;
                    let has = false;
                    //nb, we actually want "and" of ['a', 'a'] to require two 'a's, so using array rather than Set
                    //in any case, we need to remove found categories for a given row as we go,
                    //otherwise we end up with a false positive if the same category is used twice
                    
                    const catsToFind = ao ? catsArr.slice(0) : cats; //nb slice is faster than [...spread]
                    
                    for (let n=st; n < st+int; n++) {
                        if (data[n]===65535) {
                            //... if this row only has one item, we can't possibly have two of the same item,
                            //but I think we still want to match if we're looking for a single category...
                            //e.g. if we have 'Count Periostin' and we click the 'Periostin / Periostin' element of an interaction matrix.
                            //In that case, `catsToFind` will have fewer items than `catsArr`
                            //(testing with cellpairs where I want but should be right for arbitrary multitext??)
                            if (ao && singleCat && catsToFind.length < catsArr.length) has = true;
                            break;
                        }
                        if (ao) {
                            const index = catsToFind.indexOf(data[n]);
                            if (index !== -1) {
                                catsToFind.splice(index, 1);
                                if (catsToFind.length === 0) {
                                    has = true;
                                    break;
                                }
                            }
                        } else {
                            if (catsToFind.has(data[n])){
                                has=true;
                                break
                            }
                        }
                    }
                    if (has) {
                        if (localFilter[i]===1){
                            if (--filter[i] === 0){
                                    parent.filterSize++;
                                }
                            }
                            localFilter[i]=0;              
                    } else {                  
                        if (localFilter[i]===0){
                            if(++filter[i]===1){
                                parent.filterSize--;
                            }
                        }              
                        localFilter[i]=1
                    }
                }
            } //end multitext
            else{
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
        }  
    }

    filterMultipleCategories(args,columns){
        const categories =args;
        const parent = this.parent;
        const data= [];
        const indexes=[];
        for (let i=0;i<columns.length;i++){
            const col = parent.columnIndex[columns[i]];
            data.push(col.data);
            indexes.push(col.values.indexOf(args[i]))
        }
        const filter = parent.filterArray;
        const localFilter= this.filterArray;
        const len =parent.size;
        for (let i=0;i<len;i++){
            let con = true;
            for (let n=0;n<data.length;n++){
                if (data[n][i]!==indexes[n]){
                    con=false;
                    break;
                }
            }
            if (con){
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

    getSankeyData(callback,columns,config={}){
        const col1 = this.parent.columnIndex[columns[0]];
        const col2 = this.parent.columnIndex[columns[1]];
        config.values=col1.values;
        config.values2 = col2.values;

        config.method=config.method || "sankey";
        if (config.method==="proportion"){
            config.cats = config.category == null?null:this.parent.columnIndex[columns[2]].buffer;
        }
        const t = performance.now();
           
        const action = (e)=>{
            console.log(`calc sankey ${col1.name} ${col2.name} : ${performance.now()-t}`);
            callback(e.data);
            this.worker.removeEventListener("message",action)    
        }
        
        this.worker.addEventListener("message",action);
        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            col1.buffer,
            config ,
            col2.buffer
         ]);
    }



    getCategories(callback,column,config={}){
        const col = this.parent.columnIndex[column];
        config.values=col.values;
        config.datatype=col.datatype;
        config.stringLength= col.stringLength;
        const t = performance.now();
        const action = (e)=>{
            //console.log(`calc categories ${col.name} : ${performance.now()-t}`);
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
