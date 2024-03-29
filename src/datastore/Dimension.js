class Dimension{


    constructor(parent){ 
        this.filterBuffer = new SharedArrayBuffer(parent.size);       
        this.parent=parent;
        this.filterArray= new Uint8Array(this.filterBuffer);
        this.filterArguments =null;
        this.filterIndexes = null;
        this.filterMethod = null;
    }

    removeFilter(notify=true){  
        if (!this.filterMethod){
            return;
        }  
        const filter = this.parent.filterArray;
        const localFilter= this.filterArray;
        const parent = this.parent;
        const len = this.parent.size;
        delete this.filterArguments;
        delete this.filterColumns;
        delete this.filterMethod;
        for (let i=0;i<len;i++){
            if (localFilter[i]===1){
                if(--filter[i]==0){
                    parent.filterSize++;
                }
                localFilter[i]=0;           
            }
        }
        if (this.bgfArray){     
            for (let i=0;i<len;i++){          
                if (this.bgfArray[i]===0){
                    //need to filter in         
                    if(--filter[i]===0){
                        parent.filterSize++;
                    } 
                    localFilter[i]=2;                          
                }
            }
        }
        if(notify){
            this.parent._callListeners("filtered",this); 
        }    
    }

    getLocalFilter(){
        return this.filterArray;
    }

    //needs to be removed - produces strange behaviour
    filterOnIndex(indexSet){
        const filter = this.parent.filterArray;
        const len = this.parent.size;
        const localFilter= this.filterArray;
        const parent = this.parent;
        for (let i=0;i<len;i++){
          
            if (indexSet.has(i)){
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
    
    /**
    * Filters  the data
    * @param {string} method- The name of the filter method
    * @param {string[]} columns - a list of column ids used in the filtering
    * @param {object} args - any extra arguments for filtering
    * @param {boolean} [notify=true] - notify any listeners in the dataStore that the 
    * data has been filtered
    */
    filter(method,columns,args,notify=true){
        const t = performance.now();
        this.filterColumns =columns;
        this.filterArguments = args;
        this.filterMethod = method;
        this[method](args,columns);
        if (this.bgfArray){
            const filter = this.parent.filterArray;
            const len = this.parent.size;
            const localFilter= this.filterArray;
            const parent = this.parent;
            for (let i=0;i<len;i++){          
                if (this.bgfArray[i]===0){
                  
                        if(++filter[i]===1){
                            parent.filterSize--;
                        }          
                    localFilter[i]=2;
                }
            }
        }
        if (notify){
            this.parent._callListeners("filtered",this);
        }
        console.log(`method${method}: ${performance.now()-t}`);
    }


    /**
    * updates the filter if data in the supplied columns has been altered
    * does not propagate (call any listeners)
    * @param {string[]} columns - a list of columns whose data has changes
    * @returns {boolean}  True if refiltering has taken place 
    */
    reFilterOnDataChanged(columns){
        if (this.filterMethod && this.filterColumns){
            for (let c of this.filterColumns){
                if (typeof c === "string" && columns.indexOf(c)!==-1){
                    this.filter(this.filterMethod,this.filterColumns,this.filterArguments,false);
                    return true;
                }
            }        
        }
        return false;
    }

    /**
    * Experimental - should be invoked if the datastore has increased in size
    * Will increase size of local filters 
    */
    updateSize(){
        let newBuff =  new SharedArrayBuffer(this.parent.size);
        let newArr = new Uint8Array(newBuff);
        newArr.set(this.filterArray);
        this.filterBuffer =newBuff;  
        this.filterArray= newArr;
        if (this.bgfArray){
            const bd= this.bgfData;
            this.setBackgroundFilter(b.column.b.cat);
        }
        if (this.filterMethod){          
            this.filter(this.filterMethod,this.filterColumns,this.filterArguments,false);
        }
    }


    /**
    * sets a permenant filter on the chart
    * @param {string} column- The column of the filter
    * @param {string}  cat - the category in the column to filter on
    */
    setBackgroundFilter(column,cat){  
        const col = this.parent.columnIndex[column];
        const ci =   col.values.indexOf(cat);
        const data = col.data;
        this.bgfData={
            column:column,
            cat:cat
        };
        this.bgfArray= new Uint8Array(this.parent.size);
        for (let i=0;i<this.parent.size;i++){
            if (data[i]===ci){
                this.bgfArray[i]=1;
                
            }
            else{
                this.filterArray[i]=2;
            }
        }
    }

    clearBackGroundFilter(){
        this.bgfArray=null;
        this.bgfData=null;
    }

    destroy(notify=true){    
        this.removeFilter(notify);
        this.filterBuffer=null;
        this.localFilter=null;
        this.parent.dimensions.splice(this.parent.dimensions.indexOf(this),1);
    }
}

Dimension.types={base_dimension:Dimension};
export default Dimension;