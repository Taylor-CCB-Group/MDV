class Dimension{
    constructor(column,parent){
        if (typeof column === "string"){
            this.column=[column]
        }
        else{
            this.column=column;
        }

        this.filterBuffer = new SharedArrayBuffer(parent.size);       
        this.parent=parent;
        this.filterArray= new Uint8Array(this.filterBuffer);
    }
    removeFilter(){     
        const filter = this.parent.filterArray;
        const localFilter= this.filterArray;
        const parent = this.parent;
        for (let i=0;i<this.parent.size;i++){
            if (localFilter[i]===1){
                if(--filter[i]==0){
                    parent.filterSize++;
                }
                localFilter[i]=0;
                
            }
        }
        this.parent._callListeners("filtered",this); 
    }

    getLocalFilter(){
        return this.filterArray;
    }

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
        this.parent._callListeners("filtered",this);
    

    }

 

    destroy(){    
        this.removeFilter();
        this.filterBuffer=null;
        this.localFilter=null;
        this.parent.dimensions.splice(this.parent.dimensions.indexOf(this),1)
    }


}

Dimension.types={};

export default Dimension;