
class DataModel {
    constructor(dataStore,config={}){
        this.dataStore = dataStore;
        this.data = new Int32Array(0);
        let len  = 0;
        for (let n=0;n<len;n++){
            this.data[n]=n;
        }
        this.size=len;
        if (config.autoupdate == null || config.autoupdate===true){
            dataStore.addListener("sss",(type)=>{
                if (type==="filtered"){
                    this.updateModel();
                }
            });
        }
       
        this.listeners={}
    }

    getLength(){
        return this.size;
    }

    addListener(id,func){
        this.listeners[id]=func;
    }

    sort(column,order){
        const col = this.dataStore.getRawColumn(column);
        if (order=== "desc"){
            this.data.sort((a,b)=>{
                let va= col[a];
                let vb = col[b];
                va= isNaN(va)?Number.MAX_VALUE:va;
                vb= isNaN(vb)?Number.MAX_VALUE:vb;
                return va-vb;
            });
        }
        else{
            this.data.sort((a,b)=>{
                let va= col[a];
                let vb = col[b];
                va= isNaN(va)?Number.MIN_VALUE:va;
                vb= isNaN(vb)?Number.MIN_VALUE:vb;
                return vb-va;
            });
        }
      


    }

    updateModel(){
        const arr = this.dataStore.filterArray;
        this.data =new Int32Array(this.dataStore.filterSize);
        const data = this.data;
        const len = this.dataStore.size;
        let index=0;
        for (let n=0;n<len;n++){
            if(arr[n]===0){
                data[index++]=n
            }
        }
        this.size=this.dataStore.filterSize;
        for (let id in this.listeners){
            this.listeners[id]()
        }
    }

    getItem(index){       
        return  this.dataStore.getRowAsObject(this.data[index]);
    }


}

export {DataModel}