import { hierarchy as d3_hierarchy}  from "d3-hierarchy";

const distances ={
	euclidean: (v1, v2) => {
      let total = 0;
      for (let i = 0; i < v1.length; i++) {
         total += (Number.isNaN(v2[i])?0:v2[i] -v1[i]) ** 2;      
      }
      return Math.sqrt(total);
   },
   manhattan: (v1, v2) => {
     let total = 0;
     for (let i = 0; i < v1.length ; i++) {
        total += Math.abs(v2[i] - v1[i]);      
     }
     return total;
   },
   max: (v1, v2) => {
     let max = 0;
     for (let i = 0; i < v1.length; i++) {
        max = Math.max(max , Math.abs(v2[i] - v1[i]));      
     }
     return max;
   }
}

class HierarchicalClustering{
	constructor(distance, linkage, threshold) {
   		this.distance = distances[distance];
   		if (!this.distance){
   			this.distance=distances["euclidean"];
   		}
   		this.linkage = linkage?linkage:"average";
		this.threshold = threshold === undefined ? Number.POSITIVE_INFINITY : threshold;
	}

   cluster(items, snapshotPeriod, snapshotCb) {
      this.clusters = [];
      this.dists = [];  // distances between each pair of clusters
      this.mins = []; // closest cluster for each cluster
      this.index = []; // keep a hash of all clusters by key
 
      
		  for (let i = 0; i < items.length; i++) {
			 const cluster = {
				value: items[i],
				key: i,
				index: i,
				size: 1
			 };
			 this.clusters[i] = cluster;
			 this.index[i] = cluster;
			 this.dists[i] = [];
			 this.mins[i] = 0;
		  }
      

      for (let i = 0; i < this.clusters.length; i++) {
         for (let j = 0; j <= i; j++) {
            const dist = (i === j) ? Number.POSITIVE_INFINITY : 
               this.distance(this.clusters[i].value, this.clusters[j].value);
            this.dists[i][j] = dist;
            this.dists[j][i] = dist;

            if (dist < this.dists[i][this.mins[i]]) {
               this.mins[i] = j;               
            }
         }
      }

      let merged = this.mergeClosest();
      let i = 0;
      while (merged) {
        if (snapshotCb && (i++ % snapshotPeriod) === 0) {
           snapshotCb(this.clusters);           
        }
        merged = this.mergeClosest();
      }
    
      this.clusters.forEach((cluster) => {
        // clean up metadata used for clustering
        cluster.key = undefined;
        cluster.index = undefined;
      });
      this.node_order=[];
      this.order=0;
      this.orderClusters(this.clusters[0]);
      return this.clusters;
   }

   mergeClosest() {
      // find two closest clusters from cached mins
      let minKey = 0;
      let min = Number.POSITIVE_INFINITY;
      for (let i = 0; i < this.clusters.length; i++) {
         const key = this.clusters[i].key;
         const dist = this.dists[key][this.mins[key]];
         if (dist < min) {
            minKey = key;
            min = dist;
         }
      }
      if (min >= this.threshold) {
         return false;         
      }

      const c1 = this.index[minKey];
      const c2 = this.index[this.mins[minKey]];

      // merge two closest clusters
      const merged = {
         children: [c1, c2],
         key: c1.key,
         size: c1.size + c2.size
      };

      this.clusters[c1.index] = merged;
      this.clusters.splice(c2.index, 1);
      this.index[c1.key] = merged;

      // update distances with new merged cluster
      for (let i = 0; i < this.clusters.length; i++) {
         const ci = this.clusters[i];
         let dist;
         if (c1.key === ci.key) {
            dist = Number.POSITIVE_INFINITY;            
         }
         else if (this.linkage === "single") {
            dist = this.dists[c1.key][ci.key];
            if (this.dists[c1.key][ci.key] > this.dists[c2.key][ci.key]) {
               dist = this.dists[c2.key][ci.key];
            }
         }
         else if (this.linkage === "complete") {
            dist = this.dists[c1.key][ci.key];
            if (this.dists[c1.key][ci.key] < this.dists[c2.key][ci.key]) {
               dist = this.dists[c2.key][ci.key];              
            }
         }
         else if (this.linkage === "average") {
            dist = (this.dists[c1.key][ci.key] * c1.size
                   + this.dists[c2.key][ci.key] * c2.size) / (c1.size + c2.size);
         }
         else {
            dist = this.distance(ci.value, c1.value);            
         }

         this.dists[c1.key][ci.key] = this.dists[ci.key][c1.key] = dist;
      }

    
      // update cached mins
      for (let i = 0; i < this.clusters.length; i++) {
         const key1 = this.clusters[i].key;        
         if (this.mins[key1] === c1.key || this.mins[key1] === c2.key) {
            let min = key1;
            for (let j = 0; j < this.clusters.length; j++) {
               const key2 = this.clusters[j].key;
               if (this.dists[key1][key2] < this.dists[key1][min]) {
                  min = key2;                  
               }
            }
            this.mins[key1] = min;
         }
         this.clusters[i].index = i;
      }
      // clean up metadata used for clustering
      c1.key = undefined; c2.key = undefined;
      c1.index = undefined; c2.index = undefined;

      return true;
   }

   orderClusters(obj){
   	    for (const item in obj){
   	    	if (item ==="children"){
   	    		this.orderClusters(obj[item][0]);
                this.orderClusters(obj[item][1])
   	    	}
   	    	if (item==="value"){
                obj[item]._order=this.order++;
   	    		this.node_order.push(obj[item]._id);		
   	    	}
   	    }
   }
}

function getHierarchicalNodes(data,config={}){
    const hc = new HierarchicalClustering(config.distance,config.linkage,config.threshold);
    hc.cluster(data);
    try {
       const nodes = d3_hierarchy(hc.clusters[0]);
       for (const node of nodes.descendants()){
           const d = node.data;
           d.size = undefined;
           if (d.value){
               d.id= d.value._id;
               d.order = d.value._order;
               d.value = undefined;
           }
       }
       return {nodes:nodes,order:hc.node_order}
    } catch (e){
        console.error(e);
        const nodes = {children:[], data: {children:[]}, depth: 0, height: 0, parent: null};
        return { nodes, order:[] };
    }
}


function parseNewick(tree){
    let r={};
    for(let e=[],s=tree.split(/\s*(;|\(|\)|,|:)\s*/),t=0;t<s.length;t++){
        const n=s[t];
        switch(n){
            case"(": {const c={};r.children=[c],e.push(r),r=c;break;
            }
            case",": {const c={};e[e.length-1].children.push(c),r=c;break;
            }
            case")":r=e.pop();break;case":":break;
            default: {const h=s[t-1];")"===h||"("===h||","===h?r.name=n:":"===h&&(r.length=Number.parseFloat(n));
            }
        }
    }
    return d3_hierarchy(r);    
}

export {HierarchicalClustering,getHierarchicalNodes,parseNewick};

