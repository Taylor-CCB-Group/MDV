import { hierarchy as d3_hierarchy}  from "d3-hierarchy";

const distances ={
	euclidean: function(v1, v2) {
      var total = 0;
      for (var i = 0; i < v1.length; i++) {
         total += Math.pow(isNaN(v2[i])?0:v2[i] -v1[i], 2);      
      }
      return Math.sqrt(total);
   },
   manhattan: function(v1, v2) {
     var total = 0;
     for (var i = 0; i < v1.length ; i++) {
        total += Math.abs(v2[i] - v1[i]);      
     }
     return total;
   },
   max: function(v1, v2) {
     var max = 0;
     for (var i = 0; i < v1.length; i++) {
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
		this.threshold = threshold == undefined ? Infinity : threshold;
	}

   cluster(items, snapshotPeriod, snapshotCb) {
      this.clusters = [];
      this.dists = [];  // distances between each pair of clusters
      this.mins = []; // closest cluster for each cluster
      this.index = []; // keep a hash of all clusters by key
 
      
		  for (var i = 0; i < items.length; i++) {
			 var cluster = {
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
      

      for (var i = 0; i < this.clusters.length; i++) {
         for (var j = 0; j <= i; j++) {
            var dist = (i == j) ? Infinity : 
               this.distance(this.clusters[i].value, this.clusters[j].value);
            this.dists[i][j] = dist;
            this.dists[j][i] = dist;

            if (dist < this.dists[i][this.mins[i]]) {
               this.mins[i] = j;               
            }
         }
      }

      var merged = this.mergeClosest();
      var i = 0;
      while (merged) {
        if (snapshotCb && (i++ % snapshotPeriod) == 0) {
           snapshotCb(this.clusters);           
        }
        merged = this.mergeClosest();
      }
    
      this.clusters.forEach(function(cluster) {
        // clean up metadata used for clustering
        delete cluster.key;
        delete cluster.index;
      });
      this.node_order=[];
      this.order=0;
      this.orderClusters(this.clusters[0]);
      return this.clusters;
   }

   mergeClosest() {
      // find two closest clusters from cached mins
      var minKey = 0, min = Infinity;
      for (var i = 0; i < this.clusters.length; i++) {
         var key = this.clusters[i].key,
             dist = this.dists[key][this.mins[key]];
         if (dist < min) {
            minKey = key;
            min = dist;
         }
      }
      if (min >= this.threshold) {
         return false;         
      }

      var c1 = this.index[minKey],
          c2 = this.index[this.mins[minKey]];

      // merge two closest clusters
      var merged = {
         children: [c1, c2],
         key: c1.key,
         size: c1.size + c2.size
      };

      this.clusters[c1.index] = merged;
      this.clusters.splice(c2.index, 1);
      this.index[c1.key] = merged;

      // update distances with new merged cluster
      for (var i = 0; i < this.clusters.length; i++) {
         var ci = this.clusters[i];
         var dist;
         if (c1.key == ci.key) {
            dist = Infinity;            
         }
         else if (this.linkage == "single") {
            dist = this.dists[c1.key][ci.key];
            if (this.dists[c1.key][ci.key] > this.dists[c2.key][ci.key]) {
               dist = this.dists[c2.key][ci.key];
            }
         }
         else if (this.linkage == "complete") {
            dist = this.dists[c1.key][ci.key];
            if (this.dists[c1.key][ci.key] < this.dists[c2.key][ci.key]) {
               dist = this.dists[c2.key][ci.key];              
            }
         }
         else if (this.linkage == "average") {
            dist = (this.dists[c1.key][ci.key] * c1.size
                   + this.dists[c2.key][ci.key] * c2.size) / (c1.size + c2.size);
         }
         else {
            dist = this.distance(ci.value, c1.value);            
         }

         this.dists[c1.key][ci.key] = this.dists[ci.key][c1.key] = dist;
      }

    
      // update cached mins
      for (var i = 0; i < this.clusters.length; i++) {
         var key1 = this.clusters[i].key;        
         if (this.mins[key1] == c1.key || this.mins[key1] == c2.key) {
            var min = key1;
            for (var j = 0; j < this.clusters.length; j++) {
               var key2 = this.clusters[j].key;
               if (this.dists[key1][key2] < this.dists[key1][min]) {
                  min = key2;                  
               }
            }
            this.mins[key1] = min;
         }
         this.clusters[i].index = i;
      }
      // clean up metadata used for clustering
      delete c1.key; delete c2.key;
      delete c1.index; delete c2.index;

      return true;
   }

   orderClusters(obj){
   	    for (let item in obj){
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
    const nodes =d3_hierarchy(hc.clusters[0]);
    for (let node of nodes.descendants()){
        let d = node.data;
        delete d.size;
        if (d.value){
            d.id= d.value._id;
            d.order = d.value._order;
            delete d.value;
        }
    }
    return {nodes:nodes,order:hc.node_order}
}


function parseNewick(tree){
    let r={};
    for(let e=[],s=tree.split(/\s*(;|\(|\)|,|:)\s*/),t=0;t<s.length;t++){
        let n=s[t];
        switch(n){
            case"(":var c={};r.children=[c],e.push(r),r=c;break;
            case",":var c={};e[e.length-1].children.push(c),r=c;break;
            case")":r=e.pop();break;case":":break;
            default:var h=s[t-1];")"===h||"("===h||","===h?r.name=n:":"===h&&(r.length=parseFloat(n));
        }
    }
    return d3_hierarchy(r);    
}

export {HierarchicalClustering,getHierarchicalNodes,parseNewick};

