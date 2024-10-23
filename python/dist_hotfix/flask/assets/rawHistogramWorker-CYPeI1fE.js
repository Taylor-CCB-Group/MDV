(function(){"use strict";self.onmessage=async e=>{const{isInt32:i,data:c,min:a,max:l,bins:t}=e.data,f=i?Int32Array:Float32Array,o=new f(c),r=new Array(t).fill(0),y=(l-a)/t;for(let n=0;n<o.length;n++){const s=Math.floor((o[n]-a)/y);s>=0&&s<t&&r[s]++}self.postMessage(r)}})();
//# sourceMappingURL=rawHistogramWorker-CYPeI1fE.js.map
