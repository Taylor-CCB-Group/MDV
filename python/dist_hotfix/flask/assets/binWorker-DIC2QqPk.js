(function(){"use strict";const l=t=>{const f=t.data[2][1]==="int32"?Int32Array:Float32Array,r=new f(t.data[2][0]),n=t.data[3],d=(n.max-n.min)/n.bins,o=n.max,i=n.min,g=r.length,e=new Array(n.bins+1).fill(0),m=new Uint8Array(t.data[0]),c=new Uint8Array(t.data[1]);for(let a=0;a<g;a++){if(c[a]!==0&&c[a]!==m[a])continue;let s=r[a];s=s>o?o:s<i?i:s,e[Math.floor((s-i)/d)]++}return e};self.onmessage=t=>{t.data.length===void 0||typeof t.data=="string"||self.postMessage(l(t))}})();
//# sourceMappingURL=binWorker-DIC2QqPk.js.map