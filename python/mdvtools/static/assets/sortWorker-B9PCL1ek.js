(function(){"use strict";const f={int32:{arr:Int32Array,type:"number"},double:{arr:Float32Array,type:"number"},integer:{arr:Int32Array,type:"number"},text16:{arr:Uint16Array,type:"text"},text:{arr:Uint8Array,type:"text"},unique:{arr:Uint8Array,type:"unique"},multitext:{arr:Uint8Array,type:"multitext"}};onmessage=function(n){const{orderBuffer:c,columns:u}=n.data,r=new Uint32Array(c),s=u.map(t=>{const e=f[t.datatype],o=new e.arr(t.buffer);let i=null;return e.type==="text"?i=p({data:o,desc:t.desc,values:t.values}):e.type==="unique"?i=y({data:o,desc:t.desc,size:t.stringLength}):e.type==="multitext"?i=g({data:o,desc:t.desc,size:t.stringLength}):i=d({data:o,desc:t.desc}),i});r.sort((t,e)=>{for(let o of s){const i=o(t,e);if(i!==0)return i}return 0}),this.postMessage("done")};function d({data:n,desc:c}){const u=c?-1:1;return(r,s)=>{let t=n[r],e=n[s];return t=isNaN(t)?Number.MAX_VALUE:t,e=isNaN(e)?Number.MAX_VALUE:e,(t>e?1:t<e?-1:0)*u}}function p({data:n,values:c,desc:u}){const r=l(c,u);return(s,t)=>r[n[s]]-r[n[t]]}function y({data:n,desc:c,size:u}){const r={},s=new TextDecoder,t=order===c?-1:1;for(let e=0;e<n.length;e++){const o=n[e];r[o]=s.decode(col.slice(o*u,o*u+u))}return(e,o)=>r[e].localeCompare(r[o])*t}function g({data:n,desc:c,values:u,size:r}){const s=l(u,c);return s[65536]=65356,(t,e)=>{const o=n.slice(t*r,t*r+r),i=n.slice(e*r,e*r+r);for(let a=0;a<r;a++){const m=s[o[a]]-s[i[a]];if(m!==0)return m;if(o[a]===65536&&i[a]===65536)return 0}return 0}}function l(n,c){const u=c?-1:1,r=n.map((t,e)=>[t,e]).sort((t,e)=>t[0].localeCompare(e[0])*u),s={};for(let t=0;t<r.length;t++)s[r[t][1]]=t;return s}})();
//# sourceMappingURL=sortWorker-B9PCL1ek.js.map
