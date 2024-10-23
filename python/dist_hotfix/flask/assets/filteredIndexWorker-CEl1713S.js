(function(){"use strict";onmessage=e=>{const r=new Uint8Array(e.data.filterBuffer),n=new Uint32Array(e.data.outputBuffer);for(let t=0,a=0;t<r.length;t++)r[t]===0&&(n[a++]=t);postMessage({type:"done"})}})();
//# sourceMappingURL=filteredIndexWorker-CEl1713S.js.map
