(function(){"use strict";onmessage=function(e){const r=new Uint8Array(e.data.filterBuffer),n=new Uint32Array(e.data.outputBuffer);for(let t=0,a=0;t<r.length;t++)r[t]===0&&(n[a++]=t);postMessage({type:"done"})}})();
//# sourceMappingURL=filteredIndexWorker-DHIC0XDY.js.map
