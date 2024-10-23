(function(){"use strict";var we=typeof globalThis<"u"?globalThis:typeof window<"u"?window:typeof global<"u"?global:typeof self<"u"?self:{},me={exports:{}};/* @license
Papa Parse
v5.4.1
https://github.com/mholt/PapaParse
License: MIT
*/(function(ce,pe){(function(he,_){ce.exports=_()})(we,function he(){var _=typeof self<"u"?self:typeof window<"u"?window:_!==void 0?_:{},Y=!_.document&&!!_.postMessage,le=_.IS_PAPA_WORKER||!1,M={},fe=0,o={parse:function(t,e){var r=(e=e||{}).dynamicTyping||!1;if(p(r)&&(e.dynamicTypingFunction=r,r={}),e.dynamicTyping=r,e.transform=!!p(e.transform)&&e.transform,e.worker&&o.WORKERS_SUPPORTED){var n=function(){if(!o.WORKERS_SUPPORTED)return!1;var f=(A=_.URL||_.webkitURL||null,w=he.toString(),o.BLOB_URL||(o.BLOB_URL=A.createObjectURL(new Blob(["var global = (function() { if (typeof self !== 'undefined') { return self; } if (typeof window !== 'undefined') { return window; } if (typeof global !== 'undefined') { return global; } return {}; })(); global.IS_PAPA_WORKER=true; ","(",w,")();"],{type:"text/javascript"})))),l=new _.Worker(f),A,w;return l.onmessage=H,l.id=fe++,M[l.id]=l}();return n.userStep=e.step,n.userChunk=e.chunk,n.userComplete=e.complete,n.userError=e.error,e.step=p(e.step),e.chunk=p(e.chunk),e.complete=p(e.complete),e.error=p(e.error),delete e.worker,void n.postMessage({input:t,config:e,workerId:n.id})}var s=null;return o.NODE_STREAM_INPUT,typeof t=="string"?(t=function(f){return f.charCodeAt(0)===65279?f.slice(1):f}(t),s=e.download?new ie(e):new X(e)):t.readable===!0&&p(t.read)&&p(t.on)?s=new ae(e):(_.File&&t instanceof File||t instanceof Object)&&(s=new se(e)),s.stream(t)},unparse:function(t,e){var r=!1,n=!0,s=",",f=`\r
`,l='"',A=l+l,w=!1,a=null,E=!1;(function(){if(typeof e=="object"){if(typeof e.delimiter!="string"||o.BAD_DELIMITERS.filter(function(i){return e.delimiter.indexOf(i)!==-1}).length||(s=e.delimiter),(typeof e.quotes=="boolean"||typeof e.quotes=="function"||Array.isArray(e.quotes))&&(r=e.quotes),typeof e.skipEmptyLines!="boolean"&&typeof e.skipEmptyLines!="string"||(w=e.skipEmptyLines),typeof e.newline=="string"&&(f=e.newline),typeof e.quoteChar=="string"&&(l=e.quoteChar),typeof e.header=="boolean"&&(n=e.header),Array.isArray(e.columns)){if(e.columns.length===0)throw new Error("Option columns is empty");a=e.columns}e.escapeChar!==void 0&&(A=e.escapeChar+l),(typeof e.escapeFormulae=="boolean"||e.escapeFormulae instanceof RegExp)&&(E=e.escapeFormulae instanceof RegExp?e.escapeFormulae:/^[=+\-@\t\r].*$/)}})();var h=new RegExp(j(l),"g");if(typeof t=="string"&&(t=JSON.parse(t)),Array.isArray(t)){if(!t.length||Array.isArray(t[0]))return q(null,t,w);if(typeof t[0]=="object")return q(a||Object.keys(t[0]),t,w)}else if(typeof t=="object")return typeof t.data=="string"&&(t.data=JSON.parse(t.data)),Array.isArray(t.data)&&(t.fields||(t.fields=t.meta&&t.meta.fields||a),t.fields||(t.fields=Array.isArray(t.data[0])?t.fields:typeof t.data[0]=="object"?Object.keys(t.data[0]):[]),Array.isArray(t.data[0])||typeof t.data[0]=="object"||(t.data=[t.data])),q(t.fields||[],t.data||[],w);throw new Error("Unable to serialize unrecognized input");function q(i,y,S){var b="";typeof i=="string"&&(i=JSON.parse(i)),typeof y=="string"&&(y=JSON.parse(y));var I=Array.isArray(i)&&0<i.length,x=!Array.isArray(y[0]);if(I&&n){for(var T=0;T<i.length;T++)0<T&&(b+=s),b+=D(i[T],T);0<y.length&&(b+=f)}for(var u=0;u<y.length;u++){var d=I?i.length:y[u].length,v=!1,O=I?Object.keys(y[u]).length===0:y[u].length===0;if(S&&!I&&(v=S==="greedy"?y[u].join("").trim()==="":y[u].length===1&&y[u][0].length===0),S==="greedy"&&I){for(var g=[],L=0;L<d;L++){var C=x?i[L]:L;g.push(y[u][C])}v=g.join("").trim()===""}if(!v){for(var m=0;m<d;m++){0<m&&!O&&(b+=s);var B=I&&x?i[m]:m;b+=D(y[u][B],m)}u<y.length-1&&(!S||0<d&&!O)&&(b+=f)}}return b}function D(i,y){if(i==null)return"";if(i.constructor===Date)return JSON.stringify(i).slice(1,25);var S=!1;E&&typeof i=="string"&&E.test(i)&&(i="'"+i,S=!0);var b=i.toString().replace(h,A);return(S=S||r===!0||typeof r=="function"&&r(i,y)||Array.isArray(r)&&r[y]||function(I,x){for(var T=0;T<x.length;T++)if(-1<I.indexOf(x[T]))return!0;return!1}(b,o.BAD_DELIMITERS)||-1<b.indexOf(s)||b.charAt(0)===" "||b.charAt(b.length-1)===" ")?l+b+l:b}}};if(o.RECORD_SEP="",o.UNIT_SEP="",o.BYTE_ORDER_MARK="\uFEFF",o.BAD_DELIMITERS=["\r",`
`,'"',o.BYTE_ORDER_MARK],o.WORKERS_SUPPORTED=!Y&&!!_.Worker,o.NODE_STREAM_INPUT=1,o.LocalChunkSize=10485760,o.RemoteChunkSize=5242880,o.DefaultDelimiter=",",o.Parser=ee,o.ParserHandle=oe,o.NetworkStreamer=ie,o.FileStreamer=se,o.StringStreamer=X,o.ReadableStreamStreamer=ae,_.jQuery){var K=_.jQuery;K.fn.parse=function(t){var e=t.config||{},r=[];return this.each(function(f){if(!(K(this).prop("tagName").toUpperCase()==="INPUT"&&K(this).attr("type").toLowerCase()==="file"&&_.FileReader)||!this.files||this.files.length===0)return!0;for(var l=0;l<this.files.length;l++)r.push({file:this.files[l],inputElem:this,instanceConfig:K.extend({},e)})}),n(),this;function n(){if(r.length!==0){var f,l,A,w,a=r[0];if(p(t.before)){var E=t.before(a.file,a.inputElem);if(typeof E=="object"){if(E.action==="abort")return f="AbortError",l=a.file,A=a.inputElem,w=E.reason,void(p(t.error)&&t.error({name:f},l,A,w));if(E.action==="skip")return void s();typeof E.config=="object"&&(a.instanceConfig=K.extend(a.instanceConfig,E.config))}else if(E==="skip")return void s()}var h=a.instanceConfig.complete;a.instanceConfig.complete=function(q){p(h)&&h(q,a.file,a.inputElem),s()},o.parse(a.file,a.instanceConfig)}else p(t.complete)&&t.complete()}function s(){r.splice(0,1),n()}}}function F(t){this._handle=null,this._finished=!1,this._completed=!1,this._halted=!1,this._input=null,this._baseIndex=0,this._partialLine="",this._rowCount=0,this._start=0,this._nextChunk=null,this.isFirstChunk=!0,this._completeResults={data:[],errors:[],meta:{}},(function(e){var r=ge(e);r.chunkSize=parseInt(r.chunkSize),e.step||e.chunk||(r.chunkSize=null),this._handle=new oe(r),(this._handle.streamer=this)._config=r}).call(this,t),this.parseChunk=function(e,r){if(this.isFirstChunk&&p(this._config.beforeFirstChunk)){var n=this._config.beforeFirstChunk(e);n!==void 0&&(e=n)}this.isFirstChunk=!1,this._halted=!1;var s=this._partialLine+e;this._partialLine="";var f=this._handle.parse(s,this._baseIndex,!this._finished);if(!this._handle.paused()&&!this._handle.aborted()){var l=f.meta.cursor;this._finished||(this._partialLine=s.substring(l-this._baseIndex),this._baseIndex=l),f&&f.data&&(this._rowCount+=f.data.length);var A=this._finished||this._config.preview&&this._rowCount>=this._config.preview;if(le)_.postMessage({results:f,workerId:o.WORKER_ID,finished:A});else if(p(this._config.chunk)&&!r){if(this._config.chunk(f,this._handle),this._handle.paused()||this._handle.aborted())return void(this._halted=!0);f=void 0,this._completeResults=void 0}return this._config.step||this._config.chunk||(this._completeResults.data=this._completeResults.data.concat(f.data),this._completeResults.errors=this._completeResults.errors.concat(f.errors),this._completeResults.meta=f.meta),this._completed||!A||!p(this._config.complete)||f&&f.meta.aborted||(this._config.complete(this._completeResults,this._input),this._completed=!0),A||f&&f.meta.paused||this._nextChunk(),f}this._halted=!0},this._sendError=function(e){p(this._config.error)?this._config.error(e):le&&this._config.error&&_.postMessage({workerId:o.WORKER_ID,error:e,finished:!1})}}function ie(t){var e;(t=t||{}).chunkSize||(t.chunkSize=o.RemoteChunkSize),F.call(this,t),this._nextChunk=Y?function(){this._readChunk(),this._chunkLoaded()}:function(){this._readChunk()},this.stream=function(r){this._input=r,this._nextChunk()},this._readChunk=function(){if(this._finished)this._chunkLoaded();else{if(e=new XMLHttpRequest,this._config.withCredentials&&(e.withCredentials=this._config.withCredentials),Y||(e.onload=V(this._chunkLoaded,this),e.onerror=V(this._chunkError,this)),e.open(this._config.downloadRequestBody?"POST":"GET",this._input,!Y),this._config.downloadRequestHeaders){var r=this._config.downloadRequestHeaders;for(var n in r)e.setRequestHeader(n,r[n])}if(this._config.chunkSize){var s=this._start+this._config.chunkSize-1;e.setRequestHeader("Range","bytes="+this._start+"-"+s)}try{e.send(this._config.downloadRequestBody)}catch(f){this._chunkError(f.message)}Y&&e.status===0&&this._chunkError()}},this._chunkLoaded=function(){e.readyState===4&&(e.status<200||400<=e.status?this._chunkError():(this._start+=this._config.chunkSize?this._config.chunkSize:e.responseText.length,this._finished=!this._config.chunkSize||this._start>=function(r){var n=r.getResponseHeader("Content-Range");return n===null?-1:parseInt(n.substring(n.lastIndexOf("/")+1))}(e),this.parseChunk(e.responseText)))},this._chunkError=function(r){var n=e.statusText||r;this._sendError(new Error(n))}}function se(t){var e,r;(t=t||{}).chunkSize||(t.chunkSize=o.LocalChunkSize),F.call(this,t);var n=typeof FileReader<"u";this.stream=function(s){this._input=s,r=s.slice||s.webkitSlice||s.mozSlice,n?((e=new FileReader).onload=V(this._chunkLoaded,this),e.onerror=V(this._chunkError,this)):e=new FileReaderSync,this._nextChunk()},this._nextChunk=function(){this._finished||this._config.preview&&!(this._rowCount<this._config.preview)||this._readChunk()},this._readChunk=function(){var s=this._input;if(this._config.chunkSize){var f=Math.min(this._start+this._config.chunkSize,this._input.size);s=r.call(s,this._start,f)}var l=e.readAsText(s,this._config.encoding);n||this._chunkLoaded({target:{result:l}})},this._chunkLoaded=function(s){this._start+=this._config.chunkSize,this._finished=!this._config.chunkSize||this._start>=this._input.size,this.parseChunk(s.target.result)},this._chunkError=function(){this._sendError(e.error)}}function X(t){var e;F.call(this,t=t||{}),this.stream=function(r){return e=r,this._nextChunk()},this._nextChunk=function(){if(!this._finished){var r,n=this._config.chunkSize;return n?(r=e.substring(0,n),e=e.substring(n)):(r=e,e=""),this._finished=!e,this.parseChunk(r)}}}function ae(t){F.call(this,t=t||{});var e=[],r=!0,n=!1;this.pause=function(){F.prototype.pause.apply(this,arguments),this._input.pause()},this.resume=function(){F.prototype.resume.apply(this,arguments),this._input.resume()},this.stream=function(s){this._input=s,this._input.on("data",this._streamData),this._input.on("end",this._streamEnd),this._input.on("error",this._streamError)},this._checkIsFinished=function(){n&&e.length===1&&(this._finished=!0)},this._nextChunk=function(){this._checkIsFinished(),e.length?this.parseChunk(e.shift()):r=!0},this._streamData=V(function(s){try{e.push(typeof s=="string"?s:s.toString(this._config.encoding)),r&&(r=!1,this._checkIsFinished(),this.parseChunk(e.shift()))}catch(f){this._streamError(f)}},this),this._streamError=V(function(s){this._streamCleanUp(),this._sendError(s)},this),this._streamEnd=V(function(){this._streamCleanUp(),n=!0,this._streamData("")},this),this._streamCleanUp=V(function(){this._input.removeListener("data",this._streamData),this._input.removeListener("end",this._streamEnd),this._input.removeListener("error",this._streamError)},this)}function oe(t){var e,r,n,s=Math.pow(2,53),f=-s,l=/^\s*-?(\d+\.?|\.\d+|\d+\.\d+)([eE][-+]?\d+)?\s*$/,A=/^((\d{4}-[01]\d-[0-3]\dT[0-2]\d:[0-5]\d:[0-5]\d\.\d+([+-][0-2]\d:[0-5]\d|Z))|(\d{4}-[01]\d-[0-3]\dT[0-2]\d:[0-5]\d:[0-5]\d([+-][0-2]\d:[0-5]\d|Z))|(\d{4}-[01]\d-[0-3]\dT[0-2]\d:[0-5]\d([+-][0-2]\d:[0-5]\d|Z)))$/,w=this,a=0,E=0,h=!1,q=!1,D=[],i={data:[],errors:[],meta:{}};if(p(t.step)){var y=t.step;t.step=function(u){if(i=u,I())b();else{if(b(),i.data.length===0)return;a+=u.data.length,t.preview&&a>t.preview?r.abort():(i.data=i.data[0],y(i,w))}}}function S(u){return t.skipEmptyLines==="greedy"?u.join("").trim()==="":u.length===1&&u[0].length===0}function b(){return i&&n&&(T("Delimiter","UndetectableDelimiter","Unable to auto-detect delimiting character; defaulted to '"+o.DefaultDelimiter+"'"),n=!1),t.skipEmptyLines&&(i.data=i.data.filter(function(u){return!S(u)})),I()&&function(){if(!i)return;function u(v,O){p(t.transformHeader)&&(v=t.transformHeader(v,O)),D.push(v)}if(Array.isArray(i.data[0])){for(var d=0;I()&&d<i.data.length;d++)i.data[d].forEach(u);i.data.splice(0,1)}else i.data.forEach(u)}(),function(){if(!i||!t.header&&!t.dynamicTyping&&!t.transform)return i;function u(v,O){var g,L=t.header?{}:[];for(g=0;g<v.length;g++){var C=g,m=v[g];t.header&&(C=g>=D.length?"__parsed_extra":D[g]),t.transform&&(m=t.transform(m,C)),m=x(C,m),C==="__parsed_extra"?(L[C]=L[C]||[],L[C].push(m)):L[C]=m}return t.header&&(g>D.length?T("FieldMismatch","TooManyFields","Too many fields: expected "+D.length+" fields but parsed "+g,E+O):g<D.length&&T("FieldMismatch","TooFewFields","Too few fields: expected "+D.length+" fields but parsed "+g,E+O)),L}var d=1;return!i.data.length||Array.isArray(i.data[0])?(i.data=i.data.map(u),d=i.data.length):i.data=u(i.data,0),t.header&&i.meta&&(i.meta.fields=D),E+=d,i}()}function I(){return t.header&&D.length===0}function x(u,d){return v=u,t.dynamicTypingFunction&&t.dynamicTyping[v]===void 0&&(t.dynamicTyping[v]=t.dynamicTypingFunction(v)),(t.dynamicTyping[v]||t.dynamicTyping)===!0?d==="true"||d==="TRUE"||d!=="false"&&d!=="FALSE"&&(function(O){if(l.test(O)){var g=parseFloat(O);if(f<g&&g<s)return!0}return!1}(d)?parseFloat(d):A.test(d)?new Date(d):d===""?null:d):d;var v}function T(u,d,v,O){var g={type:u,code:d,message:v};O!==void 0&&(g.row=O),i.errors.push(g)}this.parse=function(u,d,v){var O=t.quoteChar||'"';if(t.newline||(t.newline=function(C,m){C=C.substring(0,1048576);var B=new RegExp(j(m)+"([^]*?)"+j(m),"gm"),z=(C=C.replace(B,"")).split("\r"),W=C.split(`
`),Q=1<W.length&&W[0].length<z[0].length;if(z.length===1||Q)return`
`;for(var P=0,k=0;k<z.length;k++)z[k][0]===`
`&&P++;return P>=z.length/2?`\r
`:"\r"}(u,O)),n=!1,t.delimiter)p(t.delimiter)&&(t.delimiter=t.delimiter(u),i.meta.delimiter=t.delimiter);else{var g=function(C,m,B,z,W){var Q,P,k,R;W=W||[",","	","|",";",o.RECORD_SEP,o.UNIT_SEP];for(var re=0;re<W.length;re++){var c=W[re],ue=0,J=0,ne=0;k=void 0;for(var $=new ee({comments:z,delimiter:c,newline:m,preview:10}).parse(C),G=0;G<$.data.length;G++)if(B&&S($.data[G]))ne++;else{var Z=$.data[G].length;J+=Z,k!==void 0?0<Z&&(ue+=Math.abs(Z-k),k=Z):k=Z}0<$.data.length&&(J/=$.data.length-ne),(P===void 0||ue<=P)&&(R===void 0||R<J)&&1.99<J&&(P=ue,Q=c,R=J)}return{successful:!!(t.delimiter=Q),bestDelimiter:Q}}(u,t.newline,t.skipEmptyLines,t.comments,t.delimitersToGuess);g.successful?t.delimiter=g.bestDelimiter:(n=!0,t.delimiter=o.DefaultDelimiter),i.meta.delimiter=t.delimiter}var L=ge(t);return t.preview&&t.header&&L.preview++,e=u,r=new ee(L),i=r.parse(e,d,v),b(),h?{meta:{paused:!0}}:i||{meta:{paused:!1}}},this.paused=function(){return h},this.pause=function(){h=!0,r.abort(),e=p(t.chunk)?"":e.substring(r.getCharIndex())},this.resume=function(){w.streamer._halted?(h=!1,w.streamer.parseChunk(e,!0)):setTimeout(w.resume,3)},this.aborted=function(){return q},this.abort=function(){q=!0,r.abort(),i.meta.aborted=!0,p(t.complete)&&t.complete(i),e=""}}function j(t){return t.replace(/[.*+?^${}()|[\]\\]/g,"\\$&")}function ee(t){var e,r=(t=t||{}).delimiter,n=t.newline,s=t.comments,f=t.step,l=t.preview,A=t.fastMode,w=e=t.quoteChar===void 0||t.quoteChar===null?'"':t.quoteChar;if(t.escapeChar!==void 0&&(w=t.escapeChar),(typeof r!="string"||-1<o.BAD_DELIMITERS.indexOf(r))&&(r=","),s===r)throw new Error("Comment character same as delimiter");s===!0?s="#":(typeof s!="string"||-1<o.BAD_DELIMITERS.indexOf(s))&&(s=!1),n!==`
`&&n!=="\r"&&n!==`\r
`&&(n=`
`);var a=0,E=!1;this.parse=function(h,q,D){if(typeof h!="string")throw new Error("Input must be a string");var i=h.length,y=r.length,S=n.length,b=s.length,I=p(f),x=[],T=[],u=[],d=a=0;if(!h)return N();if(t.header&&!q){var v=h.split(n)[0].split(r),O=[],g={},L=!1;for(var C in v){var m=v[C];p(t.transformHeader)&&(m=t.transformHeader(m,C));var B=m,z=g[m]||0;for(0<z&&(L=!0,B=m+"_"+z),g[m]=z+1;O.includes(B);)B=B+"_"+z;O.push(B)}if(L){var W=h.split(n);W[0]=O.join(r),h=W.join(n)}}if(A||A!==!1&&h.indexOf(e)===-1){for(var Q=h.split(n),P=0;P<Q.length;P++){if(u=Q[P],a+=u.length,P!==Q.length-1)a+=n.length;else if(D)return N();if(!s||u.substring(0,b)!==s){if(I){if(x=[],ne(u.split(r)),de(),E)return N()}else ne(u.split(r));if(l&&l<=P)return x=x.slice(0,l),N(!0)}}return N()}for(var k=h.indexOf(r,a),R=h.indexOf(n,a),re=new RegExp(j(w)+j(e),"g"),c=h.indexOf(e,a);;)if(h[a]!==e)if(s&&u.length===0&&h.substring(a,a+b)===s){if(R===-1)return N();a=R+S,R=h.indexOf(n,a),k=h.indexOf(r,a)}else if(k!==-1&&(k<R||R===-1))u.push(h.substring(a,k)),a=k+y,k=h.indexOf(r,a);else{if(R===-1)break;if(u.push(h.substring(a,R)),Z(R+S),I&&(de(),E))return N();if(l&&x.length>=l)return N(!0)}else for(c=a,a++;;){if((c=h.indexOf(e,c+1))===-1)return D||T.push({type:"Quotes",code:"MissingQuotes",message:"Quoted field unterminated",row:x.length,index:a}),G();if(c===i-1)return G(h.substring(a,c).replace(re,e));if(e!==w||h[c+1]!==w){if(e===w||c===0||h[c-1]!==w){k!==-1&&k<c+1&&(k=h.indexOf(r,c+1)),R!==-1&&R<c+1&&(R=h.indexOf(n,c+1));var ue=$(R===-1?k:Math.min(k,R));if(h.substr(c+1+ue,y)===r){u.push(h.substring(a,c).replace(re,e)),h[a=c+1+ue+y]!==e&&(c=h.indexOf(e,a)),k=h.indexOf(r,a),R=h.indexOf(n,a);break}var J=$(R);if(h.substring(c+1+J,c+1+J+S)===n){if(u.push(h.substring(a,c).replace(re,e)),Z(c+1+J+S),k=h.indexOf(r,a),c=h.indexOf(e,a),I&&(de(),E))return N();if(l&&x.length>=l)return N(!0);break}T.push({type:"Quotes",code:"InvalidQuotes",message:"Trailing quote on quoted field is malformed",row:x.length,index:a}),c++}}else c++}return G();function ne(U){x.push(U),d=a}function $(U){var ke=0;if(U!==-1){var _e=h.substring(c+1,U);_e&&_e.trim()===""&&(ke=_e.length)}return ke}function G(U){return D||(U===void 0&&(U=h.substring(a)),u.push(U),a=i,ne(u),I&&de()),N()}function Z(U){a=U,ne(u),u=[],R=h.indexOf(n,a)}function N(U){return{data:x,errors:T,meta:{delimiter:r,linebreak:n,aborted:E,truncated:!!U,cursor:d+(q||0)}}}function de(){f(N()),x=[],T=[]}},this.abort=function(){E=!0},this.getCharIndex=function(){return a}}function H(t){var e=t.data,r=M[e.workerId],n=!1;if(e.error)r.userError(e.error,e.file);else if(e.results&&e.results.data){var s={abort:function(){n=!0,te(e.workerId,{data:[],errors:[],meta:{aborted:!0}})},pause:ve,resume:ve};if(p(r.userStep)){for(var f=0;f<e.results.data.length&&(r.userStep({data:e.results.data[f],errors:e.results.errors,meta:e.results.meta},s),!n);f++);delete e.results}else p(r.userChunk)&&(r.userChunk(e.results,s,e.file),delete e.results)}e.finished&&!n&&te(e.workerId,e.results)}function te(t,e){var r=M[t];p(r.userComplete)&&r.userComplete(e),r.terminate(),delete M[t]}function ve(){throw new Error("Not implemented.")}function ge(t){if(typeof t!="object"||t===null)return t;var e=Array.isArray(t)?[]:{};for(var r in t)e[r]=ge(t[r]);return e}function V(t,e){return function(){t.apply(e,arguments)}}function p(t){return typeof t=="function"}return le&&(_.onmessage=function(t){var e=t.data;if(o.WORKER_ID===void 0&&e&&(o.WORKER_ID=e.workerId),typeof e.input=="string")_.postMessage({workerId:o.WORKER_ID,results:o.parse(e.input,e.config),finished:!0});else if(_.File&&e.input instanceof File||e.input instanceof Object){var r=o.parse(e.input,e.config);r&&_.postMessage({workerId:o.WORKER_ID,results:r,finished:!0})}}),(ie.prototype=Object.create(F.prototype)).constructor=ie,(se.prototype=Object.create(F.prototype)).constructor=se,(X.prototype=Object.create(X.prototype)).constructor=X,(ae.prototype=Object.create(F.prototype)).constructor=ae,o})})(me);var ye=me.exports;self.onmessage=ce=>{const pe=ce.data;let he=0,_={columnNames:[],columnTypes:[],secondRowValues:[],previewRowCount:0,columnCount:0};const Y=()=>new Promise((M,fe)=>{ye.parse(pe,{worker:!0,step:()=>{he++},complete:()=>{M(he)},error:o=>{fe(o)}})});new Promise((M,fe)=>{ye.parse(pe,{header:!0,dynamicTyping:!0,worker:!0,preview:2,complete:o=>{const K=o.meta.fields||[],F=[],ie=o.data[1]||{},se=o.data.length,X=K.length;K.forEach(oe=>{const j=o.data.map(H=>H[oe]),ee=new Set(j.map(H=>typeof H));if(ee.has("number")){const H=j.some(te=>Number(te)!==Math.floor(Number(te)));F.push(H?"float":"integer")}else if(ee.has("boolean"))F.push("boolean");else if(ee.has("string")){const H=j.every(te=>!Number.isNaN(Date.parse(te)));F.push(H?"date":"string")}else F.push("unknown")});const ae=K.map(oe=>ie[oe]);_={columnNames:K,columnTypes:F,secondRowValues:ae,previewRowCount:se,columnCount:X},M()},error:o=>{fe(o)}})}).then(()=>Y()).then(M=>{_.rowCount=M,self.postMessage(_)}).catch(M=>{self.postMessage({error:M.message})})}})();
//# sourceMappingURL=csvWorker-e8cMUEqV.js.map