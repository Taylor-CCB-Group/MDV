export function getRandomString(len=6, an){
    if (!len){
        len=6;
    }
    an = an&&an.toLowerCase();
    let str="", i=0, min=an=="a"?10:0, max=an=="n"?10:62;
        for(;i++<len;){
          let r = Math.random()*(max-min)+min <<0;
          str += String.fromCharCode(r+=r>9?r<36?55:61:48);
    }
    return str;
}

export function NPOT(n) {
    return Math.pow(2, Math.ceil(Math.log2(n)))
}

//https://www.freecodecamp.org/news/javascript-debounce-example/
export function debounce(fn, timeout) {
    let timer;
    return (...args) => {
        clearTimeout(timer);
        timer = setTimeout(() => { fn.apply(this, args); }, timeout);
    };
}

export function isColumnNumeric(column) {
    const t = column.datatype;
    return t === "double" || t === "int32" || t === "integer";
}

export function isColumnText(column) {
    const t = column.datatype;
    return t === "text" || t === "multitext" || t === "text16";
}