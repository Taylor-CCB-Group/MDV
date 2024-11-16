export function getRandomString(len, an) {
    if (!len) {
        len = 6;
    }
    an = an?.toLowerCase();
    let str = "";
    let i = 0;
    const min = an === "a" ? 10 : 0;
    const max = an === "n" ? 10 : 62;
    while (i++ < len) {
        let r = (Math.random() * (max - min) + min) << 0;
        str += String.fromCharCode((r += r > 9 ? (r < 36 ? 55 : 61) : 48));
    }
    return str;
}

export function NPOT(n) {
    return 2 ** Math.ceil(Math.log2(n));
}

//https://www.freecodecamp.org/news/javascript-debounce-example/
export function debounce(fn, timeout) {
    let timer;
    return (...args) => {
        clearTimeout(timer);
        timer = setTimeout(() => {
            fn.apply(this, args);
        }, timeout);
    };
}

// let's make TypeScript versions of these, and also for GuiSpec etc...
export function isColumnNumeric(column) {
    const t = column.datatype;
    return t === "double" || t === "int32" || t === "integer";
}

export function isColumnText(column) {
    const t = column.datatype;
    return t === "text" || t === "multitext" || t === "text16";
}
