function getRandomString(len=6,an){
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

export {getRandomString}