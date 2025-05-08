import datatypes from "./datatypes";

// biome-ignore lint/suspicious/noGlobalAssign: <explanation>
onmessage = (e) => {
    const { orderBuffer, columns } = e.data;
    //get access to the order buffer
    const ord = new Uint32Array(orderBuffer);
    //build the methods
    const methods = columns.map((x) => {
        const dinfo = datatypes[x.datatype];
        const data = new dinfo.arr(x.buffer);
        let meth = null;
        if (dinfo.type === "text") {
            meth = getTextSort({ data, desc: x.desc, values: x.values });
        } else if (dinfo.type === "unique") {
            meth = getUniqueSort({ data, desc: x.desc, size: x.stringLength });
        } else if (dinfo.type === "multitext") {
            meth = getMultiTextSort({
                data,
                desc: x.desc,
                size: x.stringLength,
            });
        } else {
            meth = getNumberSort({ data, desc: x.desc });
        }
        return meth;
    });
    ord.sort((a, b) => {
        for (const m of methods) {
            const r = m(a, b);
            if (r !== 0) {
                return r;
            }
        }
        return 0;
    });
    postMessage("done");
};

function getNumberSort({ data, desc }) {
    const mu = desc ? -1 : 1;
    return (a, b) => {
        let va = data[a];
        let vb = data[b];
        va = Number.isNaN(va) ? Number.MAX_VALUE : va;
        vb = Number.isNaN(vb) ? Number.MAX_VALUE : vb;
        return (va > vb ? 1 : va < vb ? -1 : 0) * mu;
    };
}

function getTextSort({ data, values, desc }) {
    //is it quicker to build map and sort on numbers
    //rather than sort by text?
    const map = getSortedValues(values, desc);
    return (a, b) => map[data[a]] - map[data[b]];
}

function getUniqueSort({ data, desc, size }) {
    const names = {};
    const tc = new TextDecoder();
    const mu = order === desc ? -1 : 1;
    for (let i = 0; i < data.length; i++) {
        const index = data[i];
        names[index] = tc.decode(col.slice(index * size, index * size + size));
    }
    return (a, b) => names[a].localeCompare(names[b]) * mu;
}

function getMultiTextSort({ data, desc, values, size }) {
    const map = getSortedValues(values, desc);
    map[65536] = 65356;
    return (a, b) => {
        const a1 = data.slice(a * size, a * size + size);
        const b1 = data.slice(b * size, b * size + size);
        for (let i = 0; i < size; i++) {
            const r = map[a1[i]] - map[b1[i]];
            if (r !== 0) {
                return r;
            }
            //end reached without a difference
            if (a1[i] === 65536 && b1[i] === 65536) {
                return 0;
            }
        }
        return 0;
    };
}

function getSortedValues(values, desc) {
    const mu = desc ? -1 : 1;
    const o = values
        .map((x, i) => [x, i])
        .sort((a, b) => a[0].localeCompare(b[0]) * mu);
    const map = {};
    for (let i = 0; i < o.length; i++) {
        map[o[i][1]] = i;
    }
    return map;
}
