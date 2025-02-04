/**
 *  A collection of dataloaders and helper functions that can
 * can be used with {@link ChartManager}
 * @module DataLoaders
 *
 *
 **/

/**
 * 
 * @param data {ArrayBuffer} - an array buffer containing raw concatenated
 * column data
 * @param {object} columns - any array of column objects
 * <ul>
 *   <li> field </li>
 *   <li> datatype </li>
 *   <li> sgtype </li>
 * </ul>
 * @param size {integer} - the size of the columns
 * @returns {object[]} a list of objects containing each colums's field name
 * and data
 **/
function processArrayBuffer(data, columns, size) {
    const dataList = [];
    let offset = 0;
    for (const column of columns) {
        //default values for numbers
        let arrayType = column.datatype === "int32" ? Int32Array : Float32Array;
        //the length of the typed array
        let arr_len = size;
        //the number of bytes for each item in the column's data
        let bytes = 4;
        //set the correct values for other datatypes
        if (column.datatype === "unique" || column.datatype === "text") {
            arrayType = Uint8Array;
            bytes = 1;
            if (column.datatype === "unique") {
                bytes = column.stringLength;
                arr_len = size * bytes;
            }
        } else if (column.datatype === "multitext") {
            arrayType = Uint16Array;
            if (!column.stringLength) {
                // this seems to be the actual logic needed, at least for tag data saved by the client
                bytes = 2;
                arr_len = size;
            } else {
                // not sure if/when this applies...
                bytes = column.stringLength * 2; //stringLength is 'undefined' when text16 version of multitext is saved by client. there is no 'values' key at all.
                arr_len = size * column.stringLength;
            }
        } else if (column.datatype === "text16") {
            arrayType = Uint16Array;
            bytes = 2;
            arr_len = size;
        }
        //special way to deal with sparse data
        //assumes all sparse data is integer/double
        //the first 4 bytes specifies the number of values (n)
        //The next n*4 bytes are the indexes of these values (Uint32)
        //the next n*4 bytes are the actual values (Float32)
        if (column.sgtype === "sparse") {
            //first byte is length of data
            const l = new Uint32Array(data, offset, 1)[0];
            offset += 4;
            //get the indexes
            const indexes = new Uint32Array(data, offset, l);
            offset += l * 4;
            //get the values
            const values = new Float32Array(data, offset, l);
            offset += l * 4;
            const sab = new SharedArrayBuffer(size * 4);
            const new_arr = new Float32Array(sab);
            //fill array with missing values
            new_arr.fill(Number.NaN);
            //add the sparse data
            for (let i = 0; i < indexes.length; i++) {
                new_arr[indexes[i]] = values[i];
            }
            dataList.push({ data: sab, field: column.field });
        } else {
            const len = size * bytes;
            //get the data from the arraybuffer into a SharedArrayBuffer
            //unfortunately cannot be copied directly-  have to via a typed Array
            const arr = new arrayType(data, offset, arr_len);
            const sab = new SharedArrayBuffer(len);
            const new_arr = new arrayType(sab);
            new_arr.set(arr, 0);
            dataList.push({ data: sab, field: column.field });
            offset += len;
        }
    }
    return dataList;
}

/**
 * Gets bytes from an API. The data loader will send a post request
 * to the url with with a jsonified object containing the datasource
 * and column information
 * <pre>
 * {
 *    "data_source":"mydataource"
 *    "columns":[{"field":"x1","datatype":"integer"}]
 * }
 * </pre>
 * returns a dataloader
 * 
 * @param url {string} - The url of the api
 * @returns {function} a dataloader that can be used to construct {@link ChartManager}
 **/
function getArrayBufferDataLoader(url, decompress = false) {
    return async (columns, dataSource, size) => {
        //get the data
        try {
            const response = await fetch(url, {
                method: "POST",
                body: JSON.stringify({ columns: columns, data_source: dataSource }),
                headers: {
                    "Content-Type": "application/json",
                },
            });
            //the data is any arraybuffer containing each individual
            //column's raw data
            let data = await response.arrayBuffer();
            data = decompress ? await decompressData(data) : data;
            return processArrayBuffer(data, columns, size);
        } catch (e) {
            console.warn(`Error fetching data: ${e}`);
        }
    };
}

/**
 * Gets bytes from the server form static compressed binary files
 * The files will be in in the supplied folder , one per datsource
 * named dsname1.gz dsname2.gz etc.
 *
 * returns a dataLoader
 * 
 * @param dataSources - The url of the remote folder
 * @param folder 
 * @returns {function} a dataloader that can be used to construct {@link ChartManager}
 **/

function getLocalCompressedBinaryDataLoader(dataSources, folder) {
    const loaders = {};
    for (const ds of dataSources) {
        loaders[ds.name] = new CompressedBinaryDataLoader(
            `${folder}/${ds.name}.gz`,
            ds.size,
        );
    }
    return async (columns, dataSource, size) => {
        return await loaders[dataSource].getColumnData(columns, size);
    };
}

class CompressedBinaryDataLoader {
    constructor(url, size) {
        this.url = url;
        this.index = null;
        this.size = size;
    }
    async getColumnData(cols, size) {
        const { default: pako } = await import("pako");
        if (!this.index) {
            const iurl = this.url.replace(".gz", ".json");
            const resp = await fetch(iurl);
            this.index = await resp.json();
        }
        //unable to request multiple ranges, which would be more efficient
        //send off request for each column separately and then return all of them
        return await Promise.all(
            cols.map(async (c) => {
                let lu = c.field;
                if (c.subgroup) {
                    const arr = c.field.split("|");
                    lu = arr[0] + arr[2];
                }
                const i = this.index[lu];

                let resp = null;
                const args = {
                    headers: {
                        responseType: "arraybuffer",
                        range: `bytes=${i[0]}-${i[1]}`,
                    },
                };
                try {
                    resp = await fetch(this.url, args);
                } catch (e) {
                    //lets try again
                    resp = await fetch(this.url, args);
                }
                const bytes = await resp.arrayBuffer();
                const expectedLength = i[1] - i[0] + 1;
                if (bytes.byteLength !== expectedLength) {
                    console.warn(
                        `Expected ${expectedLength} bytes but got ${bytes.byteLength} for ${c.field}... is the server correctly configured to Accept-Ranges?`,
                    );
                }
                try {
                    const output = pako.inflate(bytes);

                    if (c.sgtype === "sparse") {
                        const b = output.buffer;
                        const l = new Uint32Array(b, 0, 1)[0];
                        //get the indexes
                        const indexes = new Uint32Array(b, 4, l);
                        //get the values
                        const values = new Float32Array(b, l * 4 + 4, l);
                        const sb = new SharedArrayBuffer(size * 4);
                        const new_arr = new Float32Array(sb);
                        //fill array with missing values
                        new_arr.fill(Number.NaN);
                        for (let i = 0; i < indexes.length; i++) {
                            new_arr[indexes[i]] = values[i];
                        }
                        return { data: sb, field: c.field };
                    }

                    const sb = new SharedArrayBuffer(output.length);
                    const f = new Uint8Array(sb);
                    f.set(output, 0);
                    return { data: sb, field: c.field };
                } catch (e) {
                    console.warn(`Error inflating ${c.field}: ${e}`);
                }
            }),
        );
    }
}

async function decompressData(buffer) {
    const { default: pako } = await import("pako");
    return pako.inflate(buffer).buffer;
}

export {
    getArrayBufferDataLoader,
    processArrayBuffer,
    getLocalCompressedBinaryDataLoader,
    decompressData,
};
