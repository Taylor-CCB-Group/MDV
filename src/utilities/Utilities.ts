import type { CategoricalDataType, DataColumn, DataType, NumberDataType } from "@/charts/charts";
import { isDatatypeCategorical, isDatatypeNumeric } from "@/lib/utils";

export function getRandomString(len=6) {//, an) {
    if (!len) {
        len = 6;
    }
    // an = an?.toLowerCase();
    let str = "";
    let i = 0;
    // const min = an === "a" ? 10 : 0;
    // const max = an === "n" ? 10 : 62;
    const min = 0;
    const max = 62;
    while (i++ < len) {
        let r = (Math.random() * (max - min) + min) << 0;
        str += String.fromCharCode((r += r > 9 ? (r < 36 ? 55 : 61) : 48));
    }
    return str;
}

export function NPOT(n: number) {
    return 2 ** Math.ceil(Math.log2(n));
}

/**
 * Debounce function. Return a function that, as long as it continues to be invoked, 
 * will not be called until {@param timeout}ms have passed since the last call.
 * based on https://www.freecodecamp.org/news/javascript-debounce-example/
 */
export function debounce<A=unknown>(fn: (args: A) => void, timeout: number) {
    let timer;
    return (...args) => {
        clearTimeout(timer);
        timer = setTimeout(() => {
            //@ts-ignore
            fn.apply(this, args);
        }, timeout);
    };
}

/**
 * Type predicate to check if a column is numeric.
 */
export function isColumnNumeric(column: DataColumn<DataType>): column is DataColumn<NumberDataType> {
    const t = column.datatype;
    // return t === "double" || t === "int32" || t === "integer";
    return isDatatypeNumeric(t);
}

/**
 * Type predicate to check if a column is categorical.
 */
export function isColumnText(column: DataColumn<DataType>): column is DataColumn<CategoricalDataType> {
    const t = column.datatype;
    // return t === "text" || t === "multitext" || t === "text16";
    return isDatatypeCategorical(t);
}
