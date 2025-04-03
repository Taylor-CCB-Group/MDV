const datatypes = {
    int32: { arr: Int32Array, type: "number" },
    double: { arr: Float32Array, type: "number" }, // why is this called double?
    integer: { arr: Int32Array, type: "number" },
    text16: { arr: Uint16Array, type: "text" },
    text: { arr: Uint8Array, type: "text" },
    unique: { arr: Uint8Array, type: "unique" },
    multitext: { arr: Uint8Array, type: "multitext" },
} as const;
/**
 * The are the names used to refer to the types of data can be stored in a column.
 */
// export type DataType = keyof typeof datatypes;
/**
 * todo
 * Associates a {@link DataType} with the corresponding `TypedArray` type that will be used to store the data.
 */
// export type DataStructureTypes = {
//     [T in DataType]: typeof datatypes[T]["arr"];
// }
// type TextData = DataStructureTypes["text"];
export default datatypes;
