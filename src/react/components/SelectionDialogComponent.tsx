import { useDimensionFilter, useParamColumns, useRangeFilter } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Slider } from "@mui/material";
import { useEffect, useState } from "react";

const TextComponent = ({column} : {column: DataColumn<CategoricalDataType>}) => {
    const dim = useDimensionFilter(column);
    const [value, setValue] = useState("");
    useEffect(() => {
        // dim.filter("filterCategories", [column.name], [value], true);
        dim.removeFilter();
    }, [dim, value, column.name]);
    return (
        <div className="bg-red-500">
            WIP...
            <input type="text" value={value} onChange={(e) => setValue(e.target.value)} />
        </div>
    );
}

const UniqueComponent = ({column} : {column: DataColumn<"unique">}) => {
    return (
        <div>
            <h2>Unique:</h2>
            {column.name} : {column.datatype}
        </div>
    );
}

const NumberComponent = ({column} : {column: DataColumn<NumberDataType>}) => {
    const { min, setMin, max, setMax } = useRangeFilter(column);

    return (
        <div>
            <Slider 
            value={[min, max]}
            min={column.minMax[0]}
            max={column.minMax[1]}
            onChangeCommitted={(_, newValue) => {
                setMin(newValue[0]);
                setMax(newValue[1]);
            }}
            />
            <div>
                <input className="" type="number" value={min} onChange={(e) => setMin(Number(e.target.value))} />
                <input className="float-right" type="number" value={max} onChange={(e) => setMax(Number(e.target.value))} />
            </div>
        </div>
    );
}

const Components: {
    [K in DataType]: React.FC<{column: DataColumn<K>}>;
} = {
    integer: NumberComponent,
    double: NumberComponent,
    int32: NumberComponent,
    text: TextComponent,
    text16: TextComponent,
    multitext: TextComponent,
    unique: UniqueComponent,
}

function AbstractComponent<T extends DataType>({column} : {column: DataColumn<T>}) {
    const Component = Components[column.datatype] as React.FC<{column: DataColumn<T>}>;
    return (
        <div>
            <h2>{column.name}</h2>
            <Component column={column} />
        </div>
    );
}


export default function SelectionDialogComponent() {
    const cols = useParamColumns();
    
    return (
        <div className="p-5">
        {cols.map((col) => <AbstractComponent key={col.field} column={col} />)}
        </div>
    );
}
