import { useDimensionFilter, useParamColumns, useRangeFilter } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Autocomplete, Checkbox, Chip, Slider, TextField } from "@mui/material";
import { useCallback, useEffect, useState } from "react";

import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

const TextComponent = ({column} : {column: DataColumn<CategoricalDataType>}) => {
    const dim = useDimensionFilter(column);
    const { values } = column;
    const [value, setValue] = useState<string[]>([]);
    useEffect(() => {
        dim.filter("filterCategories", [column.name], value, true);
    }, [dim, value, column.name]);
    const [hasFocus, setHasFocus] = useState(false);
    const toggleOption = useCallback((option: string) => {
        if (value.includes(option)) {
            setValue(value.filter((v) => v !== option));
            return;
        }
        setValue([...value, option]);
    }, [value]);
    const selectAll = useCallback(() => {
        setValue([...values]);
    }, [values]);
    const toggleSelection = useCallback(() => {
        const newValues = values.filter((v) => !value.includes(v));
        setValue(newValues);
    }, [value, values]);

    return (
        <Autocomplete 
            multiple
            size="small"
            options={values}
            value={value}
            onChange={(_, newValue) => setValue(newValue)}
            onFocus={() => setHasFocus(true)}
            onBlur={() => setHasFocus(false)}
            renderInput={(props) => {
                const { key, ...p } = props as typeof props & {
                    key: string;
                }; //questionable mui types?
                return <>
                    {hasFocus && <Button onClick={selectAll}>All</Button>}
                    {hasFocus && <Button onClick={toggleSelection}>Toggle</Button>}
                    <TextField key={key} {...p} />
                </>
            }}
            renderOption={(props, option) => {
                const { key, ...optionProps } = props as typeof props & {
                    key: string;
                }; //questionable mui types?
                return (
                    <li
                        key={key}
                        {...optionProps}
                        onClick={() => toggleOption(option)}
                    >
                        <Checkbox
                            icon={icon}
                            checkedIcon={checkedIcon}
                            style={{ marginRight: 8 }}
                            checked={value.includes(option)}
                            onClick={() => toggleOption(option)}
                        />
                        {option}
                    </li>
                );
            }}
            renderTags={(value, getTagProps) => {
                //seems to be a material-ui bug with not properly handling key / props...
                //https://stackoverflow.com/questions/75818761/material-ui-autocomplete-warning-a-props-object-containing-a-key-prop-is-be
                const chips = value.map((option, index) => {
                    const { key, ...tagValues } = getTagProps({ index });
                    return <Chip key={key} {...tagValues} label={option} />;
                });
                //checking against length of chips is not ideal
                //similar to MultiSelect.js config.maxShow (default 4)
                //consider 'maxChars' perhaps so that we respond to the actual amount of text
                //or I suppose we could have some CSS container query or something...
                if (chips.length > 8) {
                    return <div>{chips.length} selected</div>
                }
                //todo improve styling - bad from UX perspective at the moment when it overflows.
                return <div className="max-h-32 overflow-auto">{chips}</div>;
            }}
        />
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
            size="small"
            value={[min, max]}
            min={column.minMax[0]}
            max={column.minMax[1]}
            onChange={(_, newValue) => {
                setMin(newValue[0]);
                setMax(newValue[1]);
            }}
            />
            <div>
                <TextField size="small" className="max-w-20" type="number" value={min} onChange={(e) => setMin(Number(e.target.value))} />
                <TextField size="small" className="max-w-20 float-right" type="number" value={max} onChange={(e) => setMax(Number(e.target.value))} />
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
    ///XXX we should be able to save back to config - and use it better!!!
    return (
        <div className="p-5">
        {cols.map((col) => <AbstractComponent key={col.field} column={col} />)}
        </div>
    );
}
