import { useConfig, useDimensionFilter, useParamColumns, useRangeFilter } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Autocomplete, Button, Checkbox, Chip, Slider, TextField } from "@mui/material";
import { useCallback, useEffect, useState } from "react";

import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import type { SelectionDialogConfig, CategoryFilter } from "./SelectionDialogReact";
import { observer } from "mobx-react-lite";
import { action } from "mobx";

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

type P<K extends DataType> = {
    column: DataColumn<K>;
    // filter not passed in as a possibly null prop...
    // each component will extract it from the config.filters object, and write it back there.
    // mobx will handle the reactivity, and it should get serialized appropriately.
    // filter?: K extends CategoricalDataType ? CategoryFilter : [number, number];
};


const TextComponent = ({column} : P<CategoricalDataType>) => {
    const dim = useDimensionFilter(column);
    const filters = useConfig<SelectionDialogConfig>().filters;
    const { values } = column;
    const filter = filters[column.field] as CategoryFilter | null;
    const value = filter?.category || [];
    const setValue = useCallback((newValue: string[]) => {
        const newFilter = filter || { category: [] };
        newFilter.category = newValue;
        action(() => filters[column.field] = newFilter)();
    }, [filters, column.field, filter]);
    // react to changes in value and update the filter
    useEffect(() => {
        dim.filter("filterCategories", [column.field], value, true);
    }, [dim, column.field, value]);
    const [hasFocus, setHasFocus] = useState(false);
    const toggleOption = useCallback((option: string) => {
        if (value.includes(option)) {
            setValue(value.filter((v) => v !== option));
            return;
        }
        setValue([...value, option]);
    }, [setValue, value]);
    const selectAll = useCallback(() => {
        setValue([...values]);
    }, [values, setValue]);
    const toggleSelection = useCallback(() => {
        const newValues = values.filter((v) => !value.includes(v));
        setValue(newValues);
    }, [values, value, setValue]);
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
    //todo...
    return (
        <div>
            <TextField size="small" />
        </div>
    );
}

const NumberComponent = ({column} : {column: DataColumn<NumberDataType>}) => {
    const filters = useConfig<SelectionDialogConfig>().filters;
    const { min, setMin, max, setMax, step } = useRangeFilter(column);
    // we need to recognise serialised values from the config object...
    // maybe the hook with useState isn't ideal here...
    useEffect(() => {
        if (filters[column.field] && Array.isArray(filters[column.field])) {
            const [min, max] = filters[column.field] as [number, number];
            setMin(min);
            setMax(max);
        }
    }, [filters, column.field, setMin, setMax]);
    return (
        <div>
            <Slider 
            size="small"
            value={[min, max]}
            min={column.minMax[0]}
            max={column.minMax[1]}
            step={step}
            onChange={(_, newValue) => {
                setMin(newValue[0]);
                setMax(newValue[1]);
                action(() => filters[column.field] = [newValue[0], newValue[1]])();
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
    [K in DataType]: React.FC<P<K>>;
} = {
    integer: observer(NumberComponent),
    double: observer(NumberComponent),
    int32: observer(NumberComponent),
    text: observer(TextComponent),
    text16: observer(TextComponent),
    multitext: observer(TextComponent),
    unique: observer(UniqueComponent),
}

function AbstractComponent<K extends DataType>({column} : P<K>) {
    const Component = Components[column.datatype] as React.FC<P<K>>;
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
