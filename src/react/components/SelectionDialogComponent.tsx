import { useConfig, useDimensionFilter, useParamColumns } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Autocomplete, Button, Checkbox, Chip, Slider, TextField } from "@mui/material";
import { useCallback, useEffect, useState } from "react";

import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import type { SelectionDialogConfig, CategoryFilter } from "./SelectionDialogReact";
import { observer } from "mobx-react-lite";
import { action } from "mobx";
import { useChart } from "../context";

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

type Props<K extends DataType> = {
    column: DataColumn<K>;
    // filter not passed in as a possibly null prop...
    // each component will extract it from the config.filters object, and write it back there.
    // mobx will handle the reactivity, and it should get serialized appropriately.
    // filter?: K extends CategoricalDataType ? CategoryFilter : [number, number];
};


const TextComponent = ({column} : Props<CategoricalDataType>) => {
    const dim = useDimensionFilter(column);
    const filters = useConfig<SelectionDialogConfig>().filters;
    const { values } = column;
    const filter = filters[column.field] as CategoryFilter | null;
    const value = filter?.category || [];
    const setValue = useCallback((newValue: string[]) => {
        const newFilter = filter || { category: [] };
        action(() => {
            newFilter.category = newValue;
            filters[column.field] = newFilter;
        })();
    }, [filters, column.field, filter]);
    // react to changes in value and update the filter
    useEffect(() => {
        if (filter === null) {
            dim.removeFilter();
            return;
        }
        dim.filter("filterCategories", [column.field], value, true);
    }, [dim, column.field, value, filter]);
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

const UniqueComponent = ({column} : Props<"unique">) => {
    //todo...this is also not handled by original version.
    return (
        <div>
            Filtering of unique values not yet implemented.
        </div>
    );
}

/**
 * This was exposed as a more general-purpose hook with useState, but moved to here to handle
 * state with mobx in the config.filters object.
 */
function useRangeFilter(column: DataColumn<NumberDataType>) {
    const filter = useDimensionFilter(column);
    const filters = useConfig<SelectionDialogConfig>().filters;
    const value = (filters[column.field] || column.minMax) as [number, number];
    const isInteger = column.datatype.match(/int/);
    const step = isInteger ? 1 : 0.01;
    useEffect(() => {
        // filter.removeFilter();
        const [min, max] = value;
        filter.filter("filterRange", [column.name], { min, max }, true);
    }, [column, filter, value]);
    return { value, step };
}

const NumberComponent = ({column} : Props<NumberDataType>) => {
    const filters = useConfig<SelectionDialogConfig>().filters;
    const { value, step } = useRangeFilter(column);
    const setValue = useCallback((newValue: [number, number]) => {
        action(() => filters[column.field] = newValue)();
    }, [filters, column.field]);
    return (
        <div>
            <Slider 
            size="small"
            value={value}
            min={column.minMax[0]}
            max={column.minMax[1]}
            step={step}
            onChange={(_, newValue) => setValue(newValue as [number, number]) }
            />
            <div>
                <TextField size="small" className="max-w-20" type="number" 
                value={value[0]} 
                onChange={(e) => setValue([Number(e.target.value), value[1]])} />
                <TextField size="small" className="max-w-20 float-right" type="number" 
                value={value[1]} 
                onChange={(e) => setValue([value[0], Number(e.target.value)])} />
            </div>
        </div>
    );
}

const Components: {
    [K in DataType]: React.FC<Props<K>>;
} = {
    integer: observer(NumberComponent),
    double: observer(NumberComponent),
    int32: observer(NumberComponent),
    text: observer(TextComponent),
    text16: observer(TextComponent),
    multitext: observer(TextComponent), //should have different behaviour in the future
    unique: observer(UniqueComponent),
}

const AbstractComponent = observer(function AbstractComponent<K extends DataType>({column} : Props<K>) {
    const Component = Components[column.datatype] as React.FC<Props<K>>;
    //todo: consider reset (& delete / active toggle?) for each filter
    const filters = useConfig<SelectionDialogConfig>().filters;
    const hasFilter = filters[column.field] !== null;
    const clearFilter = action(() => filters[column.field] = null);
    return (
        <div>
            <h3>{column.name}</h3>
            {hasFilter && <Button onClick={clearFilter}>Clear</Button>}
            <Component column={column} />
        </div>
    );
});

/**
 * This will control the behaviour of the reset menuIcon in the chart header - not rendered with react.
 */
function useResetButton() {
    const chart = useChart();
    const filters = useConfig<SelectionDialogConfig>().filters;
    const hasFilter = Object.values(filters).some((f) => f !== null);
    useEffect(() => {
        console.log("hasFilter changed (in hook): ", hasFilter);
        chart.resetButton.style.display = hasFilter ? "inline" : "none";
    }, [hasFilter, chart.resetButton]);
}

export default function SelectionDialogComponent() {
    const cols = useParamColumns();
    useResetButton();
    // todo ability to dynamically add or remove columns
    // maybe arranged in a hierarchy with a tree view?
    return (
        <div className="p-3">
        {cols.map((col) => <AbstractComponent key={col.field} column={col} />)}
        </div>
    );
}
