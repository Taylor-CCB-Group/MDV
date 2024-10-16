import { useConfig, useDimensionFilter, useParamColumns } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Accordion, AccordionDetails, AccordionSummary, Autocomplete, Button, Checkbox, Chip, IconButton, Slider, TextField, type TextFieldProps, Typography } from "@mui/material";
import { createFilterOptions } from '@mui/material/Autocomplete';
import { type MouseEvent, useCallback, useEffect, useState, useMemo } from "react";

import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';
import DoneAllIcon from '@mui/icons-material/DoneAll';
import CachedIcon from '@mui/icons-material/Cached';
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import type { SelectionDialogConfig, CategoryFilter, MultiTextFilter, UniqueFilter, RangeFilter } from "./SelectionDialogReact";
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

function useFilterConfig<K extends DataType>(column: DataColumn<K>) {
    const { datatype, field } = column;
    const filters = useConfig<SelectionDialogConfig>().filters;
    const filter = filters[field] as (K extends "multitext" ? MultiTextFilter
        : K extends "unique" ? UniqueFilter
        : K extends CategoricalDataType ? CategoryFilter
        : RangeFilter) | null;
    return filter;
}

/** Modified version of TextField that allows a `customEndAdornment`
 * along with the standard endAdornment passed in `InputProps`.
 */
const TextFieldExtended = (props: TextFieldProps & { customEndAdornment?: JSX.Element }) => {
    const { InputProps, customEndAdornment, ...rest } = props;
    const inputProps = {
        ...InputProps,
        endAdornment: (
            <>{customEndAdornment} {InputProps.endAdornment}</>
        )
    };
    return <TextField {...rest} InputProps={inputProps} />;
}

const filterOptions = createFilterOptions<any>({ limit: 100 });
const TextComponent = observer(({ column }: Props<CategoricalDataType>) => {
    const dim = useDimensionFilter(column);
    const filters = useConfig<SelectionDialogConfig>().filters;
    const { values } = column;
    const filter = useFilterConfig(column);
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
        if (filter === null || filter.category.length === 0) {
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
            filterOptions={filterOptions}
            onChange={(_, newValue) => setValue(newValue)}
            onFocus={() => setHasFocus(true)}
            onBlur={() => setHasFocus(false)}
            renderInput={(props) => {
                const { key, ...p } = props as typeof props & {
                    key: string;
                }; //questionable mui types?
                return <>
                    <TextFieldExtended key={key} {...p}
                        customEndAdornment={(
                            <>
                                <Button onClick={toggleSelection}>Toggle</Button>
                                <IconButton onClick={selectAll}><DoneAllIcon /></IconButton>
                            </>
                        )}
                    />
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
});

const MultiTextComponent = observer(({ column }: Props<"multitext">) => {
    // todo: think about what to do with null config for filter
    console.log("multitext selection dialog has missing features for 'operand' and other logic");    
    // !!! - uncommenting this stuff makes the entire chart disappear when the filter is removed
    // const config = useFilterConfig(column);
    // const operand = config.operand || "or";
    return (
        <>
            {/* operand: {operand} */}
            <TextComponent column={column} />
        </>
    )
});

const UniqueComponent = observer(({ column }: Props<"unique">) => {
    const dim = useDimensionFilter(column);
    const filters = useConfig<SelectionDialogConfig>().filters;
    const filter = useFilterConfig(column);
    const [initialValue] = useState(filter);
    useEffect(() => {
        if (filter === null) {
            dim.removeFilter();
            return;
        }
        dim.filter("filterUnique", [column.field], filter, true);
        //return () => dim.removeFilter(); //handled by useDimensionFilter
    }, [dim, column.field, filter]);
    const [localFilter, setLocalFilter] = useState(filter || "");
    return (
        <TextField size="small" defaultValue={initialValue}
        onChange={(e) => setLocalFilter(e.target.value)}
        onKeyDown={(e) => {
            if (e.key === "Enter") {
                action(() => filters[column.field] = localFilter)();
            }
        }} />
    );
});

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

const NumberComponent = observer(({ column }: Props<NumberDataType>) => {
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
                onChange={(_, newValue) => setValue(newValue as [number, number])}
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
});

const Components: {
    [K in DataType]: React.FC<Props<K>>;
} = {
    integer: NumberComponent,
    double: NumberComponent,
    int32: NumberComponent,
    text: TextComponent,
    text16: TextComponent,
    multitext: MultiTextComponent,
    unique: UniqueComponent,
}

const AbstractComponent = observer(function AbstractComponent<K extends DataType>({ column }: Props<K>) {
    const Component = Components[column.datatype] as React.FC<Props<K>>;
    //todo: consider reset (& delete / active toggle?) for each filter
    const filters = useConfig<SelectionDialogConfig>().filters;
    const f = filters[column.field];
    // todo: what about category filters with empty array?
    const hasFilter = (f !== null);
    const [defaultExpanded] = useState(hasFilter);
    const clearFilter = useCallback(
        action((e: MouseEvent) => {
            e.stopPropagation();
            filters[column.field] = null;
        }),
        []); //xxx: don't need deps here because the mobx references are stable?
    return (
        <Accordion defaultExpanded={defaultExpanded}>
            <AccordionSummary
                expandIcon={<ArrowDropDownIcon />}
            >
                <div className="flex items-center h-4">
                    <Typography variant="subtitle1">{column.name}</Typography>
                    {hasFilter && <IconButton onClick={clearFilter}>
                        <CachedIcon fontSize="small" />
                    </IconButton>}
                </div>
            </AccordionSummary>
            <AccordionDetails>
                <Component column={column} />
            </AccordionDetails>
        </Accordion>
    );
});

/**
 * This will control the behaviour of the reset menuIcon in the chart header - not rendered with react.
 */
function useResetButton() {
    const chart = useChart();
    const filters = useConfig<SelectionDialogConfig>().filters;
    // todo: what about category filters with empty array?
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
        <div className="p-3 absolute w-[100%] h-[100%] overflow-auto">
            {cols.map((col) => <AbstractComponent key={col.field} column={col} />)}
        </div>
    );
}