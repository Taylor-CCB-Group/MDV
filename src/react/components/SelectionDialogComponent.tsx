import { useCloseOnIntersection, useConfig, useDimensionFilter, useParamColumnsExperimental } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Accordion, AccordionDetails, AccordionSummary, Autocomplete, Checkbox, Chip, IconButton, TextField, Typography } from "@mui/material";
import { createFilterOptions } from '@mui/material/Autocomplete';
import { type MouseEvent, useCallback, useEffect, useState, useMemo, useRef, useId } from "react";

import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';
import DoneAllIcon from '@mui/icons-material/DoneAll';
import CachedIcon from '@mui/icons-material/Cached';
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import SwapHorizIcon from '@mui/icons-material/SwapHoriz';
import HighlightOffIcon from '@mui/icons-material/HighlightOff';
import type { SelectionDialogConfig, CategoryFilter, MultiTextFilter, UniqueFilter, RangeFilter } from "./SelectionDialogReact";
import { observer } from "mobx-react-lite";
import { action, runInAction } from "mobx";
import { useChart, useDataStore } from "../context";
import ColumnSelectionComponent from "./ColumnSelectionComponent";
import { useDebounce } from "use-debounce";
import { useHighlightedForeignRowsAsColumns, useRowsAsColumnsLinks } from "../chartLinkHooks";
import { ErrorBoundary } from "react-error-boundary";
import ErrorDisplay from "@/charts/dialogs/ErrorDisplay";
import { TextFieldExtended } from "./TextFieldExtended";
import { isArray } from "@/lib/utils";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";
import { Histogram, type set2d, useRangeFilter } from "./HistogramComponent";

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
const filterOptions = createFilterOptions<any>({ limit: 100 });
const TextComponent = observer(({ column }: Props<CategoricalDataType>) => {
    const [open, setOpen] = useState(false);
    const ref = useRef<HTMLInputElement>(null);
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
    useCloseOnIntersection(ref, () => setOpen(false));
    // react to changes in value and update the filter
    useEffect(() => {
        // filter could be undefined - but we previously checked for null causing component to crash
        if ((!filter) || (filter.category.length === 0)) {
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
            open={open}
            onOpen={() => setOpen(true)}
            onClose={() => setOpen(false)}
            ref={ref}
            renderInput={(props) => {
                const { key, ...p } = props as typeof props & {
                    key: string;
                }; //questionable mui types?
                return <>
                    <TextFieldExtended key={key} {...p}
                        customEndAdornment={(
                            <>
                                <IconButton size="small" aria-label="toggle selection" onClick={toggleSelection}><SwapHorizIcon fontSize="inherit" /></IconButton>
                                <IconButton size="small" aria-label="select all" onClick={selectAll}><DoneAllIcon fontSize="inherit" /></IconButton>
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
        // PaperComponent={(paperProps) => <Paper ref={ref} {...paperProps} />}
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


const NumberComponent = observer(({ column }: Props<NumberDataType>) => {
    const filters = useConfig<SelectionDialogConfig>().filters;
    // revisiting this, still want a thing that is based on `column`, 
    // but the flow of state - i.e. `histogramData` - is changing.
    const rangeProps = useRangeFilter(column);
    const { value, step } = rangeProps;
    const [min, max] = column.minMax;
    const setValue = useCallback<set2d>((newValue) => {
        if (newValue) {
            // constrain with min, max and step
            newValue[0] = Math.round(newValue[0] / step) * step;
            newValue[1] = Math.round(newValue[1] / step) * step;
            if (newValue[0] < min) newValue[0] = min;
            if (newValue[1] > max) newValue[1] = max;
        }
        action(() => filters[column.field] = newValue)();
    }, [filters, column.field, min, max, step]);
    const low = value ? value[0] : min;
    const high = value ? value[1] : max;
    return (
        <div>
            <Histogram {...rangeProps} name={column.name} setValue={setValue} domain={column.minMax} bins={200} histoHeight={100} />
            <div>
                <TextField size="small" className="max-w-20" type="number"
                    variant="standard"
                    value={low}
                    onChange={(e) => setValue([Number(e.target.value), high])}
                />
                <TextField size="small" className="max-w-20 float-right" type="number"
                    variant="standard"
                    value={high}
                    onChange={(e) => setValue([low, Number(e.target.value)])}
                />
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
    //todo: consider reset (& invert / active toggle?) for each filter
    const config = useConfig<SelectionDialogConfig>();
    const { filters } = config;
    const f = filters[column.field];
    // todo: what about category filters with empty array?
    const hasFilter = (f !== null);
    const [isHovered, setIsHovered] = useState(false);
    const clearFilter = useCallback(
        action((e: MouseEvent) => {
            e.stopPropagation();
            // what about clearing the brush on histogram?
            filters[column.field] = null;
        }),
        []  //xxx: don't need deps here because the mobx references are stable?
    );
    const deleteFilter = useCallback((e: MouseEvent) => {
        e.stopPropagation();
        runInAction(() => {
            delete filters[column.field];
            if (!isArray(config.param)) throw new Error("expected param array");
            config.param = config.param.filter((p) => p !== column.field);
        });
        console.log('Delete item');
    }, [filters, column.field, config]);
    return (
        <Accordion defaultExpanded={true}
            onMouseEnter={() => setIsHovered(true)}
            onMouseLeave={() => setIsHovered(false)}
        >
            <AccordionSummary
                expandIcon={<ArrowDropDownIcon />}
            >
                <IconButton
                    onClick={deleteFilter}
                    aria-label="delete"
                    sx={{
                        position: 'absolute',
                        right: '-18px',
                        top: '-18px',
                        opacity: isHovered ? 1 : 0,
                        visibility: isHovered ? 'visible' : 'hidden',
                        transition: 'opacity 0.3s, visibility 0.3s',
                        '&:focus': {
                            visibility: 'visible',
                            opacity: 1, // Ensure it's fully visible when focused
                        },
                    }}
                >
                    <HighlightOffIcon />
                </IconButton>
                <div className="flex items-center h-4">
                    <Typography variant="subtitle1">{column.name}
                        {hasFilter && <IconButton onClick={clearFilter}>
                            <CachedIcon fontSize="small" />
                        </IconButton>}
                        <em className="opacity-40 ml-1">({column.datatype})</em>
                    </Typography>
                </div>
            </AccordionSummary>
            <AccordionDetails>
                <ErrorBoundary FallbackComponent={
                    ({ error }) => <ErrorDisplay error={error} title="Unexpected Error: please report to developers." />
                }>
                    <Component column={column} />
                </ErrorBoundary>
            </AccordionDetails>
        </Accordion>
    );
});

const AddRowComponent = observer(() => {
    //! warning, as well as being disabled for now, this was not working properly
    //(something about setSelectedColumn was going wrong with other recent changes)
    const config = useConfig<SelectionDialogConfig>();
    const { filters, param } = useConfig<SelectionDialogConfig>();
    if (!isArray(param)) throw new Error("expected param array");
    const setSelectedColumn = useCallback((column: string) => {
        if (!column) return;
        if (param.includes(column)) return;
        runInAction(() => {
            // param.push(column); //doesn't trigger reactivity
            config.param = [...param, column];
            filters[column] = null;
        });
    }, [filters, param, config]);
    return (
        <div className="p-5">
            <ColumnSelectionComponent
                //@ts-expect-error !!! setSelectedColumn needs appropriate type
                setSelectedColumn={setSelectedColumn}
                placeholder="Add a filter column"
                //@ts-expect-error 'exclude' design still under review?
                exclude={param}
                // type="_multi_column:all" //not sure about this...
                type={["text", "text16", "multitext", "unique", "integer", "double", "int32"]}
            />
        </div>
    )
})

const ForeignRows = observer(() => {
    // const rlink = useRowsAsColumnsLinks();
    // //!breaking rule of hooks here, but in a way that should be ok at runtime as of now
    // //! (just testing "infinte loop with no link" fix)
    // if (rlink.length === 0) return null; //todo: 30sec video clip
    const [filter, setFilter] = useState("");
    const [max, setMax] = useState(10);
    const [debouncedFilter] = useDebounce(filter, 300);
    const rlink = useRowsAsColumnsLinks();
    const fcols = useHighlightedForeignRowsAsColumns(max, debouncedFilter);
    if (!rlink[0]) return null;
    const { linkedDs, link } = rlink[0];
    return (
        <div className="p-3">
            <Typography variant="h6" sx={{ marginBottom: '0.5em' }}>Columns associated with selected '{linkedDs.name}':</Typography>
            <TextField size="small" label="Filter" variant="outlined" onChange={e => setFilter(e.target.value)} />
            <TextField size="small" className="max-w-20 float-right" type="number"
                label="Max" variant="outlined"
                value={max}
                onChange={(e) => setMax(Number(e.target.value))}
            />

            {fcols.map(col => <AbstractComponent key={col.field} column={col} />)}
        </div>
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
    // we also need to respond to the 'reset all' button
    const ds = useDataStore();
    const id = useId();
    useEffect(() => {
        const k = `SelectionDialog-${id}`;
        const resetAll = (type: string, data: string) => {
            if (type !== "filtered" || data !== "all_removed") return;
            runInAction(() => {
                for (const key in filters) {
                    delete filters[key];
                }
            });
        };
        ds.addListener(k, resetAll);
        return () => ds.removeListener(k);
    }, [ds, id, filters]);
}

const SelectionDialogComponent = () => {
    //!! this component doesn't update with HMR and introducing another wrapper component makes things worse
    //(currently changes here aren't reflected in the browser, but the rest of the components are
    //if we wrap this, then any change causes whole page to reload)
    //! since this is outside the ErrorBoundary, problems in this hook are not well handled
    //we could consider returning some kind of `Result` object from this hook...
    const cols = useParamColumnsExperimental();
    useResetButton();
    const showAddRow = false;
    return (
        <div className="p-3 absolute w-[100%] h-[100%] overflow-x-hidden overflow-y-auto">
            {cols.map((col) => <AbstractComponent key={col.field} column={col} />)}
            {showAddRow && <ErrorBoundary FallbackComponent={
                ({ error }) =>
                (
                    <div className="col-span-3 w-full">
                        <ErrorComponentReactWrapper
                            error={{ message: error.message, stack: error.stack }}
                            //todo assign proper meta data
                            // extraMetaData={}
                            title="Error displaying 'AddRowComponent'. Click to view details"
                        />
                    </div>
                )
            }>
                <AddRowComponent />
            </ErrorBoundary>}
        </div>
    );
};
export default SelectionDialogComponent;
