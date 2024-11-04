import { useConfig, useDimensionFilter, useParamColumnsExperimental } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Accordion, AccordionDetails, AccordionSummary, Autocomplete, Checkbox, Chip, IconButton, Slider, TextField, type TextFieldProps, Typography, Select } from "@mui/material";
import { createFilterOptions } from '@mui/material/Autocomplete';
import { type MouseEvent, useCallback, useEffect, useState, useMemo } from "react";

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
import { useChart } from "../context";
import ColumnSelectionComponent from "./ColumnSelectionComponent";
import type RangeDimension from "@/datastore/RangeDimension";
import { useDebounce } from "use-debounce";
import { useHighlightedForeignRowsAsColumns, useRowsAsColumnsLinks } from "../chartLinkHooks";

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
    const filter = useDimensionFilter(column) as RangeDimension;
    const filters = useConfig<SelectionDialogConfig>().filters;
    const value = (filters[column.field] || column.minMax) as [number, number];
    const isInteger = column.datatype.match(/int/);
    const step = isInteger ? 1 : 0.01;
    const [debouncedValue] = useDebounce(value, 10);

    // Effect to manage the filter state
    useEffect(() => {
        const value = debouncedValue;
        if (value[0] === column.minMax[0] && value[1] === column.minMax[1]) {
            filter.removeFilter();
            return;
        }
        const [min, max] = value;
        filter.filter("filterRange", [column.field], { min, max }, true);
    }, [column, filter, debouncedValue]);

    const [histogram, setHistogram] = useState<number[]>([]);
    const [histogramGpu, setHistogramGpu] = useState<number[]>([]);
    // this could be a more general utility function - expect to extract soon
    const queryHistogram = useCallback(async () => {
        // waste of life trying to use Dimension class.
        // filter.getBinsAsync(column.name, { bins: 100 }).then((histogram) => {
        //     setHistogram(histogram);
        // });
        const worker = new Worker(new URL("../../datastore/rawHistogramWorker.ts", import.meta.url));
        worker.onmessage = (event) => {
            const { histCpu, histGpu } = event.data;
            setHistogram(histCpu);
            setHistogramGpu(histGpu);
            worker.terminate();
        };
        const isInt32 = column.datatype === "int32";
        const originalData = column.data as Float32Array | Int32Array;
        const data = new SharedArrayBuffer(originalData.length * 4);
        new Float32Array(data).set(originalData);
        worker.postMessage({ data, min: column.minMax[0], max: column.minMax[1], bins: 100, isFloat: isInt32 });
    }, [column]);

    const [low, high] = value;
    const [min, max] = column.minMax;
    const lowFraction = (low - min) / (max - min);
    const highFraction = (high - min) / (max - min);

    return { value, step, histogram, histogramGpu, lowFraction, highFraction, queryHistogram };
}
type RangeProps = ReturnType<typeof useRangeFilter>;
const Histogram = observer(({ histogram: data, histogramGpu, lowFraction, highFraction, queryHistogram }: RangeProps) => {
    const prefersDarkMode = window.mdv.chartManager.theme === "dark";
    const width = 99;
    const height = 40;
    const lineColor = prefersDarkMode ? '#fff' : '#000';
    // Find max value for vertical scaling
    const maxValue = Math.max(...data, ...histogramGpu);

    // Define the padding and scaling factor
    const padding = 2;
    const xStep = data.length / (width + 1); // Space between points
    const yScale = (height - 2 * padding) / maxValue; // Scale based on max value

    const lowX = lowFraction * width;
    const highX = highFraction * width;

    useEffect(() => {
        queryHistogram();
    }, [queryHistogram]);

    // Generate the points for the polyline
    // ??? useMemo was wrong ????
    const points = useMemo(() => data.map((value, index) => {
        const x = index * xStep;
        const y = height - padding - value * yScale;
        return `${x},${y}`;
    }).join(' '), [data, xStep, yScale]);

    return (
        <svg width={'100%'} height={height} 
        viewBox={`0 0 ${width} ${height}`}
        preserveAspectRatio="none"
        onClick={queryHistogram}
        >
            {/* Background polyline (the simple line connecting data points) */}
            <polyline
                points={points}
                fill="none"
                stroke={lineColor}
                strokeWidth="1.5"
                // many thanks to ChatGPT for the following line (and the rest of the component
                // but this would have been a real pain to figure out on my own)
                vectorEffect="non-scaling-stroke" // Keeps the stroke width consistent
            />
            <polyline
                points={gpuPoints}
                fill="none"
                stroke={prefersDarkMode ? '#f00' : '#f88'}
                strokeWidth="1.5"
                vectorEffect="non-scaling-stroke"
            />
            {/* Highlighted range */}
            <rect
                x={0}
                y={0}
                width={lowX}
                height={height}
                fill={prefersDarkMode ? '#333' : '#888'}
                fillOpacity="0.8"
            />
            <rect
                x={highX}
                y={0}
                width={width - highX}
                height={height}
                fill={prefersDarkMode ? '#333' : '#888'}
                fillOpacity="0.8"
            />
        </svg>
    );
});

const NumberComponent = observer(({ column }: Props<NumberDataType>) => {
    const filters = useConfig<SelectionDialogConfig>().filters;
    const rangeProps = useRangeFilter(column);
    const { value, step } = rangeProps;
    const setValue = useCallback((newValue: [number, number]) => {
        action(() => filters[column.field] = newValue)();
    }, [filters, column.field]);
    return (
        <div>
            <Histogram {...rangeProps} />
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
    const config = useConfig<SelectionDialogConfig>();
    const { filters } = config;
    const f = filters[column.field];
    // todo: what about category filters with empty array?
    const hasFilter = (f !== null);
    const [defaultExpanded] = useState(hasFilter);
    const [isHovered, setIsHovered] = useState(false);
    const clearFilter = useCallback(
        action((e: MouseEvent) => {
            e.stopPropagation();
            filters[column.field] = null;
        }),
        []  //xxx: don't need deps here because the mobx references are stable?
    );
    const deleteFilter = useCallback((e: MouseEvent) => {
        e.stopPropagation();
        runInAction(() => {
            delete filters[column.field];
            config.param = config.param.filter((p) => p !== column.field);
        });
        console.log('Delete item');
    }, [filters, column.field, config]);
    return (
        <Accordion defaultExpanded={defaultExpanded}
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
                <Component column={column} />
            </AccordionDetails>
        </Accordion>
    );
});

const AddRowComponent = observer(() => {
    const config = useConfig<SelectionDialogConfig>();
    const { filters, param } = useConfig<SelectionDialogConfig>();
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
        setSelectedColumn={setSelectedColumn} 
        placeholder="Add a filter column"
        exclude={param}
        />
        </div>
    )
})

const ForeignRows = () => {
    const [filter, setFilter] = useState("");
    const [max, setMax] = useState(10);
    const [debouncedFilter] = useDebounce(filter, 300);
    const rlink = useRowsAsColumnsLinks();
    const fcols = useHighlightedForeignRowsAsColumns(max, debouncedFilter);
    if (!rlink) return null;
    const { linkedDs, link } = rlink;
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
}

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

const SelectionDialogComponent = () => {
    //!! this component doesn't update with HMR and introducing another wrapper component makes things worse
    //(currently changes here aren't reflected in the browser, but the rest of the components are
    //if we wrap this, then any change causes whole page to reload)
    const cols = useParamColumnsExperimental();
    useResetButton();
    return (
        <div className="p-3 absolute w-[100%] h-[100%] overflow-auto">
            {cols.map((col) => <AbstractComponent key={col.field} column={col} />)}
            <AddRowComponent />
            <ForeignRows />
        </div>
    );
};
export default SelectionDialogComponent;