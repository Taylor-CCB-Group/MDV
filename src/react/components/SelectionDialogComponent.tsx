import { useCloseOnIntersection, useConfig, useDimensionFilter, useOrderedParamColumns, usePasteHandler } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Accordion, AccordionDetails, AccordionSummary, Autocomplete, Box, Button, Checkbox, Chip, Divider, IconButton, Paper, type PaperProps, TextField, Typography } from "@mui/material";
import { type AutocompleteRenderGetTagProps, createFilterOptions } from '@mui/material/Autocomplete';
import { type MouseEvent, useCallback, useEffect, useState, useMemo, useRef, useId } from "react";

import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';
import CachedIcon from '@mui/icons-material/Cached';
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import HighlightOffIcon from '@mui/icons-material/HighlightOff';
import DragIndicatorIcon from '@mui/icons-material/DragIndicator';
import type { SelectionDialogConfig, CategoryFilter, MultiTextFilter, UniqueFilter, RangeFilter } from "./SelectionDialogReact";
import { observer } from "mobx-react-lite";
import { action, runInAction } from "mobx";
import { useChart, useDataStore } from "../context";
import ColumnSelectionComponent from "./ColumnSelectionComponent";
import type RangeDimension from "@/datastore/RangeDimension";
import { useDebounce } from "use-debounce";
import { useHighlightedForeignRowsAsColumns, useRowsAsColumnsLinks } from "../chartLinkHooks";
import * as d3 from 'd3';
import { ErrorBoundary } from "react-error-boundary";
import DebugErrorComponent from "@/charts/dialogs/DebugErrorComponent";
import { TextFieldExtended } from "./TextFieldExtended";
import { isArray } from "@/lib/utils";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";
import {
    DndContext,
    closestCenter,
    PointerSensor,
    useSensor,
    useSensors,
    type DragEndEvent,
} from '@dnd-kit/core';
import {
    arrayMove,
    SortableContext,
    verticalListSortingStrategy,
} from '@dnd-kit/sortable';
import {
    useSortable,
} from '@dnd-kit/sortable';
import { CSS } from '@dnd-kit/utilities';
import { RowsAsColsQuery } from "@/links/link_utils";
import { AUTOCOMPLETE_OPTIONS_LIMIT, AUTOCOMPLETE_TAGS_LIMIT } from "@/lib/constants";



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

const TextComponent = observer(({ column }: Props<CategoricalDataType>) => {
    const [open, setOpen] = useState(false);
    const [selectAll, setSelectAll] = useState(false);
    const ref = useRef<HTMLInputElement>(null);
    const dim = useDimensionFilter(column);
    const conf = useConfig<SelectionDialogConfig>();
    const { values } = column;
    const filter = useFilterConfig(column);
    const value = filter?.category || [];
    const setValue = useCallback((newValue: string[]) => {
        const newFilter = filter || { category: [] };
        action(() => {
            newFilter.category = newValue;
            conf.filters[column.field] = newFilter;
        })();
    }, [conf.filters, column.field, filter]);

    useCloseOnIntersection(ref, () => setOpen(false));

    useEffect(() =>{
        // todo consider having ui for this state per-dimension rather than globally
        dim.setNoClear(conf.noClearFilters);
    },[dim, conf.noClearFilters]);
 
    // react to changes in value and update the filter
    useEffect(() => {
        
        // filter could be undefined - but we previously checked for null causing component to crash
        if ((!filter) || (filter.category.length === 0)) {
            dim.removeFilter();
            return;
        }
        //probably there is better place to set this
        
        dim.filter("filterCategories", [column.field], value, true);
    }, [dim, column.field, value, filter]);

    const toggleOption = useCallback((option: string) => {
        if (value.includes(option)) {
            setValue(value.filter((v) => v !== option));
            return;
        }
        setValue([...value, option]);
    }, [setValue, value]);

    const handleSelectAll = useCallback(() => {
        if (selectAll) {
            setValue([]);
            setSelectAll(false);
        } else {
            setValue(values);
            setSelectAll(true);
        }
    }, [selectAll, setValue, values]);

    const handleValueChange = useCallback((newValue: string | string[] | null) => {
        // If newValue is null, set value to empty array
        if (!newValue) {
            setValue([]);
        } else {
            const valueArray = Array.isArray(newValue) ? newValue : [newValue];
            setValue(valueArray);
        }
    }, [setValue]);

    const handlePaste = usePasteHandler({
        options: values,
        multiple: true,
        currentValue: value,
        setValue: handleValueChange,
        getLabel: (option: string) => option,
        getValue: (option: string) => option,
    });

    return (
        <Autocomplete
            multiple
            size="small"
            options={useMemo(() => {
                const selectedSet = new Set(value);
                const selected = values.filter((v) => selectedSet.has(v));
                const unselected = values.filter((v) => !selectedSet.has(v));
                return [...selected, ...unselected];
              }, [values, value])}
            value={value}
            onChange={(_, value) => {
                if(isArray(value) && value.length === 0) {
                    setSelectAll(false);
                }
                setValue(value);
            }}
            open={open}
            onOpen={() => setOpen(true)}
            onClose={() => setOpen(false)}
            ref={ref}
            renderInput={(props) => {
                const { key, InputProps, ...p } = props as typeof props & {
                    key: string;
                }; //questionable mui types?
                return <>
                    <TextFieldExtended
                        key={key}
                        {...p}
                        placeholder={value.length === 0 ? "Type or paste to search and select": ""}
                        slotProps={{
                            input: {
                                ...InputProps,
                                onPaste: handlePaste,
                                sx: {
                                    '& .MuiInputBase-input::placeholder': {
                                    fontSize: '0.75rem'
                                    },
                                },
                            },
                        }}
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
                if (!isArray(value)) return null;

                if (value.length > AUTOCOMPLETE_TAGS_LIMIT) {
                    return (
                        <div>
                            <Typography
                                variant="button"
                                color="textSecondary"
                            >
                                {value.length} selected
                            </Typography>
                        </div>
                    )
                }
                
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

                //todo improve styling - bad from UX perspective at the moment when it overflows.
                return <div className="max-h-32 overflow-auto pr-3">{chips}</div>;
            }}
            // Getting warnings in console if slotProps is used for rendering Paper Component
            PaperComponent={useMemo(() => (paperProps: PaperProps) => {
                const { children, ...restPaperProps } = paperProps;
                const showMessage = values.length > AUTOCOMPLETE_OPTIONS_LIMIT;
                return (
                  <Paper {...restPaperProps}>
                    <Box
                        onMouseDown={(e) => e.preventDefault()}
                        py={1}
                        px={2}
                        sx={{textAlign: "center"}}
                    >
                        <Button onClick={handleSelectAll} color="inherit" fullWidth>
                            {selectAll ? "Unselect All" : "Select All"}
                        </Button>
                    </Box>
                    <Divider />
                    {showMessage && (
                        <>
                            <Box 
                                py={1}
                                px={2}
                                sx={{
                                    textAlign: "center",
                                    fontSize: "0.8rem",
                                    color: "text.secondary",
                                    backgroundColor: "var(--menu_bar_color)",
                                    border: "2px solid var(--fade_background_color)"
                                }}
                            >
                                Showing first {AUTOCOMPLETE_OPTIONS_LIMIT} of {values.length} options.
                            </Box>
                            <Divider />
                        </>

                    )}
                    {children}
                  </Paper>
                );
              }, [selectAll, handleSelectAll, values])}
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
    const conf = useConfig<SelectionDialogConfig>()
    //nb - we may want to allow value to be null, rather than defaulting to minMax
    //relates to e.g. clearBrush function
    // const value = (filters[column.field] || column.minMax) as [number, number];
    const fVal = conf.filters[column.field] as RangeFilter | null;
    // if (!isArray(fVal)) throw new Error("Expected range filter to be an array");
    const value = fVal;
    // const value = fVal;
    const isInteger = column.datatype.match(/int/);
    const { minMax } = column;
    const step = useMemo(() => {
        if (isInteger) return 1;
        // not sure this is totally correct - but there was a problem with very small ranges
        // this should be better...
        const small = Math.abs(minMax[1] - minMax[0]);
        return small < 0.001 ? small/1000 : 0.001;
    }, [isInteger, minMax]);
    const [debouncedValue] = useDebounce(value, 10);
     useEffect(() =>{
        filter.setNoClear(conf.noClearFilters);
    },[filter, conf.noClearFilters]);
    // Effect to manage the filter state
    useEffect(() => {
        const value = debouncedValue;
        if (!value) {
            filter.removeFilter();
            return;
        }
        const [min, max] = value;
        filter.filter("filterRange", [column.field], { min, max }, true);
    }, [column, filter, debouncedValue]);

    const [histogram, setHistogram] = useState<number[]>([]);
    // this could be a more general utility function - expect to extract soon
    const queryHistogram = useCallback(async () => {
        // waste of life trying to use Dimension class.
        // filter.getBinsAsync(column.name, { bins: 100 }).then((histogram) => {
        //     setHistogram(histogram);
        // });
        const worker = new Worker(new URL("../../datastore/rawHistogramWorker.ts", import.meta.url));
        worker.onmessage = (event) => {
            setHistogram(event.data);
            worker.terminate();
        };
        const isInt32 = column.datatype === "int32";
        const originalData = column.data as Float32Array | Int32Array;
        const data = new SharedArrayBuffer(originalData.length * 4);
        new Float32Array(data).set(originalData);
        worker.postMessage({ data, min: column.minMax[0], max: column.minMax[1], bins: 100, isFloat: isInt32 });
    }, [column]);

    const [low, high] = value || column.minMax;
    const [min, max] = column.minMax;
    const lowFraction = (low - min) / (max - min);
    const highFraction = (high - min) / (max - min);

    return { value, step, histogram, lowFraction, highFraction, queryHistogram };
}
// type set2d = ReturnType<typeof useState<[number, number]>>[1];
type Range = [number, number];
type set2d = (v: Range | null) => void; //nb, setting undefined can actually be problematic
type RangeProps = ReturnType<typeof useRangeFilter> & {
    setValue: set2d,
    minMax: Range,
    // probably want to review how these are specified / controlled
    histoWidth: number, //number of bins
    histoHeight: number, //height of the histogram
};

const useBrushX = (
    ref: React.RefObject<SVGSVGElement>,
    { value, setValue, minMax, histoWidth, histoHeight }: RangeProps //consider different typing here
) => {
    const brushRef = useRef<(ReturnType<typeof d3.brushX>) | null>(null);
    // we need to be able to respond to changes in value - but without causing an infinite loop
    // or having the brush reset on every render
    const [initialValue] = useState(value);

    useEffect(() => {
        if (!ref.current) return;

        const svg = d3.select(ref.current);
        // Set up brush
        const brush = d3.brushX()
            .handleSize(1)
            .extent([[0, -2], [histoWidth, histoHeight+2]])
            .on("brush end", (event) => {
                if (event.selection) {
                    const [start, end] = event.selection.map((x: number) => {
                        if (!ref.current) {
                            console.error("No ref.current in brush event handler");
                            return 0;
                        }
                        const { width } = ref.current.getBoundingClientRect();
                        // Normalize x-coordinate to [minMax[0], minMax[1]]
                        const r = width / histoWidth;
                        const normalizedX = r * x / width;
                        return minMax[0] + normalizedX * (minMax[1] - minMax[0]);
                    });
                    setValue([start, end]);
                } else {
                    // warning - the null value here does behave distinctly differently from undefined
                    // e.g. as of this writing, the reset button will be glitchy if we don't use null here
                    setValue(null); // null - reset to full range if brush is cleared
                }
            });

        brushRef.current = brush;

        // Apply the brush to the SVG
        const brushGroup = svg.append("g").attr("class", "brush").call(brush);
        // Initialize brush selection based on the initial value
        if (initialValue) {
            const [start, end] = initialValue.map(
                (v) => ((v - minMax[0]) / (minMax[1] - minMax[0])) * histoWidth
            );
            brushGroup.call(brush.move, [start, end]); // Move the brush to the initial selection
        }

        // Apply `vectorEffect` directly to handles
        brushGroup.selectAll(".selection").attr("vector-effect", "non-scaling-stroke");
        brushGroup.selectAll(".handle").attr("vector-effect", "non-scaling-stroke");

        // Cleanup on unmount
        return () => {
            svg.select(".brush").remove();
        };
    }, [ref, setValue, minMax, histoWidth, histoHeight, initialValue]);

    const [debouncedValue] = useDebounce(value, 100, { 
        equalityFn: (a, b) => {
            //although the type of input argument is [number, number] | null - they are undefined when component is unmounted
            //! which causes an exception here which breaks the whole chart
            //so rather than checking === null, we check for falsy values
            if (!a && !b) return true;
            if (!a || !b) return false;
            return a[0] === b[0] && a[1] === b[1];
        }
    });
    const setBrushValue = useCallback<set2d>((v) => {
        if (!brushRef.current || !ref.current) return;
        const svg = d3.select(ref.current);

        if (!v) {
            // throw new Error("this is actually ok, but I want to test the error handling");
            //@ts-ignore life is too short
            svg.select(".brush").call(brushRef.current.move, null);
            return;
        }
        const [start, end] = v;
        const x0 = (start - minMax[0]) / (minMax[1] - minMax[0]) * histoWidth;
        const x1 = (end - minMax[0]) / (minMax[1] - minMax[0]) * histoWidth;
        //@ts-ignore life is too short
        svg.select(".brush").call(brushRef.current.move, [x0, x1]);
    }, [minMax, histoWidth, ref]);//why doesn't biome think we need brushRef?
    useEffect(() => {
        setBrushValue(debouncedValue);
    }, [debouncedValue, setBrushValue]);
};
const Histogram = observer((props: RangeProps) => {
    const { histogram: data, queryHistogram, value } = props;
    const { histoWidth, histoHeight } = props;
    const ref = useRef<SVGSVGElement>(null);
    useBrushX(ref, props);
    const prefersDarkMode = window.mdv.chartManager.theme === "dark";
    const width = histoWidth;
    const height = histoHeight;
    const lineColor = prefersDarkMode ? '#fff' : '#000';
    // Find max value for vertical scaling
    const maxValue = Math.max(...data);

    // Define the padding and scaling factor
    const padding = 2;
    const xStep = data.length / (width + 1); // Space between points
    const yScale = (height - 2 * padding) / maxValue; // Scale based on max value

    const [hasQueried, setHasQueried] = useState(false);
    useEffect(() => {
        if (!ref.current) return;
        const observer = new IntersectionObserver((entries) => {
            if (entries[0].isIntersecting && !hasQueried) {
                setHasQueried(true);
                queryHistogram();
            }
        }, { rootMargin: '0px 0px 100px 0px' });
        observer.observe(ref.current);
        // queryHistogram();
        return () => observer.disconnect();
    }, [queryHistogram, hasQueried]);

    // Generate the points for the polyline
    // ??? useMemo was wrong ????
    const points = useMemo(() => data.map((value, index) => {
        const x = index * xStep;
        const y = height - padding - value * yScale;
        return `${x},${y}`;
    }).join(' '), [data, xStep, yScale, height]);
    const v = value || props.minMax;
    return (
        <>
        <svg width={'100%'} height={height}
        viewBox={`0 0 ${width} ${height}`}
        preserveAspectRatio="none"
        ref={ref}
        cursor="move"
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
            {/* d3.brushX will add more elements as a side-effect, handled in hook */}
        </svg>
        {/* <p className="flex justify-between"><em>{`${v[0].toFixed(2)}<`}</em> <em>{`<${v[1].toFixed(2)}`}</em></p> */}
        </>
    );
});

const NumberComponent = observer(({ column }: Props<NumberDataType>) => {
    const filters = useConfig<SelectionDialogConfig>().filters;
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
            <Histogram {...rangeProps} setValue={setValue} minMax={column.minMax} histoWidth={99} histoHeight={100} />
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
    const {
        attributes,
        listeners,
        setNodeRef,
        transform,
        transition,
        isDragging,
    } = useSortable({ id: column.field });

    const style = {
        transform: CSS.Translate.toString(transform),
        transition,
        opacity: isDragging ? 0.5 : 1,
    };

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
            if (!isArray(config.param)) throw new Error("expected param array");
            
            // Filter the config.param by removing the column
            const newParam = config.param.filter(p => {
                if (typeof p === "string") {
                    return p !== column.field;
                }
                if (p instanceof RowsAsColsQuery) {
                    return false;
                }
                return true;
            });
            
            // Update the config.param with new array
            config.param = newParam;
            
            // Delete the column field in config.filters
            delete config.filters[column.field];
            if (config.order) {
                // Delete the column field in config.order
                delete config.order[column.field];
            }
        });
    }, [config, column]);
    return (
        <Accordion defaultExpanded={true}
            onMouseEnter={() => setIsHovered(true)}
            onMouseLeave={() => setIsHovered(false)}
            style={style}
            ref={setNodeRef}
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
                <div className="flex items-center h-4 w-full">
                    <IconButton
                        {...attributes}
                        {...listeners}
                        size="small"
                        sx={{
                            opacity: isHovered ? 0.7 : 0.3,
                            transition: 'opacity 0.3s',
                            cursor: 'grab',
                            '&:active': {
                                cursor: 'grabbing',
                            },
                            '&:hover': {
                                opacity: 1,
                            },
                        }}
                    >
                        <DragIndicatorIcon fontSize="small" />
                    </IconButton>
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
                    ({ error }) => <DebugErrorComponent error={error} title="Unexpected Error: please report to developers." />
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
});

type RLink = ReturnType<typeof useRowsAsColumnsLinks>[0];

const LinkComponent = observer(({ rlink, linkIndex }: { rlink: RLink, linkIndex: number }) => {
    const [filter, setFilter] = useState("");
    const [max, setMax] = useState(10);
    const [debouncedFilter] = useDebounce(filter, 300);
    // if we re-instate this, perhaps temporarily, it could be a useful way to test showing multiple links and subgroups.
    const fcols = useHighlightedForeignRowsAsColumns(max, debouncedFilter, linkIndex);
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
});

const ForeignRows = observer(() => {
    const rlink = useRowsAsColumnsLinks();
    return (
        <>
        {rlink.map((link, i) => <LinkComponent key={link.link.name} rlink={link} linkIndex={i} />)}
        </>
    )
});

/**
 * This will control the behaviour of the reset menuIcon in the chart header - not rendered with react.
 */
function useResetButton() {
    const chart = useChart();
    const conf = useConfig<SelectionDialogConfig>();
    // todo: what about category filters with empty array?
    const hasFilter = Object.values(conf.filters).some((f) => f !== null);
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
            //actually removeAllFilters already clears all the filters (if they are 
            //not tagged with noclear in single action for performance
            //filters will be cleared again here but it is not much of an issue
            if (!conf.noClearFilters){
                runInAction(() => {
                    for (const key in conf.filters) {
                        delete conf.filters[key];
                    }
                });
            }
        };
        ds.addListener(k, resetAll);
        return () => ds.removeListener(k);
    }, [ds, id, conf.filters]);
}

const SelectionDialogComponent = () => {
    //!! this component doesn't update with HMR and introducing another wrapper component makes things worse
    //(currently changes here aren't reflected in the browser, but the rest of the components are
    //if we wrap this, then any change causes whole page to reload)
    //! since this is outside the ErrorBoundary, problems in this hook are not well handled
    //we could consider returning some kind of `Result` object from this hook...
    const orderedParams = useOrderedParamColumns();
    const config = useConfig<SelectionDialogConfig>();
    useResetButton();
    const showAddRow = false;

    const sensors = useSensors(
        useSensor(PointerSensor, {
            activationConstraint: {
                // Need 8px of movement to start dragging
                distance: 8,
            }
        })
    );

    const onDragEnd = useCallback((event: DragEndEvent) => {
        // active - dragged item, over - dropped item
        const { active, over } = event;

        // Check if dragged and dropped item are different and other sanity checks
        if (active && over && active.id !== over?.id) {
            const oldIndex = orderedParams.findIndex(p => p.field === active.id);
            const newIndex = orderedParams.findIndex(p => p.field === over.id);

            if (oldIndex !== -1 && newIndex !== -1) {
                // Update the config.order with the new order
                const newOrderedFields = arrayMove(orderedParams.map(p => p.field), oldIndex, newIndex);
                const newOrder: Record<string, number> = {};
                newOrderedFields.forEach((field, index) => {
                    newOrder[field] = index;
                });
                
                runInAction(() => {
                    config.order = newOrder;
                });
            }
        }
    }, [config, orderedParams]);
    return (
        <div className="p-3 absolute w-[100%] h-[100%] overflow-x-hidden overflow-y-auto">
            <DndContext
                sensors={sensors}
                collisionDetection={closestCenter}
                onDragEnd={onDragEnd}
            >
                <SortableContext 
                    items={orderedParams.map(col => col.field)}
                    strategy={verticalListSortingStrategy}
                >
                    {orderedParams.map((col) => <AbstractComponent key={col.field} column={col} />)}
                </SortableContext>
            </DndContext>
            {showAddRow && <ErrorBoundary FallbackComponent={
                ({ error }) => 
                    (
                        <div className="col-span-3 w-full"> 
                            <ErrorComponentReactWrapper 
                                error={{message: error.message, stack: error.stack}} 
                                //todo assign proper meta data
                                // extraMetaData={} 
                                title="Error displaying 'AddRowComponent'. Click to view details"
                            />
                        </div>
                        )
            }>
                <AddRowComponent />
                <ForeignRows />
            </ErrorBoundary>}
        </div>
    );
};
export default SelectionDialogComponent;
