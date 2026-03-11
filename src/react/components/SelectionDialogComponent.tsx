import { useCloseOnIntersection, useConfig, useDimensionFilter, useOrderedParamColumns, usePasteHandler, useSimplerFilteredIndices } from "../hooks";
import type { CategoricalDataType, NumberDataType, DataColumn, DataType } from "../../charts/charts";
import { Accordion, AccordionDetails, AccordionSummary, Autocomplete, Box, Button, Checkbox, Chip, Divider, IconButton, Paper, type PaperProps, TextField, Typography } from "@mui/material";
import { type MouseEvent, useCallback, useEffect, useState, useMemo, useId, useRef } from "react";
import { useQuery } from "@tanstack/react-query";

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
import { useHighlightedForeignRowsAsColumns, useRowsAsColumnsLinks } from "../chartLinkHooks";
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
import { useHighlightedIndices } from "../selectionHooks";
import HistogramWidget, {
    type HistogramLayer,
    type HistogramScaleType,
} from "./HistogramWidget";
import {
    getNumericColumnData,
    getSharedNumericColumnData,
} from "@/lib/columnTypeHelpers";
import { createHistogram, queryHistogramWorker } from "@/react/utils/histogram";
import { useDebounce } from "use-debounce";



const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;
const HISTOGRAM_BINS = 100;

type Props<K extends DataType> = {
    column: DataColumn<K>;
    // filter not passed in as a possibly null prop...
    // each component will extract it from the config.filters object, and write it back there.
    // mobx will handle the reactivity, and it should get serialized appropriately.
    // filter?: K extends CategoricalDataType ? CategoryFilter : [number, number];
};

/**
 * Retrieve the current filter configuration for a given data column.
 *
 * @param column - The data column whose filter configuration to read
 * @returns The filter object associated with `column.field` for the column's datatype (`CategoryFilter`, `RangeFilter`, `MultiTextFilter`, or `UniqueFilter`), or `null` if no filter is set
 */
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
    const filter = useDimensionFilter(column);
    const conf = useConfig<SelectionDialogConfig>()
    //nb - we may want to allow value to be null, rather than defaulting to minMax
    //relates to e.g. clearBrush function
    // const value = (filters[column.field] || column.minMax) as [number, number];
    const fVal = useFilterConfig(column);
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

    const filteredIndices = useSimplerFilteredIndices();
    const highlightedIndices = useHighlightedIndices();
    const numericData = getNumericColumnData(column);
    const overallHistogramQuery = useQuery<number[]>({
        queryKey: [
            "selection-dialog-overall-histogram",
            column.field,
            column.datatype,
            column.minMax[0],
            column.minMax[1],
            numericData.length,
        ],
        queryFn: ({ signal }) => {
            const { data, arrayType, byteOffset, length } = getSharedNumericColumnData(numericData);
            return queryHistogramWorker({
                data,
                min: column.minMax[0],
                max: column.minMax[1],
                bins: HISTOGRAM_BINS,
                arrayType,
                byteOffset,
                length,
            }, signal);
        },
        enabled: false,
        retry: false,
        staleTime: Number.POSITIVE_INFINITY,
    });
    const overallHistogram = overallHistogramQuery.data ?? [];
    const queryHistogram = useCallback(async () => {
        await overallHistogramQuery.refetch();
    }, [overallHistogramQuery]);

    const filteredHistogram = useMemo(
        () => createHistogram(numericData, column.minMax[0], column.minMax[1], HISTOGRAM_BINS, filteredIndices),
        [numericData, column.minMax, filteredIndices],
    );
    const highlightedHistogram = useMemo(
        () => createHistogram(numericData, column.minMax[0], column.minMax[1], HISTOGRAM_BINS, highlightedIndices),
        [numericData, column.minMax, highlightedIndices],
    );

    const [low, high] = value || column.minMax;
    const [min, max] = column.minMax;
    const lowFraction = (low - min) / (max - min);
    const highFraction = (high - min) / (max - min);

    return {
        value,
        step,
        overallHistogram,
        filteredHistogram,
        highlightedHistogram,
        overallHistogramError: overallHistogramQuery.error,
        lowFraction,
        highFraction,
        queryHistogram,
    };
}
// type set2d = ReturnType<typeof useState<[number, number]>>[1];
type Range = [number, number];
type ScaleMode = "auto" | HistogramScaleType;
type set2d = (v: Range | null) => void; //nb, setting undefined can actually be problematic
type RangeProps = ReturnType<typeof useRangeFilter> & {
    setValue: set2d,
    minMax: Range,
    // probably want to review how these are specified / controlled
    histoWidth: number, //number of bins
    histoHeight: number, //height of the histogram
};
const SCALE_MODE_ORDER: ScaleMode[] = ["auto", "linear", "log"];

const nextScaleMode = (mode: ScaleMode): ScaleMode =>
    SCALE_MODE_ORDER[(SCALE_MODE_ORDER.indexOf(mode) + 1) % SCALE_MODE_ORDER.length];

const resolveAutoXScale = (domain: Range): HistogramScaleType => {
    const [min, max] = domain;
    if (!Number.isFinite(min) || !Number.isFinite(max) || min === max) {
        return "linear";
    }
    const shiftedMin = Math.min(Math.abs(min), Math.abs(max)) < 1e-9
        ? 1e-9
        : Math.max(1e-9, Math.min(Math.abs(min), Math.abs(max)));
    const shiftedMax = Math.max(Math.abs(min), Math.abs(max), shiftedMin);
    return shiftedMax / shiftedMin > 500 ? "log" : "linear";
};

const resolveAutoYScale = (histogram: number[]): HistogramScaleType => {
    const nonZero = histogram.filter((value) => value > 0);
    if (nonZero.length < 2) return "linear";
    const min = Math.min(...nonZero);
    const max = Math.max(...nonZero);
    const mean = nonZero.reduce((sum, value) => sum + value, 0) / nonZero.length;
    return max / min > 50 || max / Math.max(1, mean) > 10 ? "log" : "linear";
};

const Histogram = observer((props: RangeProps) => {
    const { overallHistogram, filteredHistogram, highlightedHistogram, overallHistogramError, queryHistogram } = props;
    const { histoWidth, histoHeight } = props;
    const [xScaleMode, setXScaleMode] = useState<ScaleMode>("auto");
    const [yScaleMode, setYScaleMode] = useState<ScaleMode>("auto");
    const prefersDarkMode = window.mdv.chartManager.theme === "dark";
    const overallColor = prefersDarkMode ? "rgba(255,255,255,0.35)" : "rgba(15,23,42,0.22)";
    const filteredColor = prefersDarkMode ? "rgba(96,165,250,0.8)" : "rgba(37,99,235,0.78)";
    const highlightColor = prefersDarkMode ? "rgba(251,191,36,0.95)" : "rgba(217,119,6,0.95)";
    const emptyHistogram = useMemo(() => new Array(HISTOGRAM_BINS).fill(0), []);
    const backgroundData = overallHistogram.length > 0 ? overallHistogram : emptyHistogram;
    const filteredData = filteredHistogram.length > 0 ? filteredHistogram : emptyHistogram;
    const highlightedData = highlightedHistogram.length > 0 ? highlightedHistogram : emptyHistogram;
    const resolvedXScale = useMemo(
        () => (xScaleMode === "auto" ? resolveAutoXScale(props.minMax) : xScaleMode),
        [props.minMax, xScaleMode],
    );
    const resolvedYScale = useMemo(
        () => (yScaleMode === "auto" ? resolveAutoYScale(backgroundData) : yScaleMode),
        [backgroundData, yScaleMode],
    );
    useEffect(() => {
        if (!overallHistogramError) return;
        console.error("Failed to query overall histogram", overallHistogramError);
    }, [overallHistogramError]);
    const layers = useMemo<HistogramLayer[]>(() => [
        {
            id: "overall",
            data: backgroundData,
            color: overallColor,
            variant: "bars",
            widthFactor: 0.8,
        },
        {
            id: "filtered",
            data: filteredData,
            color: filteredColor,
            variant: "bars",
            inset: 0.15,
            widthFactor: 0.7,
            radius: 0.4,
        },
        {
            id: "highlighted",
            data: highlightedData,
            color: highlightColor,
            variant: "markers",
            widthFactor: 0.35,
        },
    ], [backgroundData, filteredData, highlightedData, overallColor, filteredColor, highlightColor]);
    const brush = useMemo(() => ({
        value: props.value,
        setValue: props.setValue,
        minMax: props.minMax,
    }), [props.value, props.setValue, props.minMax]);
    const handleVisibleOnce = useCallback(() => {
        void queryHistogram();
    }, [queryHistogram]);
    return (
        <>
        <div className="mb-1 flex items-center justify-end gap-1 text-[10px] opacity-75">
            <button
                type="button"
                className="rounded border px-1.5 py-0.5"
                onClick={() => setXScaleMode((mode) => nextScaleMode(mode))}
            >
                X:{xScaleMode === "auto" ? resolvedXScale : xScaleMode}
            </button>
            <button
                type="button"
                className="rounded border px-1.5 py-0.5"
                onClick={() => setYScaleMode((mode) => nextScaleMode(mode))}
            >
                Y:{yScaleMode === "auto" ? resolvedYScale : yScaleMode}
            </button>
        </div>
        <HistogramWidget
            layers={layers}
            width={histoWidth}
            height={histoHeight}
            bins={HISTOGRAM_BINS}
            xScaleType={resolvedXScale}
            yScaleType={resolvedYScale}
            brush={brush}
            onVisibleOnce={handleVisibleOnce}
        />
        <div className="mt-1 flex flex-wrap gap-x-3 gap-y-1 text-[11px] opacity-75">
            <span className="inline-flex items-center gap-1">
                <span className="inline-block h-2 w-2 rounded-sm" style={{ backgroundColor: overallColor }} />
                overall
            </span>
            <span className="inline-flex items-center gap-1">
                <span className="inline-block h-2 w-2 rounded-sm" style={{ backgroundColor: filteredColor }} />
                filtered
            </span>
            <span className="inline-flex items-center gap-1">
                <span className="inline-block h-2 w-2 rounded-full" style={{ backgroundColor: highlightColor }} />
                highlighted
            </span>
        </div>
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
    }, [ds, id, conf.filters, conf.noClearFilters]);
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
