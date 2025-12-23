import { observer } from "mobx-react-lite";
import { useMemo, useId, useCallback, useState, useRef, useEffect } from "react";
import type { AnyGuiSpec, DropDownValues, GuiSpec, GuiSpecType, Disposer } from "../../charts/charts";
import { action, makeAutoObservable } from "mobx";
import { ErrorBoundary } from "react-error-boundary";
import {
    Accordion,
    AccordionContent,
    AccordionItem,
    AccordionTrigger,
} from "@/components/ui/accordion";
import { v4 as uuid } from "uuid";
import { Button, Chip, FormControl, FormControlLabel, Radio, RadioGroup, Slider, Typography } from "@mui/material";
import Checkbox from "@mui/material/Checkbox";
import TextField from "@mui/material/TextField";
import Autocomplete from "@mui/material/Autocomplete";
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import { ChartProvider } from "../context";
import ColumnSelectionComponent from "./ColumnSelectionComponent";
import { inferGenericColumnSelectionProps } from "@/lib/columnTypeHelpers";
import { isArray, matchString, notEmpty, parseDelimitedString } from "@/lib/utils";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";
import { useCloseOnIntersection, usePasteHandler } from "../hooks";

export const MLabel = observer(({ props, htmlFor }: { props: AnyGuiSpec, htmlFor?: string }) => (
    <Typography fontSize="small" sx={{alignSelf: "center", justifySelf: "end", textAlign: "right", paddingRight: 2}}>
        {props.label} 
        {/* <em className="opacity-30"> ({props.type})</em> */}
    </Typography>
    // todo fix justifySelf - it's not working as expected
    // <FormControlLabel
    //     sx={{ justifySelf: "right" }}
    //     control={<Typography>{props.label}</Typography>}
    //     label=""
    //     htmlFor={htmlFor}
    // />
));

const TextComponent = ({ props }: { props: GuiSpec<"text"> }) => (
    <>
        <MLabel props={props} />
        <TextField
            value={props.current_value}
            onChange={action((e) => {
                props.current_value = e.target.value;
                if (props.func) props.func(e.target.value);
            })}
            className="w-full"
        />
    </>
);

const TextBoxComponent = ({ props }: { props: GuiSpec<"textbox"> }) => (
    <>
        <MLabel props={props} />
        <div />
        <textarea
            value={props.current_value}
            onChange={action((e) => {
                props.current_value = e.target.value;
                if (props.func) props.func(e.target.value);
            })}
            className="w-full col-span-2"
        />
    </>
);

const SliderComponent = ({ props }: { props: GuiSpec<"slider"> }) => (
    <>
        <MLabel props={props} />
        <Slider
            value={props.current_value}
            min={props.min || 0}
            max={props.max || 1}
            step={props.step || 0.01} //todo - infer default step from min/max
            onChange={action((_, value) => {
                // const value = Number.parseFloat(e.target.value);
                if (Array.isArray(value))
                    throw new Error("slider callback should have multiple values");
                props.current_value = value;
                if (props.func) props.func(value);
            })}
            size="small"
        />
    </>
);

const SpinnerComponent = ({ props }: { props: GuiSpec<"spinner"> }) => (
    <>
        <MLabel props={props} />
        <input
            className="p-1"
            type="number"
            value={props.current_value}
            min={props.min || 0}
            max={props.max}
            step={props.step || 1}
            onChange={action((e) => {
                const value = (props.current_value = Number.parseInt(
                    e.target.value,
                ));
                if (props.func) props.func(value);
            })}
        />
    </>
);
export type ColumnSelectionSpec = GuiSpec<"column"> | GuiSpec<"multicolumn">;

/**
 * 
 * 
 * The situation with trying to marshall types of column widgets, for settings vs addChartDialog,
 * is making this all a bit of a mess. There should be a simpler way to express this stuff - which
 * may reviewing the ways widgets are defined in different contexts.
 */
function isMultiSpec(props: ColumnSelectionSpec): props is GuiSpec<"multicolumn"> {
    return props.type === "multicolumn";
}


/**
 * Wrap the ColumnSelectionComponent in a setting GUI component.
 *
 * nb, for some weird reason if this is defined in ColumnSelectionComponent.tsx HMR doesn't work...
 */
export const ColumnSelectionSettingGui = observer(({ props }: { props: ColumnSelectionSpec }) => {
    // multiple etc... this is the level above ColumnSelectionComponent I need to be focusing on.
   const multiple = isMultiSpec(props);

   const filter = props.columnType;
    // not only is the filter not working, but we need to decide how to express "multiple"
    const props2 = useMemo(() => inferGenericColumnSelectionProps({
        // fixing this stuff is high priority
        //@ts-expect-error ColumnSelection setSelectedColumn(v: never)???
        setSelectedColumn: action((v) => {
            props.current_value = v;
            props.func?.(v);
        }),
        type: filter, //is it ok for this to be optional?
        multiple,
        current_value: props.current_value //! ideally this would also react to changes in the mobx store, why doesn't it?
        // current_value: props.current_value... maybe want to be more mobx-y about this
    }), [filter, props, props.current_value, multiple]);
    return (
        <>
            <MLabel props={props} />
            <ColumnSelectionComponent {...props2} />
        </>
    );
});

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

// nb this is not the same as `GuiSpec<"dropdown" | "multidropdown">`
export type DropdownSpec = GuiSpec<"dropdown"> | GuiSpec<"multidropdown">;
// type DropdownSpec = GuiSpec<"dropdown" | "multidropdown">;
function getOptionAsObjectHelper(values: DropDownValues) {
    const useObjectKeys = values.length === 3;
    const [_valueObjectArray, labelKey, valueKey] = useObjectKeys ? values : ["X", "X", "X"];
    if ((labelKey === undefined) || (valueKey === undefined)) throw new Error("DropdownAutocompleteComponent requires labelKey and valueKey for useObjectKeys");
    type Option = { label: string, value: string, original: any };
    // might put some more elaborate types here on original at some point
    return (original: string | any): Option => {
        const isString = typeof original === "string";
        if (isString === useObjectKeys) throw new Error("DropdownAutocompleteComponent requires values to be either all strings or all objects");
        if (isString) {
            return { label: original, value: original, original };
        }
        const label: string = original[labelKey];
        const value: string = original[valueKey];
        return { label, value, original };
    };
}
function useDropdownOptions(props: DropdownSpec) {
    // this check should be redundant with current types, but leaving it in for now
    if (!props.values) throw new Error("DropdownAutocompleteComponent requires props.values - probable logic error in type definitions");
    const toOption = useCallback(getOptionAsObjectHelper(props.values), []);
    type OptionType = ReturnType<typeof toOption>;
    const NO_OPTIONS = [{ label: "No options (error - sorry...)", value: "", original: "" }];
    // useMemo here was problematic because the values array is a stable reference...
    // *mobx can be a pain* or at least it goes against the grain of react and functional programming
    //! is always recomputing this when we hit this hook a problem? should test with large datasets
    // actually seemed to be ok, but useMemo works if we're more careful about the dependencies
    const options = useMemo(() => props.values?.[0].map(toOption) || NO_OPTIONS, [props.values[0], toOption]);
    
    // bit of a faff with sometimes getting a one-item array, sometimes a single item...
    const getSingleOption = useCallback(
        (option: OptionType | OptionType[] | undefined) => {
            if (!option) return NO_OPTIONS[0];
            const a = isArray(option);
            if (a && option.length > 1) {
                console.warn("ideally we shouldn't have to deal with arrays at all here, but we only expect one value when we do");
            }
            return (isArray(option) ? option[0] : option)
        },
        []
    );
    // coerce to single option
    const single = getSingleOption;
    // get label for option
    // -- do we need `| OptionType[]` here - perhaps `getSingleOption` can go away as well?
    const label = useCallback((option: OptionType | OptionType[]) => single(option).label, [single]);
    // const label = useCallback((option: OptionType) => single(option).label, [single]);
    // get value for option
    const val = useCallback((option: OptionType | OptionType[]) => single(option).value, [single]);

    // ------
    // deal with cases where the options have changed (e.g. values from a different column
    // when a different "contour parameter" is selected by another GUI element)
    // so available values may be incompatible with props.current_value
    // ------
    // if 'multiple' is true, make sure we get an array even if one / zero values selected.
    // this should now be handled by type-predicate isMultidropdown
    // const multiple = isMultidropdown(props);
    const v = props.current_value;
    const validVal = useCallback(
        (v: string) => options.some((item) => item.value === v),
        [options],
    );
    const isVArray = isArray(v);
    const allValid = isVArray ? v.every(validVal) : validVal(v);
    const okValue = allValid ? v : isVArray ? v.filter(validVal) : null;
    //map from 'value' string to option object
    const okOption = (isArray(okValue)
        ? okValue.map((v) => options.find((o) => o.value === v))
        : [options.find((o) => o.value === v)]).filter(notEmpty);
    // : options.find((o) => o.value === v); //not-multiple...

    return {options, single, label, val, okOption};
}
function isMultidropdown(props: DropdownSpec): props is GuiSpec<"multidropdown"> {
    const multi = props.type === "multidropdown";
    if (multi !== isArray(props.current_value)) throw new Error("current_value should be an array if (and only if) type is multidropdown");
    return multi;
}
function inferDropdownType<T extends "dropdown" | "multidropdown">(props: GuiSpec<T>): GuiSpec<T> {
    return props;
}
function getCurrentValue(props: DropdownSpec) {
    const multi = isMultidropdown(props);
    if (multi) return props.current_value; //ts correctly infers this is an array
    // return props.current_value; //.. fails to infer this is a string
    if (isArray(props.current_value)) throw new Error("current_value should be an array if and only if type is multidropdown");
    return props.current_value;
}
export const DropdownAutocompleteComponent = observer(({
    props//: propsAmbiguous,
}: { props: DropdownSpec }) => {
    // todo review 'virtualization' for large lists
    const id = useId();
    const [open, setOpen] = useState(false);
    const ref = useRef<HTMLInputElement>(null);
    // const props = inferDropdownType(propsAmbiguous);
    const multiple = isMultidropdown(props);
    // the props.values may be a tuple of [valueObjectArray, textKey, valueKey], or an array of length 1 - [string[]]
    if (!props.values) throw new Error("DropdownAutocompleteComponent requires props.values");
    //todo handle multitext / tags properly.

    const {options, single, label, val, okOption} = useDropdownOptions(props);

    useCloseOnIntersection(ref, () => setOpen(false));

    type Option = typeof options[number]; // hopefully we can improve `{ original: any }`
    // still not entirely sure about the type for onChange...
    type OVal = Option | Option[] | (Option | Option[])[] | null;

    const handleValueChange = useCallback((newValue: Option | Option[] | null) => {
        if (multiple) {
            const valueArray = Array.isArray(newValue) ? newValue : (newValue ? [newValue] : []);
            const valueStrings = valueArray.map(val);
            props.current_value = valueStrings;
            props.func?.(valueStrings);
        } else {
            const valueSingle = Array.isArray(newValue) ? newValue[0] : newValue;
            if (valueSingle) {
                const valueString = val(valueSingle);
                props.current_value = valueString;
                props.func?.(valueString);
            }
        }
    }, [multiple, props, val]);

    const handlePaste = usePasteHandler({
        options,
        currentValue: okOption.filter((option): option is Option => option !== undefined),
        setValue: handleValueChange,
        multiple,
        getValue: (option: Option) => option,
        getLabel: (option: Option) => label(option),
    });

    return (
        <>
            <MLabel htmlFor={id} props={props} />
            <Autocomplete
                className="w-full"
                multiple={multiple}
                size="small"
                id={id}
                options={options}
                limitTags={5}
                disableCloseOnSelect={multiple}
                getOptionLabel={label}
                value={multiple ? 
                    okOption.filter((a) => a !== undefined) 
                    : okOption.length > 0 ? okOption : null
                }
                open={open}
                onOpen={() => setOpen(true)}
                onClose={() => setOpen(false)}
                ref={ref}
                onChange={action((_, value: OVal) => {
                    //added type annotation above because mobx seems to fluff the inference to `never`
                    if (value === null) return;
                    if (isArray(value)) {
                        if (!multiple) throw new Error("shouldn't get an array here");
                        const vals = value.map(val);
                        props.current_value = vals;
                        props.func?.(vals);
                        return;
                    }
                    if (multiple) throw new Error("shouldn't get a single value here");
                    const v = val(value);
                    props.current_value = v;
                    if (props.func) props.func(v);
                })}
                isOptionEqualToValue={(option, value) => (single(option).original === single(value).original)}
                renderOption={(props, option, { selected }) => {
                    if (!option) return null; //FFS let's switch to Rust or something
                    // we could potentially render something richer here - like color swatches or icons
                    // if we had a richer sense of the data in context
                    // ^^ e.g. if we had a `GuiSpec<'category'>` it could have it's own internal logic for
                    // managing internal state, and also use a richer component here.
                    const { key, ...optionProps } = props as typeof props & {
                        key: string;
                    }; //questionable mui types?
                    if (multiple)
                        return (
                            <li key={key} {...optionProps}>
                                <Checkbox
                                    icon={icon}
                                    checkedIcon={checkedIcon}
                                    style={{ marginRight: 8 }}
                                    checked={selected}
                                />
                                {label(option)}
                            </li>
                        );
                    return (
                        <li key={key} {...optionProps}>
                            {label(option)}
                        </li>
                    );
                }}
                renderTags={(tagValue, getTagProps) => {
                    //seems to be a material-ui bug with not properly handling key / props...
                    //https://stackoverflow.com/questions/75818761/material-ui-autocomplete-warning-a-props-object-containing-a-key-prop-is-be
                    return tagValue.map((option, index) => !option ? null : (
                        <Chip
                            {...getTagProps({ index })}
                            key={val(option)}
                            label={label(option)}
                        />
                    ));
                }}
                renderInput={(params) => {
                    const { InputProps } = params;
                    return (
                        <TextField
                            {...params}
                            slotProps={{
                                input: {
                                    ...InputProps,
                                    onPaste: handlePaste,
                                }
                            }}
                            // label="Checkboxes"
                            // placeholder={props.label}
                        />
                    )
                }}
            />
        </>
    );
});
// removed unused DropdownComponent...

const CheckboxComponent = ({ props }: { props: GuiSpec<"check"> }) => (
    <>
        <MLabel props={props} />
        <Checkbox
            checked={props.current_value || false}
            onChange={action((e) => {
                props.current_value = e.target.checked;
                if (props.func) props.func(e.target.checked);
            })}
            size="small"
            sx={{ padding: 0 }}
        />
    </>
);

const RadioButtonComponent = ({
    props,
}: { props: GuiSpec<"radiobuttons"> }) => {
    const choices = useMemo(
        () => props.choices?.map((v) => ({ v, id: uuid() })) || [],
        [props.choices],
    );
    return (
        <>
            <MLabel props={props} />
            <FormControl margin="dense">
                <RadioGroup
                    row
                    // aria-labelledby=""
                    value={props.current_value}
                    onChange={action((e) => {
                        props.current_value = e.target.value;
                        if (props.func) props.func(e.target.value);
                    })}
                >
                    {choices.map(({v, id}) => {
                        const cid = `${id}-${v[0]}`;
                        return (
                            <FormControlLabel key={cid}
                            control={<Radio size="small"/>} label={v[0]} value={v[1]}
                            style={{minWidth: "unset"}}
                            />
                        )
                    })}
                </RadioGroup>
            </FormControl>
            <div className="ciview-radio-group" style={{display: "none"}}>
                {choices.map(({ v, id }) => (
                    <span key={id}>
                        <span className="m-1">{v[0]}</span>
                        <input
                            type="radio"
                            value={v[1]}
                            // biome-ignore lint/suspicious/noDoubleEquals: number == string is ok here
                            checked={v[1] == props.current_value}
                            onChange={action((e) => {
                                props.current_value = e.currentTarget.value;
                                if (props.func)
                                    props.func(e.currentTarget.value);
                            })}
                        />
                    </span>
                ))}
            </div>
        </>
    );
};

const DoubleSliderComponent = ({
    props,
}: { props: GuiSpec<"doubleslider"> }) => (
    <>
        <MLabel props={props} />
        <Slider
            value={props.current_value}
            min={props.min || 0}
            max={props.max || 1}
            onChange={action((_, value) => {
                if (!Array.isArray(value)) {
                    throw "doubleslider callback should have multiple values";
                }
                props.current_value = value as [number, number];
                if (props.func) props.func(value as [number, number]);
            })}
            size="small"
            // sx={{ padding: 0, paddingRight: 10 }}
        />
    </>
);

const ButtonComponent = ({ props }: { props: GuiSpec<"button"> }) => (
    <>
        <MLabel props={props} />
        <Button
            variant="contained"
            onClick={() => {
                //@ts-expect-error button `func(v: never)` - passing `undefined as never` works, is there a nicer way to write this / define GuiValueTypes?
                if (props.func) props.func();
            }}
        >
            {props.label}
        </Button>
    </>
);

const FolderComponent = ({ props }: { props: GuiSpec<"folder"> }) => {
    // add uuid to each setting to avoid key collisions
    const settings = useMemo(
        () => props.current_value.map((setting) => ({ setting, id: uuid() })),
        [props.current_value],
    );
    if (settings.length === 0) return null;
    return (
        <Accordion
            type="single"
            collapsible
            className="w-full col-span-2"
            //uncomment to expand by default
            // defaultValue={props.label}
        >
            <AccordionItem value={props.label}>
                <AccordionTrigger>{props.label}</AccordionTrigger>
                <AccordionContent>
                    {settings.map(({ setting, id }) => (
                        <AbstractComponent key={id} props={setting} />
                    ))}
                </AccordionContent>
            </AccordionItem>
        </Accordion>
    );
};

/**
 * Recursively collects all disposers from a GuiSpec tree.
 * Handles nested folders and collects disposers from all specs in the tree.
 */
export function collectDisposers(specs: AnyGuiSpec[]): Disposer[] {
    const disposers: Disposer[] = [];
    
    function traverse(spec: AnyGuiSpec) {
        // Collect disposers from this spec if it has any
        if (spec._disposers && Array.isArray(spec._disposers)) {
            disposers.push(...spec._disposers);
        }
        
        // If this is a folder, recursively traverse its contents
        if (spec.type === "folder" && Array.isArray(spec.current_value)) {
            spec.current_value.forEach(traverse);
        }
    }
    
    specs.forEach(traverse);
    return disposers;
}

const Components: {
    [K in GuiSpecType]: React.FC<{ props: GuiSpec<K> }>;
} = {
    text: observer(TextComponent),
    textbox: observer(TextBoxComponent),
    slider: observer(SliderComponent),
    spinner: observer(SpinnerComponent),
    dropdown: DropdownAutocompleteComponent,
    // consider having component specifically for column/category selection <<<<
    // the column selection can make use of column groups
    // category selection can have some logic for multitext / tags
    // both can have the ability to reactively update their options based on the current data
    // 'multidropdown': observer(DropdownComponent),
    multidropdown: DropdownAutocompleteComponent,
    check: observer(CheckboxComponent),
    radiobuttons: observer(RadioButtonComponent),
    doubleslider: observer(DoubleSliderComponent),
    button: observer(ButtonComponent),
    folder: observer(FolderComponent),
    column: ColumnSelectionSettingGui,
    multicolumn: ColumnSelectionSettingGui,
} as const;

// how close is this to something we could use from AddChartDialog?
/**
 * A component for rendering a GUI setting.
 */
export const AbstractComponent = observer(
    ({ props }: { props: AnyGuiSpec | GuiSpec<GuiSpecType> }) => {
        // would like to lose this `as` cast - maybe a newer/future typescript might manage it better?
        const Component = Components[props.type] as React.FC<{
            props: typeof props;
        }>;
        if (!(props.type in Components)) {
            //todo use zod to validate the props object
            console.error(`Unknown component type: '${props.type}'`);
            const errorObj = { props, error: `Unknown component type: ${props.type}` };
            return (
                <div className="col-span-3 w-full"> 
                    <ErrorComponentReactWrapper title={errorObj.error} error={{message: errorObj.error}} extraMetaData={errorObj.props} />
                </div>
            );
        }
        return (
            <div className="grid p-2 pr-4 justify-items-start" style={{ gridTemplateColumns: "1fr 2fr" }}>
                <ErrorBoundary
                    FallbackComponent={({error}) =>
                        (
                        <div className="col-span-3 w-full"> 
                            <ErrorComponentReactWrapper 
                                error={{message: error.message, stack: error.stack}} 
                                extraMetaData={props} 
                                title={`Error displaying '${props.label}'. Click to view details`}
                            />
                        </div>
                        )
                    }
                >
                    <Component props={props} />
                </ErrorBoundary>
            </div>
        );
    },
);

export default observer(<T extends BaseConfig,>({ chart }: { chart: BaseChart<T> }) => {
    // Get the raw settings first so we can collect disposers from the same objects
    const rawSettings = useMemo(() => {
        return chart.getSettings();
    }, [chart]);
    
    const settings = useMemo(() => {
        // is the id just for a key in this component, or should the type passed to the component recognise it?
        // for now, I don't think there's a benefit to including it in the type.
        // FolderComponent also makes keys in a similar way that is again only relevant locally I think.
        const settings = rawSettings.map((setting) => ({ setting, id: uuid() }));
        const wrap = { settings };
        makeAutoObservable(wrap);
        return wrap.settings;
    }, [rawSettings]);
    
    // Collect and dispose all disposers when the component unmounts
    // Use the same rawSettings that are used for rendering
    useEffect(() => {
        const disposers = collectDisposers(rawSettings);
        
        // Cleanup function: dispose all collected disposers when dialog closes
        return () => {
            disposers.forEach((disposer) => {
                try {
                    disposer();
                } catch (e) {
                    console.error("Error disposing reaction:", e);
                }
            });
        };
    }, [rawSettings]);
    
    return (
        <ChartProvider chart={chart}>
            <div className="w-full max-h-[80vh]">
                {settings.map(({ setting, id }) => (
                    <AbstractComponent key={id} props={setting} />
                ))}
            </div>
        </ChartProvider>
    );
});
