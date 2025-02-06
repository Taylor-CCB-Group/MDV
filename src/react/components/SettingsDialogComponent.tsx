import { observer } from "mobx-react-lite";
import { useState, useMemo, useId, useCallback, useEffect } from "react";
import type { AnyGuiSpec, DropDownValues, GuiSpec, GuiSpecType } from "../../charts/charts";
import { action, makeAutoObservable } from "mobx";
import { ErrorBoundary } from "react-error-boundary";
import {
    Accordion,
    AccordionContent,
    AccordionItem,
    AccordionTrigger,
} from "@/components/ui/accordion";
import { v4 as uuid } from "uuid";
import { Button, Chip, Dialog, FormControl, FormControlLabel, MenuItem, Radio, RadioGroup, Select, Slider, Typography } from "@mui/material";
import Checkbox from "@mui/material/Checkbox";
import TextField from "@mui/material/TextField";
import Autocomplete from "@mui/material/Autocomplete";
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import JsonView from "react18-json-view";
import { ChartProvider } from "../context";
import ColumnSelectionComponent from "./ColumnSelectionComponent";
import { inferGenericColumnSelectionProps } from "@/lib/columnTypeHelpers";
import { g, isArray, notEmpty } from "@/lib/utils";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";

export const MLabel = observer(({ props, htmlFor }: { props: AnyGuiSpec, htmlFor?: string }) => (
    <Typography fontSize="small" sx={{alignSelf: "center", justifySelf: "end", paddingRight: 2}}>
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
                    throw "slider callback should have multiple values";
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
 * Wrap the ColumnSelectionComponent in a setting GUI component.
 *
 * nb, for some weird reason if this is defined in ColumnSelectionComponent.tsx HMR doesn't work...
 */
export const ColumnSelectionSettingGui = observer(({ props }: { props: ColumnSelectionSpec }) => {
    // multiple...
    // proably want to change the type of ColumnSelectionProps anyway...
    // perhaps we should be looking at other places where it's used & make them use this,
    // with a different evolution of the API.
    // currently this is not showing the current_value, among other missing features...
    /** this needs fixing... and we should be able to observe mobx state for column queries 
     * which should persist when the dialog (entire react root) is closed.
    */

    const setSelectedColumn = useCallback(action((v: string) => {
        props.current_value = v;
        //@ts-expect-error string is not assignable to FieldSpecs, type of v is wrong here
        props.func?.(v);
    }), []); //as of this writing, biome is right that props is not a dependency
    const filter = props.columnSelection?.filter;
    // not only is the filter not working, but we need to decide how to express "multiple"
    const props2 = useMemo(() => inferGenericColumnSelectionProps({
        // fixing this stuff is high priority
        //@ts-expect-error ColumnSelection setSelectedColumn type
        setSelectedColumn,
        //@ts-expect-error ColumnSelection `type` type
        type: filter,
        multiple: props.type === "multicolumn",
        //@ts-expect-error ColumnSelection current_value type
        current_value: props.current_value
        // current_value: props.current_value... maybe want to be more mobx-y about this
    }), [setSelectedColumn, props.type, filter, props.current_value]);
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
    if (!labelKey || !valueKey) throw "DropdownAutocompleteComponent requires labelKey and valueKey for useObjectKeys";
    type Option = { label: string, value: string, original: any };
    // might put some more elaborate types here on original at some point
    return (original: string | any): Option => {
        const isString = typeof original === "string";
        if (isString === useObjectKeys) throw "DropdownAutocompleteComponent requires values to be either all strings or all objects";
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
    if (!props.values) throw "DropdownAutocompleteComponent requires props.values - probable logic error in type definitions";
    const toOption = useCallback(getOptionAsObjectHelper(props.values), []);
    type OptionType = ReturnType<typeof toOption>;
    const NO_OPTIONS = [{ label: "No options (error - sorry...)", value: "", original: "" }];
    const options = useMemo(() => props.values?.[0].map(toOption) || NO_OPTIONS, [props.values, toOption]);
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
    // some type checking / evaluation
    // if (props.type === "multidropdown") {
    //     const t: "multidropdown" = props.type;
    //     const arr: string[] = props.current_value;
    //     const func: GuiSpec<"multidropdown">['func'] = props.func;
    // } else {
    //     const t: "dropdown" = props.type;
    //     //@ts-expect-error
    //     const arr: string = v;
    //     const arr2: string = props.current_value;
    //     if (props.func) {
    //         const func: (v: string) => void = props.func;
    //     }
    // }
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
    if (multi !== isArray(props.current_value)) throw "current_value should be an array if (and only if) type is multidropdown";
    return multi;
}
function inferDropdownType<T extends "dropdown" | "multidropdown">(props: GuiSpec<T>): GuiSpec<T> {
    return props;
}
function getCurrentValue(props: DropdownSpec) {
    const multi = isMultidropdown(props);
    if (multi) return props.current_value; //ts correctly infers this is an array
    // return props.current_value; //.. fails to infer this is a string
    if (isArray(props.current_value)) throw "current_value should be an array if and only if type is multidropdown";
    return props.current_value;
}
export const DropdownAutocompleteComponent = observer(({
    props//: propsAmbiguous,
}: { props: DropdownSpec }) => {
    // todo review 'virtualization' for large lists
    const id = useId();
    // const props = inferDropdownType(propsAmbiguous);
    const multiple = isMultidropdown(props);
    // the props.values may be a tuple of [valueObjectArray, textKey, valueKey], or an array of length 1 - [string[]]
    if (!props.values) throw "DropdownAutocompleteComponent requires props.values";
    //todo handle multitext / tags properly.

    const {options, single, label, val, okOption} = useDropdownOptions(props);

    type Option = typeof options[number]; // hopefully we can improve `{ original: any }`
    // still not entirely sure about the type for onChange...
    type OVal = Option | Option[] | (Option | Option[])[] | null;
    return (
        <>
            <MLabel htmlFor={id} props={props} />
            <Autocomplete
                className="w-full"
                multiple={multiple}
                size="small"
                id={id}
                options={options}
                disableCloseOnSelect={multiple}
                getOptionLabel={label}
                value={okOption.filter((a) => a !== undefined)}
                onChange={action((_, value: OVal) => {
                    //added type annotation above because mobx seems to fluff the inference to `never`
                    if (value === null) return;
                    if (isArray(value)) {
                        if (!multiple) throw "shouldn't get an array here";
                        const vals = value.map(val);
                        props.current_value = vals;
                        props.func?.(vals);
                        return;
                    }
                    if (multiple) throw "shouldn't get a single value here";
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
                renderInput={(params) => (
                    <TextField
                        {...params}
                    // label="Checkboxes"
                    // placeholder={props.label}
                    />
                )}
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
            //expand by default
            defaultValue={props.label}
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

const ErrorComponent = observer(({ props, label }: { props: any, label: string }) => {
    const [expanded, setExpanded] = useState(true);
    return (
        <>
        <Dialog open={expanded}
            className="border-red-500 border-2 border-solid"
            >
            <h2>Settings Dialog Error...</h2>
            <JsonView src={props} collapsed={1} />
            <Button
                onClick={() => setExpanded(!expanded)}
                >ok</Button>
        </Dialog>
        <label className="text-red-500">Error in component:</label>
        <label className="text-red-500">{label}</label>
        </>
    );
});

// how close is this to something we could use from AddChartDialog?
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
            return <ErrorComponent props={errorObj} label={props.label} />;
        }
        return (
            <div className="grid grid-cols-2 p-1 justify-items-start">
                <ErrorBoundary fallback={<ErrorComponent props={props} label={props.label} />}>
                    <Component props={props} />
                </ErrorBoundary>
            </div>
        );
    },
);

export default observer(<T extends BaseConfig,>({ chart }: { chart: BaseChart<T> }) => {
    const settings = useMemo(() => {
        // is the id just for a key in this component, or should the type passed to the component recognise it?
        // for now, I don't think there's a benefit to including it in the type.
        // FolderComponent also makes keys in a similar way that is again only relevant locally I think.
        const settings = chart //todo: fix this typecast, resolve `Chart` vs `BaseChart`...
            .getSettings()
            .map((setting) => ({ setting, id: uuid() }));
        const wrap = { settings };
        makeAutoObservable(wrap);
        return wrap.settings;
    }, [chart]);
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
