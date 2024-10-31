import { observer } from "mobx-react-lite";
import { useState, useMemo, useId, useCallback, useEffect } from "react";
import type { Chart, GuiSpec, GuiSpecType } from "../../charts/charts";
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
import type { ColumnSelectionProps } from "./ColumnSelectionComponent";
import ColumnSelectionComponent from "./ColumnSelectionComponent";

export const MLabel = observer(({ props, htmlFor }: { props: GuiSpec<GuiSpecType>, htmlFor?: string }) => (
    <Typography fontSize="small" sx={{alignSelf: "center", justifySelf: "end", paddingRight: 2}}>{props.label}</Typography>
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
            max={props.max || null}
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
/**
 * Wrap the ColumnSelectionComponent in a setting GUI component.
 * 
 * nb, for some weird reason if this is defined in ColumnSelectionComponent.tsx HMR doesn't work...
 */
export const ColumnSelectionSettingGui = observer(({ props }: { props: GuiSpec<"column"> }) => {
    // proably want to change the type of ColumnSelectionProps anyway...
    // perhaps we should be looking at other places where it's used & make them use this,
    // with a different evolution of the API.
    // currently this is not showing the current_value, among other missing features...
    const setSelectedColumn = useCallback(action((v: string) => {
        props.current_value = v;
        props.func?.(v);
    }), []); //as of this writing, biome is right that props is not a dependency

    const props2: ColumnSelectionProps = useMemo(() => ({
        setSelectedColumn,
        // type: props.type,
        multiple: false,
    }), [setSelectedColumn]);
    return (
        <>
            <MLabel props={props} />
            <ColumnSelectionComponent {...props2} />
        </>
    );
});

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

export const DropdownAutocompleteComponent = observer(({
    props,
}: { props: GuiSpec<"dropdown" | "multidropdown"> }) => {
    // todo review 'virtualization' for large lists
    const id = useId();
    const multiple = props.type === "multidropdown";
    // the props.values may be a tuple of [valueObjectArray, textKey, valueKey], or an array of length 1 - [string[]]
    if (!props.values) throw "DropdownAutocompleteComponent requires props.values";
    const useObjectKeys = props.values.length === 3;
    const [_valueObjectArray, labelKey, valueKey] = useObjectKeys ? props.values : ["X", "X", "X"];
    if (!labelKey || !valueKey) throw "DropdownAutocompleteComponent requires labelKey and valueKey for useObjectKeys";
    //todo handle multitext / tags properly.

    // todo think about how this relates to different type logic with {label, value, original}
    // const validVals = useObjectKeys ? valueObjectArray.map(item => item[valueKey]) : valueObjectArray;

    const toOption = useCallback((original) => {
        const label: string = useObjectKeys ? original[labelKey] : original;
        // is value really always a string?
        const value: string = useObjectKeys ? original[valueKey] : original;
        // const id: string = uuid();
        return { label, value, original };
    }, [useObjectKeys, labelKey, valueKey]);
    type OptionType = ReturnType<typeof toOption>;
    const options = useMemo(() => props.values[0].map(toOption), [props.values, toOption]);
    // bit of a faff with sometimes getting a one-item array, sometimes a single item...
    const getSingleOption = useCallback(
        (option: OptionType | OptionType[]) => {
            const a = Array.isArray(option);
            if (a && option.length > 1) {
                console.warn("ideally we shouldn't have to deal with arrays at all here, but we only expect one value when we do");
            }
            return (Array.isArray(option) ? option[0] : option)
        },
        []);
    const single = getSingleOption;
    const label = useCallback((option: OptionType | OptionType[]) => single(option).label, [single]);
    const val = useCallback((option: OptionType | OptionType[]) => single(option).value, [single]);

    // ------
    // deal with cases where the options have changed (e.g. values from a different column)
    // and may be incompatible with props.current_value
    // ------
    // if 'multiple' is true, make sure we get an array even if one / zero values selected.
    // second half of ternary will either be the existing array if 'multiple', or the single value otherwise.
    const v = useMemo(() => (
        multiple && !Array.isArray(props.current_value)
            ? [props.current_value]
            : props.current_value
    ), [props.current_value, multiple]);
    // check that and maybe provide a bit of a type guard
    if (multiple !== Array.isArray(v))
        throw "logical inconsistency - 'multidropdown' value should be coerced to array by now";
    const validVal = useCallback(
        (v: string) => options.some((item) => item.value === v),
        [options],
    );
    const isArray = Array.isArray(v);
    const allValid = isArray ? v.every(validVal) : validVal(v);
    const okValue = allValid ? v : isArray ? v.filter(validVal) : null;
    //map from 'value' string to option object
    const okOption = (Array.isArray(okValue)
        ? okValue.map((v) => options.find((o) => o.value === v))
        : [options.find((o) => o.value === v)])
    // : options.find((o) => o.value === v); //not-multiple...

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
                value={okOption}
                onChange={action((_, value: OptionType) => {
                    //added type annotation above because mobx seems to fluff the inference to `never`
                    if (value === null) return;
                    if (Array.isArray(value)) {
                        // && valueA.length > 1) {
                        const selected = Array.from(value).map(
                            (a: OptionType) => a.value,
                        );
                        props.current_value = selected;
                        if (props.func) props.func(selected);
                        return;
                    }
                    props.current_value = value.value;
                    if (props.func) props.func(value.value);
                })}
                isOptionEqualToValue={(option, value) => single(option).original === single(value).original}
                renderOption={(props, option, { selected }) => {
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
                    return tagValue.map((option, index) => (
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

const DropdownComponent = ({
    props,
}: { props: GuiSpec<"dropdown" | "multidropdown"> }) => {
    const id = useId();
    const [filter, setFilter] = useState("");
    const filterArray = filter.toLowerCase().split(" ");
    const multiple = props.type === "multidropdown";
    const v =
        multiple && !Array.isArray(props.current_value)
            ? [props.current_value]
            : props.current_value;
    // the props.values may be a tuple of [valueObjectArray, textKey, valueKey], or an array of length 1 - [string[]]
    const useObjectKeys = props.values.length === 3;
    const [valueObjectArray, textKey, valueKey] = props.values;
    // for some reason I can't get this to work with useMemo, but it's not particularly heavy - we also don't memoize the children of the dropdown...
    // so if we do find this is expensive, we can definitely optimize better.
    // we just want a string array to filter on to avoid throwing error e.g. if we have a current_value that's not in the dropdown because category changed.
    //todo handle multitext / tags properly.
    const validVals = useObjectKeys
        ? valueObjectArray.map((item) => item[valueKey])
        : valueObjectArray;

    const validVal = useCallback(
        (v: string) => validVals.some((item) => item === v),
        [validVals],
    );
    const isArray = Array.isArray(v);
    const allValid = isArray ? v.every(validVal) : validVal(v);
    const okValue = allValid ? v : isArray ? v.filter(validVal) : null; //not ok after changing category?

    // type E = SelectChangeEvent<string | string[]>;
    type E = { target: { value: string | string[] } }; // material-ui vs native types are different, but compatible enough to use this here
    const handleChange = action((e: E) => {
        const {
            target: { value },
        } = e;
        if (multiple && Array.isArray(value) && value.length > 1) {
            const selected = Array.from(value); // .map(o => o.value);
            props.current_value = selected;
            if (props.func) props.func(selected);
            return;
        }
        props.current_value = value;
        if (props.func) props.func(value);
    });

    return (
        <>
            <MLabel htmlFor={id} props={props} />
            <Select
                size="small"
                id={id}
                multiple={multiple}
                value={okValue}
                className="w-full"
                onChange={handleChange}
            >
                {props.values[0].map((item) => {
                    const text = useObjectKeys ? item[textKey] : item;
                    const value = useObjectKeys ? item[valueKey] : item;
                    const id = uuid();
                    if (
                        filterArray.some((f) => !text.toLowerCase().includes(f))
                    )
                        return null;
                    return (
                        <MenuItem key={id} value={value}>
                            {text}
                        </MenuItem>
                    );
                })}
            </Select>
            <div />
            <input
                type="text"
                value={filter}
                placeholder="Filter options..."
                onChange={(e) => setFilter(e.target.value)}
                className="m-1 pl-1 justify-self-center"
            />
        </>
    );
};

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
        () => props.choices.map((v) => ({ v, id: uuid() })),
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
                if (props.func) props.func(undefined);
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
    // multicolumn: ColumnSelectionSettingGui,
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
    ({ props }: { props: GuiSpec<GuiSpecType> }) => {
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

export default observer(({ chart }: { chart: Chart }) => {
    const settings = useMemo(() => {
        // is the id just for a key in this component, or should the type passed to the component recognise it?
        // for now, I don't think there's a benefit to including it in the type.
        // FolderComponent also makes keys in a similar way that is again only relevant locally I think.
        const settings = chart
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
