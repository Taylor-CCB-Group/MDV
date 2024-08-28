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
import { Button, Chip, MenuItem, Select, Slider } from "@mui/material";
import Checkbox from "@mui/material/Checkbox";
import TextField from "@mui/material/TextField";
import Autocomplete from "@mui/material/Autocomplete";
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";

const TextComponent = ({ props }: { props: GuiSpec<"text"> }) => (
    <>
        <label>{props.label}</label>
        <input
            type="text"
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
        <label>{props.label}</label>
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
        <label>{props.label}</label>
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
        />
    </>
);

const SpinnerComponent = ({ props }: { props: GuiSpec<"spinner"> }) => (
    <>
        <label>{props.label}</label>
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

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

function DropdownAutocompleteComponent({
    props,
}: { props: GuiSpec<"dropdown" | "multidropdown"> }) {
    // todo review 'virtualization' for large lists
    const id = useId();
    const multiple = props.type === "multidropdown";
    if (!multiple)
        console.warn(
            "DropdownAutocompleteComoponet with non-multi dropdown is WIP",
        );

    // the props.values may be a tuple of [valueObjectArray, textKey, valueKey], or an array of length 1 - [string[]]
    const useObjectKeys = props.values.length === 3;
    const [valueObjectArray, textKey, valueKey] = props.values;

    //todo handle multitext / tags properly.

    // todo think about how this relates to different type logic with {text, value, original}
    // const validVals = useObjectKeys ? valueObjectArray.map(item => item[valueKey]) : valueObjectArray;

    const toOption = (original) => {
        const text: string = useObjectKeys ? original[textKey] : original;
        // is value really always a string?
        const value: string = useObjectKeys ? original[valueKey] : original;
        // const id: string = uuid();
        return { text, value, original };
    };
    type OptionType = ReturnType<typeof toOption>;
    const options = props.values[0].map(toOption);

    // ------
    // deal with cases where the options have changed (e.g. values from a different column)
    // and may be incompatible with props.current_value
    // ------
    // if 'multiple' is true, make sure we get an array even if one / zero values selected.
    // second half of ternary will either be the existing array if 'multiple', or the single value otherwise.
    const v =
        multiple && !Array.isArray(props.current_value)
            ? [props.current_value]
            : props.current_value;
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
    const okOption = Array.isArray(okValue)
        ? okValue.map((v) => options.find((o) => o.value === v))
        : [options.find((o) => o.value === v)];

    return (
        <>
            <label className="align-middle" htmlFor={id}>
                {props.label}
            </label>
            <Autocomplete
                className="w-full"
                multiple={multiple}
                size="small"
                id={id}
                options={options}
                disableCloseOnSelect={multiple}
                getOptionLabel={(option) => option.text}
                value={okOption}
                onChange={action((_, value: OptionType) => {
                    //added type annotation because mobx seems to fluff the inference
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
                // isOptionEqualToValue={(option, value) => option.original === value.original}
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
                                {option.text}
                            </li>
                        );
                    return (
                        <li key={key} {...optionProps}>
                            {option.text}
                        </li>
                    );
                }}
                renderTags={(tagValue, getTagProps) => {
                    //seems to be a material-ui bug with not properly handling key / props...
                    //https://stackoverflow.com/questions/75818761/material-ui-autocomplete-warning-a-props-object-containing-a-key-prop-is-be
                    return tagValue.map((option, index) => (
                        <Chip
                            {...getTagProps({ index })}
                            key={option.value}
                            label={option.text}
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
}

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
            <label htmlFor={id}>{props.label}</label>
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
        <label>{props.label}</label>
        <input
            type="checkbox"
            checked={props.current_value}
            onChange={action((e) => {
                props.current_value = e.target.checked;
                if (props.func) props.func(e.target.checked);
            })}
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
            <label>{props.label}</label>
            {/* <RadioGroup << meh
            className="ciview-radio-group"
            defaultValue={props.current_value}
                onValueChange={action(e => {
                    props.current_value = e;
                    if (props.func) props.func(e);
                })}
            >
                {props.choices.map((v, i) => {
                    const cid = `${id}-${v[0]}`;
                    return (
                        <span key={cid}>
                            <RadioGroupItem id={cid} value={v[1]} />
                            <Label htmlFor={cid}>{v[0]}</Label>
                        </span>
                    )
                })}
            </RadioGroup> */}
            <div className="ciview-radio-group">
                {choices.map(({ v, id }) => (
                    <span key={id}>
                        <span className="m-1">{v[0]}</span>
                        <input
                            // biome-ignore lint/suspicious/noDoubleEquals: number == string is ok here
                            type="radio"
                            value={v[1]}
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
        <label>{props.label}</label>
        <Slider
            value={props.current_value}
            min={props.min || 0}
            max={props.max || 1}
            onChange={action((_, value) => {
                if (!Array.isArray(value))
                    throw "doubleslider callback should have multiple values";
                props.current_value = value as [number, number];
                if (props.func) props.func(value as [number, number]);
            })}
        />
    </>
);

const ButtonComponent = ({ props }: { props: GuiSpec<"button"> }) => (
    <>
        <label>{props.label}</label>
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
    dropdown: observer(DropdownComponent), //todo also use Autocomplete for this
    // consider having component specifically for column/category selection
    // the column selection can make use of column groups
    // category selection can have some logic for multitext / tags
    // both can have the ability to reactively update their options based on the current data
    // 'multidropdown': observer(DropdownComponent),
    multidropdown: observer(DropdownAutocompleteComponent),
    check: observer(CheckboxComponent),
    radiobuttons: observer(RadioButtonComponent),
    doubleslider: observer(DoubleSliderComponent),
    button: observer(ButtonComponent),
    folder: observer(FolderComponent),
} as const;

const ErrorComponent = ({ props }: { props: GuiSpec<GuiSpecType> }) => {
    const [expanded, setExpanded] = useState(false);
    return (
        <div
            className="border-red-500 border-2 border-solid"
            onClick={() => setExpanded(!expanded)}
        >
            <h2>Error...</h2>
            <pre>
                {JSON.stringify(props, null, 2).substring(
                    0,
                    expanded ? undefined : 100,
                )}
            </pre>
        </div>
    );
};

const AbstractComponent = observer(
    ({ props }: { props: GuiSpec<GuiSpecType> }) => {
        // would like to lose this `as` cast - maybe a newer/future typescript might manage it better?
        const Component = Components[props.type] as React.FC<{
            props: typeof props;
        }>;
        return (
            <div className="grid grid-cols-2 p-1 justify-items-start">
                <ErrorBoundary fallback={<ErrorComponent props={props} />}>
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
        <div className="w-full">
            {settings.map(({ setting, id }) => (
                <AbstractComponent key={id} props={setting} />
            ))}
        </div>
    );
});
