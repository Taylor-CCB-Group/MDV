import { observer } from "mobx-react-lite";
import { useState, useMemo, useId, useCallback } from "react";
import type { Chart, GuiSpec, GuiSpecType } from "../../charts/charts";
import { action, makeAutoObservable } from "mobx";
import { ErrorBoundary } from "react-error-boundary";
import {
    Accordion,
    AccordionContent,
    AccordionItem,
    AccordionTrigger,
} from "@/components/ui/accordion"
import { v4 as uuid } from 'uuid';
import { MenuItem, Select } from "@mui/material";

const TextComponent = ({ props }: { props: GuiSpec<'text'> }) => (
    <>
        <label>{props.label}</label>
        <input type="text" value={props.current_value} onChange={action(e => {
            props.current_value = e.target.value;
            if (props.func) props.func(e.target.value);
        })}
            className="w-full"
        />
    </>
);

const TextBoxComponent = ({ props }: { props: GuiSpec<'text'> }) => (
    <>
        <label>{props.label}</label>
        <div />
        <textarea value={props.current_value} onChange={action(e => {
            props.current_value = e.target.value;
            if (props.func) props.func(e.target.value);
        })}
            className="w-full col-span-2"
        />
    </>
);

const SliderComponent = ({ props }: { props: GuiSpec<'slider'> }) => (
    <>
        <label>{props.label}</label>
        <input type="range" value={props.current_value}
            min={props.min || 0}
            max={props.max || 1}
            step={props.step || 0.01}
            onChange={action(e => {
                const value = Number.parseFloat(e.target.value);
                props.current_value = value;
                if (props.func) props.func(value);
            })} />
    </>
);

const SpinnerComponent = ({ props }: { props: GuiSpec<'spinner'> }) => (
    <>
        <label>{props.label}</label>
        <input type="number"
            value={props.current_value}
            min={props.min || 0}
            max={props.max || null}
            step={props.step || 1}
            onChange={action(e => {
                const value = props.current_value = Number.parseInt(e.target.value);
                if (props.func) props.func(value);
            })} />
    </>
)

const DropdownComponent = ({ props }: { props: GuiSpec<'dropdown' | 'multidropdown'> }) => {
    const id = useId();
    const [filter, setFilter] = useState('');
    const filterArray = filter.toLowerCase().split(' ');
    const multiple = props.type === 'multidropdown';
    // const textFilter = (item: string | DropdownMappedValue<string, string>) => {
    //     if (typeof item === 'string') return filter.some(f => !item.toLowerCase().includes(f));
    //     const exclude = filter.some(f => !item.toString().toLowerCase().includes(f));
    //     return true;
    // }
    // if (!v) return <>DropdownComponent: no values</>; // I guess we'll get an ErrorComponent if we don't handle here... shouldn't happen though
    const v = multiple && !Array.isArray(props.current_value) ? [props.current_value] : props.current_value;
    // the props.values may be a tuple of [valueObjectArray, textKey, valueKey], or an array of length 1 - [string[]]
    const useObjectKeys = props.values.length === 3;
    const [valueObjectArray, textKey, valueKey] = props.values;
    // for some reason I can't get this to work with useMemo, but it's not particularly heavy - we also don't memoize the children of the dropdown...
    // so if we do find this is expensive, we can definitely optimize better.
    // we just want a string array to filter on to avoid throwing error e.g. if we have a current_value that's not in the dropdown because category changed.
    //todo handle multitext / tags properly.
    const validVals = useObjectKeys ? valueObjectArray.map(item => item[valueKey]) : valueObjectArray;

    const validVal = useCallback((v: string) => validVals.some(item => item === v), [validVals]);
    const isArray = Array.isArray(v);
    const allValid = isArray ? v.every(validVal) : validVal(v);
    const okValue = allValid ? v : (isArray ? v.filter(validVal) : null); //not ok after changing category?

    // type E = SelectChangeEvent<string | string[]>;
    type E = { target: { value: string | string[] } }; // material-ui vs native types are different, but compatible enough to use this here
    const handleChange = action((e: E) => {
        const {
            target: { value },
        } = e;
        if (multiple && Array.isArray(value) && value.length > 1) {
            const selected = Array.from(value);// .map(o => o.value);
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
                // https://github.com/snakesilk/react-fullscreen/issues/44
                // MenuProps={{container: () => document.getElementById('fullscreen-node')}} //todo maybe a hook to get the fullscreen node or document.body
                size="small"
                id={id}
                multiple={multiple}
                value={okValue}
                className="w-full"
                onChange={handleChange}>
                {props.values[0].map((item) => {
                    const text = useObjectKeys ? item[textKey] : item;
                    const value = useObjectKeys ? item[valueKey] : item;
                    const id = uuid();
                    if (filterArray.some(f => !text.toLowerCase().includes(f))) return null;
                    return <MenuItem key={id} value={value}>{text}</MenuItem>
                })}
            </Select>
            <div />
            <input type="text" value={filter} placeholder="Filter options..."
                onChange={(e => setFilter(e.target.value))}
                className="m-1 pl-1 justify-self-center"
            />
        </>
    )
};

const CheckboxComponent = ({ props }: { props: GuiSpec<'check'> }) => (
    <>
        <label>{props.label}</label>
        <input type="checkbox" checked={props.current_value} onChange={action(e => {
            props.current_value = e.target.checked;
            if (props.func) props.func(e.target.checked);
        })} />
    </>
);

const RadioButtonComponent = ({ props }: { props: GuiSpec<'radiobuttons'> }) => {
    const choices = useMemo(() => (
        props.choices.map(v => ({ v, id: uuid() }))
    ), [props.choices]);
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
                        <span
                            className="m-1"
                        >
                            {v[0]}
                        </span>
                        <input
                            // biome-ignore lint/suspicious/noDoubleEquals: number == string is ok here
                            type="radio" value={v[1]} checked={v[1] == props.current_value} onChange={action(e => {
                                props.current_value = e.currentTarget.value;
                                if (props.func) props.func(e.currentTarget.value);
                            })} />
                    </span>
                ))}
            </div>
        </>
    )
};

const DoubleSliderComponent = ({ props }: { props: GuiSpec<'doubleslider'> }) => (
    <>
        <label>{props.label}</label>
        <input type="range" value={props.current_value[0]} onChange={action(e => {
            const v = props.current_value[0] = Number.parseFloat(e.target.value);
            if (props.func) props.func([v, props.current_value[1]]);
        })} />
        <input type="range" value={props.current_value[1]} onChange={action(e => {
            const v = props.current_value[1] = Number.parseFloat(e.target.value);
            if (props.func) props.func([props.current_value[0], v]);
        })} />
    </>
);

const ButtonComponent = ({ props }: { props: GuiSpec<'button'> }) => (
    <>
        <button type="button" onClick={action(e => {
            if (props.func) props.func(undefined);
        })}>{props.label}</button>
    </>
);

const FolderComponent = ({ props }: { props: GuiSpec<'folder'> }) => {
    // add uuid to each setting to avoid key collisions
    const settings = useMemo(() => (
        props.current_value.map(setting => ({ setting, id: uuid() }))
    ), [props.current_value]);
    if (settings.length === 0) return null;
    return (
        <Accordion type='single' collapsible className="w-full col-span-2">
            <AccordionItem value={props.label}>
                <AccordionTrigger>{props.label}</AccordionTrigger>
                <AccordionContent>
                    {settings.map(({ setting, id }) => <AbstractComponent key={id} props={setting} />)}
                </AccordionContent>
            </AccordionItem>
        </Accordion>
    )
}

const Components: Record<GuiSpecType, React.FC<{ props: GuiSpec<GuiSpecType> }>> = {
    'text': observer(TextComponent),
    'textbox': observer(TextBoxComponent),
    'slider': observer(SliderComponent),
    'spinner': observer(SpinnerComponent),
    'dropdown': observer(DropdownComponent),
    'multidropdown': observer(DropdownComponent),
    'check': observer(CheckboxComponent),
    'radiobuttons': observer(RadioButtonComponent),
    'doubleslider': observer(DoubleSliderComponent),
    'button': observer(ButtonComponent),
    'folder': observer(FolderComponent),
} as const;

const ErrorComponent = ({ props }: { props: GuiSpec<GuiSpecType> }) => {
    const [expanded, setExpanded] = useState(false);
    return (
        <div className="border-red-500 border-2 border-solid" onClick={() => setExpanded(!expanded)}>
            <h2>Error...</h2>
            <pre>{JSON.stringify(props, null, 2).substring(0, expanded ? undefined : 100)}</pre>
        </div>
    )
}

const AbstractComponent = observer(({ props }: { props: GuiSpec<GuiSpecType> }) => {
    const Component = Components[props.type];
    return (
        <div className="grid grid-cols-2 p-1 justify-items-start">
            <ErrorBoundary fallback={<ErrorComponent props={props} />}>
                <Component props={props} />
            </ErrorBoundary>
        </div>
    )
})

export default observer(({ chart }: { chart: Chart }) => {
    const settings = useMemo(() => {
        const settings = chart.getSettings().map(setting => ({ setting, id: uuid() }));
        const wrap = { settings };
        makeAutoObservable(wrap);
        return wrap.settings;
    }, [chart]);
    return (
        <div className="w-full">
            {settings.map(({ setting, id }) => <AbstractComponent key={id} props={setting} />)}
        </div>
    )
})