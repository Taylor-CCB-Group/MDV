import { observer } from "mobx-react-lite";
import { useState, useMemo, useId } from "react";
import type { Chart, GuiSpec, GuiSpecType } from "../../charts/charts";
import { action, makeAutoObservable } from "mobx";
import { ErrorBoundary } from "react-error-boundary";
import {
    Accordion,
    AccordionContent,
    AccordionItem,
    AccordionTrigger,
} from "@/components/ui/accordion"

const TextComponent = function({props}: {props: GuiSpec<'text'>}) {
    return (
        <>
            <label>{props.label}</label>
            <input type="text" value={props.current_value} onChange={action(e => {
                props.current_value = e.target.value;
                if (props.func) props.func(e.target.value);
            })} />
        </>
    )
};

const SliderComponent = function({props}: {props: GuiSpec<'slider'>}) {
    return (
        <>
            <label>{props.label}</label>
            <input type="range" value={props.current_value} 
            min={props.min || 0}
            max={props.max || 1}
            step={props.step || 0.01}
            onChange={action(e => {
                const value = parseFloat(e.target.value);
                props.current_value = value;
                if (props.func) props.func(value);
            })} />
        </>
    )
};

const SpinnerComponent = function({props}: {props: GuiSpec<'spinner'>}) {
    return (
        <>
            <label>{props.label}</label>
            <input type="number"
            value={props.current_value} 
            min={props.min || 0}
            max={props.max || null}
            step={props.step || 1}
            onChange={action(e => {
                const value = props.current_value = parseInt(e.target.value);
                if (props.func) props.func(value);
            })} />
        </>
    )
}

const DropdownComponent = function({props}: {props: GuiSpec<'dropdown' | 'multidropdown'>}) {
    const id = useId();
    const [filter, setFilter] = useState('');
    const filterArray = filter.toLowerCase().split(' ');
    // const textFilter = (item: string | DropdownMappedValue<string, string>) => {
    //     if (typeof item === 'string') return filter.some(f => !item.toLowerCase().includes(f));
    //     const exclude = filter.some(f => !item.toString().toLowerCase().includes(f));
    //     return true;
    // }
    if (!props.values) return <>DropdownComponent: no values</>;
    const v = props.type === 'multidropdown' && !Array.isArray(props.current_value) ? [props.current_value] : props.current_value;
    return (
        <>
            <label htmlFor={id}>{props.label}</label>
            <select 
            id={id}
            multiple={props.type === 'multidropdown'}
            value={v}
            className="w-full"
            onChange={action(e => {
                props.current_value = e.target.value;
                if (props.func) props.func(e.target.value);
            })}>
                {props.values[0].map((item, i) => {
                    const s = props;
                    const text = s.values.length > 1 ? item[s.values[1]] : item;
                    const value = s.values.length > 1 ? item[s.values[2]] : item;
                    if (filterArray.some(f => !text.toLowerCase().includes(f))) return null;
                    return <option key={i} value={value}>{text}</option>
                })}
            </select>
            <div></div>
            <input type="text" value={filter} placeholder="Filter options..." 
            onChange={(e => setFilter(e.target.value))} 
            className="m-1 pl-1 justify-self-center"
            />
        </>
    )
};

const CheckboxComponent = function({props}: {props: GuiSpec<'check'>}) {
    return (
        <>
            <label>{props.label}</label>
            <input type="checkbox" checked={props.current_value} onChange={action(e => {
                props.current_value = e.target.checked;
                if (props.func) props.func(e.target.checked);
            })}/>
        </>
    )
};

const RadioButtonComponent = function({props}: {props: GuiSpec<'radiobuttons'>}) {
    const id = useId();
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
                {props.choices.map((v, i) => (
                    <>
                    <span 
                    className="m-1"
                    key={i + id}>
                        {v[0]}
                    </span>
                        <input key={i + id + 'input'}
                        type="radio" value={v[1]} checked={v[1] === props.current_value} onChange={action(e => {
                            props.current_value = e.currentTarget.value;
                            if (props.func) props.func(e.currentTarget.value);
                        })}/>
                    </>
                ))}
            </div>
        </>
    )
};

const DoubleSliderComponent = function({props}: {props: GuiSpec<'doubleslider'>}) {
    return (
        <>
            <label>{props.label}</label>
            <input type="range" value={props.current_value[0]} onChange={action(e => {
                const v = props.current_value[0] = parseFloat(e.target.value);
                if (props.func) props.func([v, props.current_value[1]]);
            })} />
            <input type="range" value={props.current_value[1]} onChange={action(e => {
                const v = props.current_value[1] = parseFloat(e.target.value);
                if (props.func) props.func([props.current_value[0], v]);
            })} />
        </>
    )
};

const ButtonComponent = function({props}: {props: GuiSpec<'button'>}) {
    return (
        <>
            <button onClick={action(e => {
                if (props.func) props.func(undefined);
            })}>{props.label}</button>
        </>
    )
};

const FolderComponent = function({props}: {props: GuiSpec<'folder'>}) {
    return (
        <Accordion type='single' collapsible className="w-full col-span-2">
            <AccordionItem value={props.label}>
                <AccordionTrigger>{props.label}</AccordionTrigger>
                <AccordionContent>
                    {props.current_value.map((setting, i) => <AbstractComponent key={i} props={setting} />)}
                </AccordionContent>
            </AccordionItem>
        </Accordion>
    )
}

const Components: Record<GuiSpecType, React.FC<{props: GuiSpec<GuiSpecType>}>> = {
    'text': observer(TextComponent),
    'textbox': observer(TextComponent),
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

const ErrorComponent = function({props}: {props: GuiSpec<GuiSpecType>}) {
    const [expanded, setExpanded] = useState(false);
    return (
        <div className="border-red-500 border-2 border-solid" onClick={e => setExpanded(!expanded)}>
            <h2>Error...</h2>
            <pre>{JSON.stringify(props, null, 2).substring(0, expanded ? undefined : 100)}</pre>
        </div>
    )
}

const AbstractComponent = observer(function({props}: {props: GuiSpec<GuiSpecType>}) {
    const Component = Components[props.type];
    return (
        <div className="grid grid-cols-2 p-1 justify-items-start">
            <ErrorBoundary fallback={<ErrorComponent props={props} />}>
                <Component props={props} />
            </ErrorBoundary>
        </div>
    )
})

export default observer(function({chart}: {chart: Chart}) {
    const settings = useMemo(() => {
        const settings = chart.getSettings();
        const wrap = {settings};
        makeAutoObservable(wrap);
        return wrap.settings;
    }, [chart]);
    return (
        <div className="w-full">
            {settings.map((setting, i) => <AbstractComponent key={i} props={setting} />)}
        </div>
    )
})