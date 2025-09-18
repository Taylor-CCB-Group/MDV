import { useCallback, useEffect, useMemo, useState } from "react";
import type DataStore from "../../datastore/DataStore.js";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import { observer } from "mobx-react-lite";
import TextField from "@mui/material/TextField";
import Autocomplete from "@mui/material/Autocomplete";

import BaseChart from "../BaseChart";
import { DataStoreContext, useDataStore } from "@/react/context.js";
import JsonView from "react18-json-view";
import { Button, Dialog, DialogActions, DialogContent, DialogTitle, Divider, IconButton, Paper, Typography } from "@mui/material";
import Grid from "@mui/material/Grid2";
import ColumnSelectionComponent from "@/react/components/ColumnSelectionComponent.js";
import { action, observable, reaction, toJS } from "mobx";
import z from "zod/v4";
import type { DataColumn, DataType, ExtraControl, GuiSpec, GuiValueTypes } from "../charts.js";
import { AbstractComponent } from "@/react/components/SettingsDialogComponent.js";
import { columnMatchesType, isMultiColumn } from "@/lib/columnTypeHelpers";
import type { Param } from "@/lib/columnTypeHelpers.js";
import { g, isArray } from "@/lib/utils.js";
import { serialiseConfig } from "../chartConfigUtils.js";
import "react18-json-view/src/style.css";
// If dark mode is needed, import `dark.css`.
import "react18-json-view/src/dark.css";
import "../../utilities/css/JsonDialogStyles.css";
import { Close as CloseIcon } from "@mui/icons-material";
import { ErrorBoundary } from "react-error-boundary";
import ErrorComponentReactWrapper from "@/react/components/ErrorComponentReactWrapper.js";

// how do I match this with the type in link_utils?
//! using zod here is an extra layer of complexity and it's not at all clear that this is the place for it
// we're not validating data coming in, we're just using it to define the shape of the data 
// - the whole schema-isation is really a separate concern that I very much want to have - but in the right place/time
// we don't want a schema for the live version of the object, only for the serialised version
// const RowsAsColsQuerySchema = z.object({});

const ParamzSchema = z.optional(z.array(z.union([z.string(), z.array(z.string())])));
type Paramz = z.infer<typeof ParamzSchema>;

const ChartConfigSchema = z.object({
    title: z.string(),
    legend: z.string(),
    param: ParamzSchema,
    type: z.string(),
    // in the original AddChartDialog, extra props are on the root object, not nested
    // so when we pass this to ChartManager, we'll need to flatten it
    extra: z.record(z.string(), z.unknown()),
    _updated: z.date(),
});
type ChartConfig = z.infer<typeof ChartConfigSchema>;


function flattenParams(param: Paramz) {
    if (!param) return []; //would quite like it to remain undefined
    //! todo change the type so that it can also include column-query... probably not defined by zod
    const result = param?.flatMap(p => isArray(p) ? [...p] : p);
    return result;
}

const ChartPreview = observer(({config}: {config: ChartConfig}) => {
    const dataStore = useDataStore();
    const chartType = useMemo(() => {
        return Object.values(BaseChart.types).find(t => t.name === config.type);
    }, [config.type]);
    // this could be an actual chart preview, but for now, the config is useful
    // will probably keep it as an option to show the JSON, but default to a chart preview
    // this should always have the result of any extra controls & init applied - which may be an error
    // todo rearrange so that the same code is used for adding actual chart
    const scratchProps = useMemo(() => {
        const scratchConfig = toJS(config);
        scratchConfig.param = flattenParams(scratchConfig.param);
        const { extra, _updated, ...props } = scratchConfig;
        // flatten the config into props
        chartType?.extra_controls?.(dataStore).forEach(control => {
            if (control.defaultVal && !extra[control.name]) {
                extra[control.name] = control.defaultVal;
            } else {
                // if (!control.defaultVal)
            }
        });
        for (const [k, v] of Object.entries(extra)) {
            //@ts-ignore
            props[k] = v;
        }
        // find the key in BaseChart.types that corresponds to this chart type
        // I don't much like this, but it seems to vaguely work for now
        for (const [k, v] of Object.entries(BaseChart.types)) {
            if (v.name === config.type) {
                props.type = k;
                break;
            }
        }
        if (chartType?.init) {
            try {
                chartType.init(props, dataStore, extra);
            } catch (e) {
                return e;
            }
        }
        return props;
    }, [config, dataStore, chartType]);
    const { _updated, ...configWithoutUpdated } = toJS(config);
    return (
        <div className="grid grid-cols-2 px-3 py-5">
            {/* <h2>Chart Preview</h2> */}
            <JsonView style={{
                    backgroundColor: 'transparent',
                    fontSize: '0.875rem',
                    fontFamily: 'ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace',
                }} src={configWithoutUpdated} collapseStringsAfterLength={10} collapsed={1} />
            <JsonView style={{
                    backgroundColor: 'transparent',
                    fontSize: '0.875rem',
                    fontFamily: 'ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace',
                }} src={scratchProps} collapsed={1} />
        </div>
    );
});

/** given an (observable) ExtraControl, return an (observable) GuiSpec
 * adapted so we can use SettingsDialogComponent widgets
 * When the returned GuiSpec is updated, the property of config corresponding to the control is updated
*/
function controlToGuiSpec<T extends keyof GuiValueTypes>(control: ExtraControl<T>, onChange: (v?: GuiValueTypes[T]) => void): GuiSpec<T> {
    const draftSpec = {
        type: control.type,
        label: control.label,
        current_value: control.defaultVal,
        // values: control.values,
    } satisfies Partial<GuiSpec<T>>;
    if (draftSpec.type === "dropdown" || draftSpec.type === "multidropdown") {
        // we run into this irritating logic with format that dropdown values can come in
        // >the props.values may be a tuple of [valueObjectArray, textKey, valueKey],
        // or an array of length 1 - [string[]]
        //const useObjectKeys = Array.isArray(control.values) && control.values.length === 3;
        if (!control.values || control.values.length === 0) throw new Error(`Dropdown control '${control.label}' has no values`);
        const values = control.values;
        // if (!control.defaultVal) draftSpec.current_value = useObjectKeys ? values[0][0][values[0][2]] : values[0][0];
        // todo - add some type predicate so we know we're assigning the right type
        //@ts-expect-error controlToGuiSpec WIP
        if (!control.defaultVal) draftSpec.current_value = values[0]["name"];
        //(draftSpec as any).values = control.values; //doesn't give runtime error, doesn't work
        (draftSpec as any).values = [values, "name", "value"];
    }
    if (draftSpec.type === "radiobuttons") {}//todo test
    //bad bad bad <<< how did we end up needing "as string", "as any"???
    //if (draftSpec.type === "check" && !control.defaultVal) (draftSpec as any).current_value = false;
    // if ((draftSpec as any).type === "checkbox" && !control.defaultVal) (draftSpec as any).current_value = false;
    if ((draftSpec.type as string) === "checkbox") (draftSpec as any).type = "check";
    const spec = observable(draftSpec);
    // if (!spec.current_value) throw new Error(`current_value should be set - no default in '${JSON.stringify(control)}'`);
    // strict typescript warns us about spec.current_value being potentially undefined, and it's dead right.
    // IMO much of the faffing around I'm doing getting this all to work could be avoided with better types.
    reaction(() => spec.current_value, onChange);
    //@ts-expect-error controlToGuiSpec WIP
    return g(spec);
}
function useGuiSpecControl<T extends keyof GuiValueTypes>(control: ExtraControl<T>, config: ChartConfig): GuiSpec<T> {
    const spec = useMemo(() => {
        const val = controlToGuiSpec(control, (v) => {config.extra[control.name] = v});
        config.extra[control.name] = val.current_value;
        return val;
    }, [control, config]);
    return spec;
}

const ExtraControlComponent = observer(({ control, config }: { control: ExtraControl<keyof GuiValueTypes>, config: ChartConfig }) => {
    const spec = useGuiSpecControl(control, config);
    return (
        <Paper elevation={6}>
            <AbstractComponent props={spec} />
        </Paper>
    );
});


const ConfigureChart = observer(({config, onDone}: {config: ChartConfig, onDone: () => void}) => {
    const dataStore = useDataStore();
    const chartTypes = useMemo(() => {
        return Object.values(BaseChart.types).filter(t => {
            if (t.allow_user_add === false) return false;
            if (t.required) {
                console.log(t.name, "required", t.required);
                if (typeof t.required === "function") {
                    return t.required(dataStore);
                }
                // nb `t.required.every(r => r in dataStore)` seems like it should work,
                // but since we're in glorious JS land, we have things like `dataStore['interactions'] = undefined`
                // which means that `'interactions' in dataStore` is true, but `!!dataStore['interactions']` is false.
                //@ts-ignore DataStore string index on things 'required' for chart types
                return t.required.every(r => !!dataStore[r]);
            }
            return true;
        });
    }, [dataStore]);
    const chartNames = chartTypes.map(t => t.name).sort((a, b) => a.localeCompare(b));
    // const [chartTypeName, setChartTypeName] = useState(chartNames[0]);
    const paramColumns = useMemo(
        () => config.param ? flattenParams(config.param)?.map(p => dataStore.columnIndex[p] as DataColumn<DataType>) : [],
    [config.param, dataStore]);
    // biome-ignore lint/correctness/useExhaustiveDependencies: need to figure out mobx/biome linting...
    const setChartTypeName = useCallback(action((chartTypeName: string) => {
        config.type = chartTypeName;
        const chartType = chartTypes.find(t => t.name === config.type);
        if (!chartType) throw new Error(`Chart type ${config.type} not found`);
        // set default values for params, use previous values if possible
        if (chartType.params) {
            // we could have some metadata on the column like `{intendedFor: "X axis"}` etc
            const previousParams = toJS(config.param) || []; // warning: mobx gives us a proxy & normal indexing doesn't work
            config.param = new Array(chartType.params.length);
            for (const [i, p] of chartType.params.entries()) {
                const previousColumn = paramColumns[i];
                if (previousColumn && previousParams[i] && columnMatchesType(previousColumn, p.type)) {
                    console.log(`[${config.type}, ${p.name}] using previous selection: '${previousParams[i]}'`);
                    config.param[i] = previousParams[i];
                } else {
                    // we could try to be a bit clever about finding a suitable column, like "x"...
                    const c = dataStore.columns.find(c => columnMatchesType(c, p.type));
                    if (!c) {
                        // we should let the user know - highlight the dropdown or something
                        console.warn(`[${config.type}, ${p.name}] No column found for type '${p.type}'`);
                        continue;
                    }
                    console.log(`[${config.type}, ${p.name}] using first available ${p.type}: '${c.name}'`);
                    // todo: figure out why this setting the value in the dropdown (it does apply to the config preview)
                    config.param[i] = c.name;
                }
            }
        } else {
            config.param = [];
        }
        config._updated = new Date();
    }), [config.type, config, chartTypes, paramColumns, dataStore]);
    const chartType = chartTypes.find(t => t.name === config.type);
    const extraControls = useMemo(() => {
        return chartType?.extra_controls?.(dataStore).map((control, i) => (
            <ExtraControlComponent key={control.name} control={control} config={config} />
        ));
    }, [chartType, dataStore, config]);


    const addChart = useCallback(() => {
        const { extra, _updated, ...props } = config;
        // flatten the config into props
        chartType?.extra_controls?.(dataStore).forEach(control => {
            if (control.defaultVal && !extra[control.name]) {
                extra[control.name] = control.defaultVal;
            } else {
                // if (!control.defaultVal)
            }
        });
        props.param = flattenParams(props.param);
        for (const [k, v] of Object.entries(extra)) {
            //@ts-ignore
            props[k] = v;
        }
        // find the key in BaseChart.types that corresponds to this chart type
        // I don't much like this, but it seems to vaguely work for now
        for (const [k, v] of Object.entries(BaseChart.types)) {
            if (v.name === config.type) {
                props.type = k;
                break;
            }
        }
        if (chartType?.init) {
            chartType.init(props, dataStore, extra);
        }

        // this is where we'd call the chart manager to add the chart
        if (!window.mdv.chartManager) throw new Error("chartManager not found");
        const finalConfig = serialiseConfig(props);
        window.mdv.chartManager.addChart(dataStore.name, finalConfig, true);
        // and then close the dialog...
        onDone();
    }, [config, dataStore, onDone, chartType]);
    return (
        <ErrorBoundary FallbackComponent={({error}) =>
            (
            <div className="flex items-center" style={{width: 600, height: 420}}> 
                <ErrorComponentReactWrapper 
                    error={{message: error.message, stack: error.stack}} 
                    extraMetaData={config} 
                    title={chartType?.name 
                        ? `Error creating '${chartType?.name}'. Click to view details`
                        : 'Something went wrong. Click to view details'}
                />
            </div>
            )
        }>
            <DialogContent dividers>
            <Grid container spacing={2} sx={{paddingBottom: 2, width: "100%"}}>
                <Grid size={6}>
                    <Grid container direction={"column"} spacing={1}
                    >
                        <Typography variant="h6">Chart Type</Typography>
                        <Autocomplete
                            options={chartNames}
                            value={config.type}
                            onChange={(_, value) => value && setChartTypeName(value)}
                            renderInput={(params) => (
                                <TextField {...params} label="Chart Type" />
                            )}
                            sx={{pt: 1}}
                        />
                        <TextField label="Title" size="small" onChange={action((e) => config.title = e.target.value)} />
                        <TextField label="Description" size="small"
                        multiline aria-multiline
                        rows={6}
                        onChange={action((e) => config.legend = e.target.value)}
                        />
                    </Grid>
                </Grid>
                <Grid size={6}>
                    <Grid container direction={"column"} spacing={1}
                    sx={{gap: 1}}
                    >
                        {chartType?.params && <Typography variant="h6">Columns</Typography>}
                        {chartType?.params?.map((p, i) => (
                            <>
                                <Typography key={`label_${p.name}`} fontSize="small">
                                    {p.name}:
                                </Typography>
                                <ColumnSelectionComponent key={p.name} placeholder={p.name}
                                multiple={isMultiColumn(p.type)}
                                // this may be changing - perhaps we always pass mobx object to ColumnSelectionComponent
                                // but we definitely need to know whether it's a multi-column or not & be able to set the value accordingly
                                //@ts-expect-error setSelectedColumn(c: never)???
                                setSelectedColumn={action((column) => {
                                    // !! in the original AddChartDialog, we use a "ChooseColumnDialog" for multi-column
                                    // that sets `this.multiColumns` which in `submit` is concatenated to `config.param`
                                    if (!config.param) throw new Error("it shouldn't be possible for config.param to be undefined here");
                                    
                                    // the type of config.param is string[] - but this could be a multi-column or virtual column query...
                                    // we need to decide at which point to apply these transformations.
                                    // - to deal with string[] we should be able to specify multiple={false} which could change the return type
                                    //... could we set up a reaction to update the config when the column changes?

                                    // if (typeof column !== "string") throw new Error("Expected string column name");
                                    config.param[i] = column; //issue now is that because ref is stable, we don't re-flatten the param
                                    // grumble grumble
                                    config._updated = new Date();
                                })}
                                type={p.type}
                                current_value={config.param?.[i]}
                                />
                            </>
                        ))}
                        {extraControls && <h2>Extra Controls</h2>}
                        {extraControls}
                        {/* <Dialog open={false}>
                            <JsonView src={chartType} collapsed/>
                        </Dialog> */}
                    </Grid>
                </Grid>
            </Grid>
            <Divider />
            {/* NOTE - before config._updated was added, this causes update when config.legend chages,
                otherwise not unless we spread config... and then nested `config.params[i] = c` updates didn't work */}
            {/* <ChartPreview config={{...config}} /> */}
            </DialogContent>
            <DialogActions sx={{justifyContent: "center"}}>
                <Button
                    variant="contained"
                    onClick={addChart}
                    disabled={!chartType}
                >
                    Add Chart
                </Button>
            </DialogActions>
        </ErrorBoundary>
    );
});

const AddChartDialogComponent = observer(
    ({ dataStore, config, onDone }: { dataStore: DataStore, config: ChartConfig, onDone: () => void }) => {
        // prompt, e.g.: "dot-plot of gene expression"
        return (
            <>
                <DataStoreContext.Provider value={dataStore}>
                    {/* <div className="align-top"> */}
                    <ConfigureChart config={config} onDone={onDone} />
                    {/* </div> */}
                    {/* NOTE - before config._updated was added, this causes update when config.legend chages,
                    otherwise not unless we spread config... and then nested `config.params[i] = c` updates didn't work */}
                    {/* <div>
                        <p>Chart Preview: "{config.legend}"</p>
                        <JsonView src={{...config}} />
                    </div> */}
                    {/* <ChartPreview config={{...config}} /> */}
                </DataStoreContext.Provider>
            </>
        )
    },
);
/** we could standardise this more for reuse... but not sure we want this to be modal, actually */
const Wrapper = (props: { dataStore: DataStore, modal: boolean, onDone: () => void }) => {
    // how to manage this state?
    // it is currently possible for multiple dialogs to be open at once, with different states
    // this should not be allowed.
    // then we could use zustand without bothering with faff of multiple store contexts
    // the only consequence of allowing multiple dialogs to be open would be that they'd share state
    // but that could actually be bad in that they are all associated with a particular dataStore
    const config = useMemo(() => observable.object((ChartConfigSchema.parse({
            title: "",
            legend: "",
            type: "",
            param: [] as Param[],
            extra: {},
            _updated: new Date(),
        })
    )), []);
    const [open, setOpen] = useState(true);
    useEffect(() => {
        if (!props.modal) return;
        const handleKeyDown = (event: KeyboardEvent) => {
            if (event.key === "Escape") {
                setOpen(false);
                props.onDone();
            }
        };
        window.addEventListener("keydown", handleKeyDown);
        return () => {
            window.removeEventListener("keydown", handleKeyDown);
        };
    }, [props.modal, props.onDone]);
    //p-4 mt-2 z-50 text-center border-2 border-dashed rounded-lg ${isDragOver ? "bg-gray-300 dark:bg-slate-800" : "bg-white dark:bg-black"} min-w-[90%]
    if (!props.modal) {
        return <AddChartDialogComponent {...props} config={config} />;
    }
    return (
        <Dialog open={open} disableEscapeKeyDown={true}
            PaperProps={{
                style: { //copied from FileUploadDialog, should maybe be part of theme.
                //     backgroundColor: "var(--fade_background_color)",
                //     backdropFilter: "blur(1px)",
                border: "1px solid var(--border_menu_bar_color)"
                },
            }}
            fullWidth
            onClose={() => {
                setOpen(false);
                props.onDone();
            }} 
        >
            <DialogTitle>
                Add Chart in "{props.dataStore.name}"
                        <IconButton
                            aria-label="close"
                                            onClick={() => {
                                                setOpen(false);
                                                props.onDone();
                                            }}
                                            sx={{
                                                position: "absolute",
                                                right: 8,
                                                top: 8,
                                            }}
                            >
                            <CloseIcon />
                        </IconButton>
                    </DialogTitle>
                    <AddChartDialogComponent {...props} config={config} />
        </Dialog>
    );
}
class AddChartDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;
    constructor(dataStore: DataStore) {
        super(
            {
                title: `Add Chart in '${dataStore.name}'`,
                width: 600,
                height: 420,
            },
            null,
        );
        const modal = true; //doesn't work on fullscreen panel as of writing, also needs some redesign of "add chart" button in particular
        this.root = createMdvPortal(
            <Wrapper dataStore={dataStore} modal={modal} onDone={() => this.close()}/>,
            this.dialog,
            this,
        );
        this.outer.style.display = modal ? "none" : "block";
    }
    close(): void {
        super.close();
        this.root.unmount();
    }
}
// https://github.com/Taylor-CCB-Group/MDV/discussions/44
// BaseDialog.experiment["AddChartDialogReact"] = AddChartDialogReact;
export default AddChartDialogReact;
