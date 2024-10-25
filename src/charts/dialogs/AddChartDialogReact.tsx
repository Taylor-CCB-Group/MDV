import { useCallback, useEffect, useMemo, useState } from "react";
import type DataStore from "../../datastore/DataStore.js";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import { observer } from "mobx-react-lite";
import TextField from "@mui/material/TextField";
import Autocomplete from "@mui/material/Autocomplete";
import Chip from "@mui/material/Chip";

import BaseChart from "../BaseChart.js";
import type { Param } from "@/charts/ChartTypes.js";
import { DataStoreContext, useDataStore } from "@/react/context.js";
import JsonView from "react18-json-view";
import { AppBar, Box, Button, Dialog, Grid, Paper, Typography } from "@mui/material";
import ColumnSelectionComponent, { columnMatchesType } from "@/react/components/ColumnSelectionComponent.js";
import { action, observable, runInAction } from "mobx";
import z from "zod";
import type { DataColumn } from "../charts.js";

const ChartConfigSchema = z.object({
    title: z.string(),
    legend: z.string(),
    // in future, allow for "virtual" / "computed" / "smart" columns
    param: z.array(z.string()),
    type: z.string(),
    // in the original AddChartDialog, extra props are on the root object, not nested
    // so when we pass this to ChartManager, we'll need to flatten it
    extra: z.record(z.unknown()),
    _updated: z.date(),
});
type ChartConfig = z.infer<typeof ChartConfigSchema>;

// may want something like this... or to obviate the need for it, ideally.
function actionWithUpdate<T extends any[]>(
    config: ChartConfig,
    fn: (...args: T) => void
): (...args: T) => void {
    return action((...args: T): void => {
        fn(...args);
        config._updated = new Date();
    });
}

const ChartPreview = observer(({config}: {config: ChartConfig}) => {
    // this could be an actual chart preview, but for now, the config is useful
    return (
        <div>
            {/* <h2>Chart Preview</h2> */}
            <JsonView src={config} collapsed />
        </div>
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
                return t.required.every(r => !!dataStore[r]);
            }
            return true;
        });
    }, [dataStore]);
    const chartNames = chartTypes.map(t => t.name).sort((a, b) => a.localeCompare(b));
    // const [chartTypeName, setChartTypeName] = useState(chartNames[0]);
    const paramColumns = useMemo(
        () => config.param.map(p => dataStore.columnIndex[p] as DataColumn),
    [config.param, dataStore]);
    const setChartTypeName = useCallback((chartTypeName: string) => {
        runInAction(() => {
            config.type = chartTypeName;
            const chartType = chartTypes.find(t => t.name === config.type);
            if (!chartType) throw new Error(`Chart type ${config.type} not found`);
            // set default values for params, use previous values if possible
            const previousParams = config.param;
            // what is the observable status of this new array?
            config.param = new Array(chartType.params.length);
            for (const [i, p] of chartType.params.entries()) {
                const previousColumn = paramColumns[i];
                if (previousParams[i] && columnMatchesType(previousColumn, p.type)) {
                    config.param[i] = previousParams[i];
                } else {
                    // we could try to be a bit clever about finding a suitable column, like "x"...
                    config.param[i] = "";
                }
            }
            config._updated = new Date();
        });
    }, [config.type, config, chartTypes, paramColumns]);
    const chartType = chartTypes.find(t => t.name === config.type);
    const extraControls = useMemo(() => {
        // I don't think these need to be observer components
        // (and anyway, this part of GUI is totally placeholder for now)
        return chartType?.extra_controls?.(dataStore).map((control, i) => (
            <Box key={control.name}>
                <Typography>{control.label}</Typography>
                <p>
                {control.name} :  {control.type}
                </p>
                <p>
                defaultVal: {control.defaultVal}
                </p>
                <div>
                    <h3>Values</h3>
                    {control.values?.map(({name, value}) => (
                        <Chip key={name} label={name} />
                    ))}
                </div>
            </Box>
        ));
    }, [chartType, dataStore]);

    const addChart = useCallback(() => {
        const { extra, _updated, ...props } = config;
        // flatten the config into props
        for (const [k, v] of Object.entries(extra)) {
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

        // this is where we'd call the chart manager to add the chart
        window.mdv.chartManager.addChart(dataStore.name, props, true);
        // and then close the dialog...
        onDone();
    }, [config, dataStore]);
    return (
        <>
            <Grid container spacing={2} sx={{margin: 1, width: '50em'}}>
                <Grid item xs={6}>
                    <Grid container direction={"column"} spacing={1}
                    sx={{gap: 1}}
                    >
                        <h2>Chart Type</h2>
                        <Autocomplete
                            options={chartNames}
                            value={config.type}
                            onChange={(_, value) => setChartTypeName(value)}
                            renderInput={(params) => (
                                <TextField {...params} label="Chart Type" />
                            )}
                        />
                        <TextField label="Title" size="small" onChange={action((e) => config.title = e.target.value)} />
                        <TextField label="Description" size="small" 
                        multiline aria-multiline
                        rows={6}
                        onChange={action((e) => config.legend = e.target.value)}
                        />
                    </Grid>
                </Grid>
                <Grid item xs={6}>
                    <Grid container direction={"column"} spacing={1}
                    sx={{gap: 1}}
                    >
                        <h2>Columns</h2>
                        {chartType?.params?.map((p, i) => (
                            <ColumnSelectionComponent key={p.name} placeholder={p.name} 
                            setSelectedColumn={action((column) => {
                                config.param[i] = column;
                                // grumble grumble
                                config._updated = new Date();
                            })}
                            type={p.type}
                            />
                        ))}
                        {extraControls}
                        <Dialog open={false}>
                            <JsonView src={chartType} collapsed/>
                        </Dialog>
                    </Grid>
                </Grid>
            </Grid>
            <AppBar position="absolute" sx={{ top: 'auto', bottom: 0 }}>
                <Button variant="contained" color="primary" onClick={addChart}>Add Chart</Button>
            </AppBar>
        </>
    );
});

const AddChartDialogComponent = observer(
    ({ dataStore, config, onDone }: { dataStore: DataStore, config: ChartConfig, onDone: () => void }) => {
        // prompt, e.g.: "dot-plot of gene expression"
        return (
            <>
                <DataStoreContext.Provider value={dataStore}>
                    <div className="align-top">
                    <ConfigureChart config={config} onDone={onDone} />
                    </div>
                    {/* NOTE - before config._updated was added, this causes update when config.legend chages,
                    otherwise not unless we spread config... and then nested `config.params[i] = c` updates didn't work */}
                    {/* <div>
                        <p>Chart Preview: "{config.legend}"</p>
                        <JsonView src={{...config}} />
                    </div> */}
                    <ChartPreview config={{...config}} />
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
            }
        };
        window.addEventListener("keydown", handleKeyDown);
        return () => {
            window.removeEventListener("keydown", handleKeyDown);
        };
    }, [props.modal]);
    //p-4 mt-2 z-50 text-center border-2 border-dashed rounded-lg ${isDragOver ? "bg-gray-300 dark:bg-slate-800" : "bg-white dark:bg-black"} min-w-[90%]
    if (!props.modal) {
        return <AddChartDialogComponent {...props} config={config} />;
    }
    return (
        <Dialog open={open} fullScreen disableEscapeKeyDown={true}>
            <div className="h-screen flex items-center justify-center">
                <Paper elevation={2} sx={{ p: 2 }}>
                    <h1>Add Chart in "{props.dataStore.name}"...</h1>
                    <AddChartDialogComponent {...props} config={config} />
                </Paper>
            </div>
        </Dialog>
    );
}
class AddChartDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;
    constructor(dataStore: DataStore) {
        super(
            {
                title: `Add Chart in '${dataStore.name}'`,
                width: 400,
                height: 420,
            },
            null,
        );
        const modal = true;
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
BaseDialog.experiment["AddChartDialogReact"] = AddChartDialogReact;
export default "AddChartDialogReact loaded";
