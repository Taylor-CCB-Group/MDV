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
import ColumnSelectionComponent from "@/react/components/ColumnSelectionComponent.js";
import { makeAutoObservable, runInAction } from "mobx";
import z from "zod";

const ChartConfigSchema = z.object({
    title: z.string(),
    legend: z.string(),
    // in future, allow for "virtual" / "computed" / "smart" columns
    param: z.array(z.string()),
    type: z.string(),
    // in the original AddChartDialog, extra props are on the root object, not nested
    // so when we pass this to ChartManager, we'll need to flatten it
    extra: z.record(z.unknown()),
});
type ChartConfig = z.infer<typeof ChartConfigSchema>;

const ConfigureChart = observer(({config}: {config: ChartConfig}) => {
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
    const setChartTypeName = useCallback((chartTypeName: string) => {
        runInAction(() => {
            config.type = chartTypeName;
        });
    }, [config.type, config]);
    const chartType = chartTypes.find(t => t.name === config.type);
    const extraControls = useMemo(() => {
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
                        <TextField label="Title" size="small" />
                        <TextField label="Description" size="small" 
                        multiline aria-multiline
                        rows={6}
                        />
                    </Grid>
                </Grid>
                <Grid item xs={6}>
                    <Grid container direction={"column"} spacing={1}
                    sx={{gap: 1}}
                    >
                        <h2>Columns</h2>
                        {chartType?.params?.map(p => (
                            <ColumnSelectionComponent key={p.name} placeholder={p.name} 
                            setSelectedColumn={(column) => console.log(column)}
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
                <Button variant="contained" color="primary" disabled>Add Chart</Button>
            </AppBar>
        </>
    );
});

const AddChartDialogComponent = observer(
    ({ dataStore, config }: { dataStore: DataStore, config: ChartConfig }) => {
        const chartTypes = BaseChart.types;
        const [selectedChartType, setSelectedChartType] = useState("");
        return (
            <>
                <DataStoreContext.Provider value={dataStore}>
                    <div className="align-top">
                    <ConfigureChart config={config} />
                    </div>
                    <div className="w-full h-[40vh]">
                        <p>Chart Preview: "{config.title}"</p>
                        <JsonView src={config} />
                    </div>
                </DataStoreContext.Provider>
            </>
        )
    },
);
/** we could standardise this more for reuse... but not sure we want this to be modal, actually */
const Wrapper = (props: { dataStore: DataStore, modal: boolean }) => {
    // how to manage this state?
    // it is currently possible for multiple dialogs to be open at once, with different states
    // this should not be allowed.
    // then we could use zustand without bothering with faff of multiple store contexts
    // the only consequence of allowing multiple dialogs to be open would be that they'd share state
    // but that could actually be bad in that they are all associated with a particular dataStore
    const config = useMemo(() => makeAutoObservable(ChartConfigSchema.parse({
            title: "",
            legend: "",
            type: "",
            param: [] as Param[],
            extra: {},
        } satisfies ChartConfig)
    ), []);
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
            <Wrapper dataStore={dataStore} modal={modal} />,
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
