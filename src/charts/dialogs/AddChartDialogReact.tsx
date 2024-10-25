import { useEffect, useMemo, useState } from "react";
import type DataStore from "../../datastore/DataStore.js";
import TagModel from "../../table/TagModel";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import { observer } from "mobx-react-lite";
import TextField from "@mui/material/TextField";
import Autocomplete from "@mui/material/Autocomplete";
import Checkbox from "@mui/material/Checkbox";
import Chip from "@mui/material/Chip";
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";

import BaseChart from "../BaseChart.js";
import { DataStoreContext, useDataStore } from "@/react/context.js";
import JsonView from "react18-json-view";
import { AppBar, Box, Button, Dialog, Grid, Typography } from "@mui/material";
import ColumnSelectionComponent from "@/react/components/ColumnSelectionComponent.js";

function ChooseChartType() {
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
    const [chartTypeName, setChartTypeName] = useState(chartNames[0]);
    const chartType = chartTypes.find(t => t.name === chartTypeName);
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
            <Grid container spacing={2} sx={{margin: 1}}>
                <Grid item xs={6}>
                    <Grid container direction={"column"} spacing={1}
                    sx={{gap: 1}}
                    >
                        <h2>Chart Type</h2>
                        <Autocomplete
                            options={chartNames}
                            value={chartTypeName}
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
}

const AddChartDialogComponent = observer(
    ({ dataStore }: { dataStore: DataStore }) => {
        const chartTypes = BaseChart.types;
        const [selectedChartType, setSelectedChartType] = useState("");
        return (
            <>
                <DataStoreContext.Provider value={dataStore}>
                    <div className="align-top">
                    <ChooseChartType />
                    </div>
                </DataStoreContext.Provider>
            </>
        )
    },
);

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
        this.root = createMdvPortal(
            <AddChartDialogComponent dataStore={dataStore} />,
            this.dialog,
            this,
        );
    }
    close(): void {
        super.close();
        this.root.unmount();
    }
}
// https://github.com/Taylor-CCB-Group/MDV/discussions/44
BaseDialog.experiment["AddChartDialogReact"] = AddChartDialogReact;
export default "AddChartDialogReact loaded";
