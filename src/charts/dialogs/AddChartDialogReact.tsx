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

function ChooseChartType() {
    const [chartTypeName, setChartTypeName] = useState("");
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
    const chartNames = chartTypes.map(t => t.name);
    const chartType = chartTypes.find(t => t.name === chartTypeName);
    
    return (
        <>
            <h2>Chart Type</h2>
            <Autocomplete
                options={chartNames}
                value={chartTypeName}
                onChange={(_, value) => setChartTypeName(value)}
                renderInput={(params) => (
                    <TextField {...params} label="Chart Type" />
                )}
            />
            <JsonView src={chartType} />
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
                    <ChooseChartType />
                </DataStoreContext.Provider>
            </>
        )
    },
);

class AddChartDialogReact extends BaseDialog {
    // tagModel: TagModel; //prefer to keep this state in react... but we do need to know what the dataStore is.
    // tagColumn: DataColumn<'text'>;
    // dataModel: DataModel;
    // tagListElement: HTMLDivElement;
    // tagInput: any;
    root: ReturnType<typeof createMdvPortal>;
    constructor(dataStore: DataStore) {
        super(
            {
                title: `Add Chart in '${dataStore.name}'`,
                width: 400,
                height: 300,
            },
            null,
        );
        // this.outer.classList.add('annotationDialog');
        // this.tagModel = new TagModel(dataStore);
        // this.tagModel = tagModel;
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
