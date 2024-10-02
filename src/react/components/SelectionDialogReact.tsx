import BaseChart from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import SelectionDialogComponent from "./SelectionDialogComponent";

type SelectionDialogConfig = {
    type: "selection_dialog_experimental";
};

class SelectionDialogReact extends BaseReactChart<SelectionDialogConfig> {
    constructor(dataStore, div, config) {
        super(dataStore, div, config, SelectionDialogComponent);
        console.log("todo add reset button & associated logic in SelectionDialogReact");
    }
}

BaseChart.types["selection_dialog_experimental"] = {
    name: "Selection Dialog (Experimental)",
    class: SelectionDialogReact,
    params: [
        {
            type: "_multi_column:all",
            name: "Columns To filter",
        },
    ],
}