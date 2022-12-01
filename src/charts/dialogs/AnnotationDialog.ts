import { DataModel } from "../../table/DataModel.js";
import { setTagOnAllValues } from "../../table/DataTagOperations";
import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";
import { DataColumn, DataSource } from "../charts.js";


export default class AnnotationDialog extends BaseDialog {
    tagColumn: DataColumn<'text'>;
    constructor(dataSource: DataSource) {
        super({
            title: "Annotate selection",
            width: 400,
            height: 200,
            columns: 2
        }, null);
        // is there a '__tags' column? if so, use it, otherwise add it.
        // do we actually want/need a DataModel here?
        const {dataStore} = dataSource;
        const dataModel = new DataModel(dataStore, {autoupdate: true});
        if (!dataStore.columnIndex['__tags']) dataModel.createColumn('__tags', null);
        this.tagColumn = dataStore.columnIndex['__tags'];

        const tagInput = createEl("input", {}, this.columns[0]); //chrome treats as password field (but without hidden value)?
        const addTagButton = createEl("button", {text: "Add tag to selection"}, this.columns[1]);
        addTagButton.addEventListener("click", () => {
            setTagOnAllValues(tagInput.value, this.tagColumn, dataModel);
        });
    }
}