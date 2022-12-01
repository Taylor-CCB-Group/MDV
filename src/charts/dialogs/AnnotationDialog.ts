import { DataModel } from "../../table/DataModel.js";
import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";
import { DataColumn, DataSource } from "../charts.js";


export default class AnnotationDialog extends BaseDialog {
    tagColumn: DataColumn<'text'>;
    constructor(dataSource: DataSource) {
        super({
            title: "Annotate selection"
        }, null);
        // is there a '__tags' column? if so, use it, otherwise add it.
        // do we actually want/need a DataModel here?
        const {dataStore} = dataSource;
        const dataModel = new DataModel(dataStore, {autoupdate: false});
        if (!dataStore.columnIndex['__tags']) dataModel.createColumn('__tags', null);
        this.tagColumn = dataStore.columnIndex['__tags'];
    }
}