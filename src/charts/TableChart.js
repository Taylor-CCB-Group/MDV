import { SlickGrid, TextEditor, Overlays } from "../table/SlickGrid.js";
import { RowSelectionModel } from "../table/RowSelectionModel.js";
import { DataModel } from "../table/DataModel.js";
import BaseChart from "./BaseChart";
import { createEl } from "../utilities/Elements.js";
import { BaseDialog } from "../utilities/Dialog.js";

class TableChart extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        this.config.type = "table_chart";
        let cols = [];
        // react to changes in config.param
        this.mobxAutorun(() => {
            const { config, dataStore } = this;
            cols = [];
            //add the index column (only if include_index is true)
            if (config?.include_index)
                cols = [
                    {
                        field: "__index__",
                        id: "__index__",
                        name: "index",
                        datatype: "integer",
                        sortable: true,
                        width: 100,
                    },
                ];
            const cw = config.column_widths || {};
            const oldCols = this.grid?.getColumns() || [];
            oldCols.forEach(c => {
                cw[c.field] = c.width;
            });
            for (const c of config.param) {
                const column = dataStore.columnIndex[c];
                const col = {
                    id: column.field,
                    field: column.field,
                    name: column.name,
                    datatype: column.datatype,
                    sortable: true,
                    width: cw[c] ? cw[c] : 100,
                };
                // TODO review PJT coerce multitext internally... for now(?)
                if (column.datatype === "multitext") col.datatype = "text";

                if (column.editable) {
                    col.editor = TextEditor;
                }
                if (column.is_url) {
                    col.formatter = (row, cell, value, columnDef, dataContext) => {
                        return `<a href="${value}" target="_blank">${value}</a>`;
                    };
                }
                cols.push(col);
            }
            window.mdv.chartManager.loadColumnSet(config.param, dataStore.name, () => {
                this.grid?.setColumns(cols);
                this.dataModel?.setColumns(config.param);
                this.grid?.render();
                try {
                    this.onDataFiltered();
                } catch {}
            })
        });
        this.options = {
            enableCellNavigation: true,
            enableColumnReorder: false,
            editable: false,
            autoEdit: true,
            enableAsyncPostRender: true,
            frozenColumn: null,
            //showHeaderRow:true,
            //headerRowHeight:40
        };
        this.dataModel = new DataModel(dataStore, { autoupdate: false });
        this.dataModel.setColumns(this.config.param);
        this.grid = new SlickGrid(
            this.contentDiv,
            this.dataModel,
            cols,
            this.options,
        );
        this.grid.onHeaderCellRendered.subscribe((e, args) => {
            const c = args.column;
            if (c.editor && this.mode === "edit") {
                const i = createEl("i", {
                    classes: ["far", "fa-edit"],
                    style: { cursor: "pointer" },
                });
                args.node.prepend(i);
                i.addEventListener("click", (e) => {
                    e.stopImmediatePropagation();
                    new EditColumnDialog(this, c.field, c.name, this.__doc__);
                });
            }
        });
        this.dataModel.addListener(this.config.id, () => {
            this.grid.invalidate();
        });
        //this is to update the header renderer
        this.grid.setColumns(this.grid.getColumns());

        this.mode = "select";
        //this._changeMode();
        this.grid.setSelectionModel(new RowSelectionModel());

        this.grid.init();

        this.grid.onSort.subscribe((e, args) => {
            this._sort(args);
        });

        this.grid.onSelectedRowsChanged.subscribe((e, args) => {
            if (this.mode === "select") {
                this._rowsSelected(args);
            }
        });

        this.editIcon = this.addMenuIcon("fas fa-edit", "Change To Edit Mode", {
            func: () => {
                this.editIcon.firstElementChild.remove();
                this.mode = this.mode === "select" ? "edit" : "select";
                const classes =
                    this.mode === "select"
                        ? ["fas", "fa-edit"]
                        : ["fas", "fa-table"];
                createEl("i", { classes: classes }, this.editIcon);
                this.editIcon.setAttribute(
                    "aria-label",
                    this.mode === "edit"
                        ? "Change To Selection Mode"
                        : "Change To Edit Mode",
                );
                this._changeMode();
            },
        });

        this.addMenuIcon("fas fa-download", "Download data", {
            func: () => {
                this.downloadData();
            },
        });

        this.filter = [];
        this.themeChanged();
        this.onDataFiltered();
        if (this.config.sort) this._sort(this.config.sort); //now done in onDataFiltered()
    }
    _sort({ columnId, sortAsc }) {
        this.dataModel.sort(columnId, sortAsc ? "asc" : "desc");
        this.grid.invalidateAllRows();
        this.grid.render();
        this.config.sort = { columnId, sortAsc };
    }

    onDataHighlighted(data) {
        if (data.source === this) {
            return;
        }
        
        // Get all highlighted indices
        let indices = [];

        if (data.indexes) {
            indices = Array.isArray(data.indexes)
                ? data.indexes
                : Object.values(data.indexes);
        }
        
        if (indices.length === 0) {
            this.grid.setSelectedRows([]);
            return;
        }

        // Find all positions in the data model
        const positions = [];
        for (const index of indices) {
            const pos = this.dataModel.data.indexOf(index);
            if (pos !== -1) {
                positions.push(pos);
            }
        }

        // Clear the selection if there are no selected rows
        if (positions.length === 0) {
            this.tempMode = this.mode;
            this.mode = "";
            this.grid.setSelectedRows([]);
            this.mode = this.tempMode;
            return;
        }

        if (positions.length > 0) {
            // Scroll to the first highlighted row
            this.grid.scrollRowIntoView(positions[0]);
            this.tempMode = this.mode;
            this.mode = "";
            // Set all highlighted rows
            this.grid.setSelectedRows(positions);
            this.mode = this.tempMode;
        }
    }

    onColumnRemoved(column) {
        if (this.config.param.indexOf(column) === -1) {
            return false;
        }
        const editor = this.grid.getCellEditor();
        if (editor) {
            editor.destroy();
        }
        this.grid.setActiveCell(null);
        let cols = this.grid.getColumns();
        cols = cols.filter((x) => x.field !== column);
        this.config.param = this.config.param.filter((x) => x !== column);
        this.dataModel.setColumns(this.config.param);
        this.grid.setColumns(cols);
        return false;
    }

    getConfig() {
        const config = super.getConfig();
        const cols = this.grid.getColumns();

        config.param = cols
            .filter((x) => x.field !== "__index__")
            .map((x) => x.field);
        config.column_widths = {};
        for (const c of cols) {
            if (c.width !== 100) {
                config.column_widths[c.field] = c.width;
            }
        }
        return config;
    }

    async downloadData() {
        const blob = await this.dataModel.getDataAsBlob();
        const save = createEl(
            "a",
            {
                download: this.dataStore.name,
                target: "_blank",
                href: window.URL.createObjectURL(blob),
            },
            this.__doc__.body,
        );
        save.click();
        save.remove();
    }

    _changeMode() {
        if (this.mode === "select") {
            if (this.addColumnIcon) {
                this.addColumnIcon.remove();
                this.addColumnIcon = null;
            }
            this.selectionModel = new RowSelectionModel();
            this.grid.setSelectionModel(this.selectionModel);
            if (this.overlay) {
                this.grid.unregisterPlugin(this.overlay);
                this.overlay = null;
            }

            this.grid.setOptions({ editable: false }, true, true);
            this.grid.setColumns(this.grid.getColumns());
        } else {
            this.grid.setSelectedRows([]);
            if (this.selectionModel) {
                this.grid.setSelectionModel(null);
                this.selectionModel = null;
            }

            this.grid.setOptions({ editable: true }, true, true);
            this.grid.setColumns(this.grid.getColumns());
            this.addColumnIcon = this.addMenuIcon("fas fa-plus", "Add Column", {
                func: (e) => {
                    new AddColumnDialog(this);
                },
            });
            this.overlay = new Overlays({});
            this.overlay.onFillUpDown.subscribe((e, args) =>
                this._updateOverlay(args),
            );

            this.grid.registerPlugin(this.overlay);
        }
    }

    //all we need is to redraw the chart
    onDataChanged(data) {
        this.grid.invalidate();
    }

    getColumnInfo() {
        return this.grid
            .getColumns()
            .filter((x) => {
                //PJT suspected "multitext" needed to be handled here,
                //now working on the basis that any 'multitext' entering this system will be described as 'text'
                return x.datatype === "text";
            })
            .map((x) => {
                return { name: x.name, field: x.field };
            })
            .sort((a, b) => a.name.localeCompare(b.name));
    }

    createColumn(text, cloneColumn = null, position = 2) {
        this.dataModel.createColumn(text, cloneColumn);
        const cols = this.grid.getColumns();
        position -= 1;
        if (position < 0) {
            position === 0;
        } else if (position >= cols.length) {
            position = cols.length - 1;
        }
        cols.splice(position, 0, {
            name: text,
            field: text,
            id: text,
            datatype: "text",
            editor: TextEditor,
            sortable: true,
        });
        this.config.param.push(text);
        this.grid.setColumns(cols);
        this.grid.invalidate();
    }

    _removeColumn(col) {
        //remove column and set as dirty plus
        this.dataModel.removeColumn(col, true, true);
    }

    _updateOverlay(args) {
        const column = this.grid.getColumns()[args.range.fromCell];

        // Ensure the column is editable
        if (!column.editor) {
            return;
        }

        this.dataModel.updateRange(
            column.field,
            args.range.fromRow,
            args.range.toRow,
            args.dir,
        );
        this.grid.invalidate();
    }
    _rowsSelected(args) {
        const indexes = [];
        for (const i of args.rows) {
            indexes.push(this.dataModel.data[i]);
        }

        this.dataStore.dataHighlighted(indexes, this);
    }

    // pjt: deprecated?
    themeChanged() {
        // console.warn('themeChanged() deprecated');
    }

    changeBaseDocument(doc) {
        super.changeBaseDocument(doc);
        this.grid.setBaseDocument(doc);
        if (this.mode === "edit") {
            this.grid.unregisterPlugin(this.overlay);
            this.grid.registerPlugin(this.overlay);
        }
    }

    setSize(x, y) {
        super.setSize(x, y);
        this.grid.resizeCanvas();
        this.grid.invalidate();
    }

    pinChart() {
        this.isPinned = true;
    }

    unpinChart() {
        this.isPinned = false;
        this.onDataFiltered();
    }

    addColumns() {}
    remove() {
        this.grid.destroy();
        super.remove();
    }

    onDataFiltered() {
        if (this.isPinned) {
            return;
        }
        this.dataModel.updateModel();
        if (this.config.sort) {
            this._sort(this.config.sort);
        }
    }
}

class EditColumnDialog extends BaseDialog {
    constructor(table, column, columnName, doc) {
        const config = {
            doc: doc,
            columns: 2,
            width: 380,
            maxHeight: 500,
            title: `Bulk Edit ${columnName}`,
        };
        super(config, { table, column });
    }
    init(content) {
        this.col = content.column;
        this.dataModel = content.table.dataModel;
        this.table = content.table;
        const d1 = createEl("div", {}, this.columns[0]);
        createEl("label", { text: "Value:" }, d1); //a11y: label for?

        this.valInput = createEl("input", {}, d1);

        this.valInput.addEventListener("keyup", (e) => {
            const k = e.which;
            if (k === 37 || k === 39 || k === 8) {
                e.stopImmediatePropagation();
                return;
            }
            const suggest = this.dataModel.getValueSuggestion(
                this.valInput.value,
                this.col,
            );
            if (suggest) {
                this.valInput.value = suggest;
            }
        });

        const d2 = createEl("div", {}, this.columns[0]);
        createEl("label", { text: "Replace With:" }, d2); //a11y: label for?

        this.replaceInput = createEl("input", {}, d2);

        this._createButton("Fill All Cells", "_all_");
        this._createButton("Fill Empty Cells", "_blank_");
        this._createButton("Replace Value", null);
        this._createButton("Delete Unused Values", "_delete_");
        this._createButton("Remove Column", "_delete_column_");
    }

    _createButton(text, replace) {
        const d = createEl(
            "div",
            { style: { padding: "5px 10px" } },
            this.columns[1],
        );
        createEl(
            "button",
            { classes: ["ciview-button"], text: text },
            d,
        ).addEventListener("click", () => {
            const r = replace || this.replaceInput.value;
            if (r === "_delete_column_") {
                this.table._removeColumn(this.col);
                this.close();
            } else {
                this.dataModel.replaceValues(this.valInput.value, r, this.col);
            }
        });
    }
}

class AddColumnDialog extends BaseDialog {
    constructor(table, doc) {
        const config = {
            doc: doc,
            footer: true,
            buttons: [
                {
                    text: "Add",
                    method: "addColumn",
                },
            ],
            width: 300,
            maxHeight: 500,
            title: "Add Column",
        };
        super(config, table);
    }
    init(table) {
        const dstyle = {
            padding: "5px",
        };
        this.table = table;
        const d1 = createEl("div", { style: dstyle }, this.dialog);
        createEl("label", { text: "Column Name" }, d1);
        this.name = createEl("input", {}, d1);
        this.name.addEventListener("keydown", (e) => {
            if (e.key === "Enter") this.addColumn();
        });
        const d2 = createEl("div", { style: dstyle }, this.dialog);

        const dc = createEl("label", { style: { display: "flex" } }, d2); //a11y: not convinced about label arrangement here
        this.cloneCheck = createEl("input", { type: "checkbox" }, dc);
        createEl("div", { text: "Copy Existing Column" }, dc);
        const cols = table.getColumnInfo();
        this.cloneCol = createEl("select", {}, d2);
        for (const c of cols) {
            createEl("option", { text: c.name, value: c.field }, this.cloneCol);
        }
        const d3 = createEl("div", { style: dstyle }, this.dialog);
        createEl("label", { text: "Position" }, d3);
        this.position = createEl(
            "input",
            { style: { width: "100px" }, value: 2 },
            d3,
        );

        setTimeout(() => this.name.focus(), 100);
    }

    addColumn() {
        if (!this.name.value) {
            this.name.focus();
            return;
        }
        const clone = this.cloneCheck.checked ? this.cloneCol.value : null;
        let pos = Number.parseInt(this.position.value);
        pos = Number.isNaN(pos) ? null : pos;
        this.table.createColumn(this.name.value, clone, pos);
        this.close();
    }
}

export default TableChart;

BaseChart.types["table_chart"] = {
    class: TableChart,
    name: "Table",
    params: [
        {
            type: "_multi_column:all",
            name: "Columns To Display",
        },
    ],
    extra_controls: () => {
        //maybe '_multi_column:all' should handle this?
        return [
            {
                type: "checkbox",
                name: "include_index",
                label: "Include Index",
            },
        ];
    },
    init: (config, ds, ec) => {
        config.include_index = ec.include_index;
    },
};
