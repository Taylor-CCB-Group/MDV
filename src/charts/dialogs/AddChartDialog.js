import { createEl, createFilterElement } from "../../utilities/Elements";
import BaseChart from "../BaseChart";
import { BaseDialog } from "../../utilities/Dialog.js";
import { ChooseColumnDialog } from "./ChooseColumnDialog";

/**
 * Creates a dialog for the user to choose a chart and its associated parameters. When chosen the
 * supplied callback will be invoked with the config of the chosen chart.
 * @param {DataStore} dataStore - the dataStore the chart will be created from.
 * @param {function} callback - A function called when the user has selected the chart and
 * its parameters. The callback is provided with the config of the chosen chart
 */
export class AddChartDialog extends BaseDialog {
    constructor(dataSource, callback) {
        const config = {
            title: "Add Chart",
            columns: 2,
            footer: true,
            width: 380,
        };
        super(config, { dataSource: dataSource, callback: callback });
    }
    init(content) {
        this.extraControls = {};
        const types = [];
        this.dataSource = content.dataSource;
        this.dataStore = content.dataSource.dataStore;
        for (const type in BaseChart.types) {
            const t = BaseChart.types[type];
            //check to see if chart has any requirements
            if (t.required) {
                if (typeof t.required === "function") {
                    if (!t.required(this.dataStore)) {
                        continue;
                    }
                }

                //is an array of parameters required in the datasource
                else {
                    let allow = true;
                    for (const r of t.required) {
                        if (!this.dataStore[r]) {
                            allow = false;
                        }
                    }
                    if (!allow) {
                        continue;
                    }
                }
            }
            if (t.allow_user_add === false) {
                continue;
            }
            types.push({
                name: t.name,
                type: type,
            });
        }

        types.sort((a, b) => a.name.localeCompare(b.name));
        this.defaultType = types[0].type;

        createEl(
            "div",
            {
                text: "Chart Type", //a11y - should be a label
                classes: ["ciview-title-div"],
            },
            this.columns[0]
        );

        this.chartType = createEl(
            "select",
            {
                styles: {
                    maxWidth: "200px",
                },
            },
            this.columns[0]
        );
        for (const item of types) {
            createEl(
                "option",
                {
                    text: item.name,
                    value: item.type,
                },
                this.chartType
            );
        }
        // createEl("div",{},this.columns[0]).append(this.chartType);
        this.chartType.addEventListener("change", (e) => {
            this.setParamDiv(this.chartType.value, content.dataStore);
        });

        createEl(
            "div",
            {
                text: "Title",
                classes: ["ciview-title-div"],
            },
            this.columns[0]
        );

        this.chartName = createEl("input", {}, this.columns[0]);

        createEl(
            "div",
            {
                text: "Description",
                classes: ["ciview-title-div"],
            },
            this.columns[0]
        );
        this.chartDescription = createEl(
            "textarea",
            { styles: { height: "100px" } },
            this.columns[0]
        );

        this.columnsHeading = createEl(
            "div",
            {
                text: "Columns",
                classes: ["ciview-title-div"],
            },
            this.columns[1]
        );
        this.paramDiv = createEl("div", {}, this.columns[1]);
        this.setParamDiv(types[0].type, content.dataStore);

        createEl(
            "button",
            {
                text: "Add",
                classes: ["ciview-button"],
            },
            this.footer
        ).addEventListener("click", () => this.submit(content.callback));
    }

    submit(callback) {
        // this logic to be used as a reference in React version (at least initially)
        const config = {
            title: this.chartName.value,
            legend: this.chartDescription.value,
            type: this.chartType.value,
            param: this.paramSelects.map((x) => x.value),
            // options: this.options ? Object.fromEntries(this.options) : undefined,
        };
        const ed = {};
        for (const name in this.extraControls) {
            const c = this.extraControls[name];
            ed[name] = c.type === "checkbox" ? c.checked : c.value;
            //mjs: sometimes its not as simple as this and more complex alterations
            //to the config are required based on the user input the dataSore's config
            config[name] = ed[name]; // pjt: is there a reason we didn't do this before?
        }
        if (this.multiColumns) {
            config.param = config.param.concat(this.multiColumns);
        }
        console.log("config from add chart dialog", config);
        const t = BaseChart.types[this.chartType.value];

        if (t.init) {
            t.init(config, this.dataSource.dataStore, ed);
        }
        callback(config);
        this.chartName.value = "";
        this.chartDescription.value = "";
        /// pjt I find this annoying... not sure why we didn't close the div before
        /// but otherwise, would rather not reset these (can be handy when testing stuff)
        // this.chartType.value= this.defaultType;
        // this.setParamDiv(this.defaultType)
        this.close();
    }

    _addMultiColumnSelect(holder, filter) {
        //get default values
        const ps = this.dataStore.getColumnList(filter);
        let text = "";
        if (ps.length > 1) {
            text = `${ps[0].name},... (1)`;
            this.multiColumns = [ps[0].field];
        }
        const dd = createEl("span", { text: text }, holder);
        createEl("i", { classes: ["fas", "fa-plus"] }, holder);
        holder.style.cursor = "pointer";
        holder.addEventListener("click", () => {
            new ChooseColumnDialog(
                this.dataStore,
                (cols) => {
                    this.multiColumns = cols;
                    let text = "";
                    if (cols.length > 0) {
                        const max = cols.length < 3 ? cols.length : 3;
                        const arr = [];
                        for (let n = 0; n < max; n++) {
                            arr.push(this.dataStore.getColumnName(cols[n]));
                        }
                        text = arr.join(",");
                        if (cols.length > 3) {
                            text += ",....";
                        }
                        text += `(${cols.length})`;
                        dd.textContent = text;
                    }
                },
                filter
            );
        });
    }

    setParamDiv(type) {
        this.paramDiv.innerHTML = "";
        const params = BaseChart.types[type].params;
        this.paramSelects = [];
        this.columnsHeading.style.display = params?.length ? "" : "none";
        if (params) {
            for (const p of params) {
                const d = createEl(
                    "div",
                    { styles: { padding: "4px" } },
                    this.paramDiv
                );
                const sp = createEl("div", { text: `${p.name}:` }, d);
                const holder = createEl("div", {}, this.paramDiv);
                if (!Array.isArray(p.type) && p.type.startsWith("_multi")) {
                    this._addMultiColumnSelect(holder, p.type.split(":")[1]);
                } else {
                    this.multiColumns = null;
                    const dd = createEl(
                        "select",
                        {
                            styles: {
                                maxWidth: "200px",
                            },
                        },
                        holder
                    );
                    const ps = this.dataStore.getColumnList(p.type);
                    const sgs = {};
                    for (const ds of this.dataStore.subgroupDataSources) {
                        sgs[ds] = createEl("optgroup", { label: ds });
                    }
                    for (const item of ps) {
                        const ele = item.subgroup
                            ? sgs[item.subgroup.dataSource]
                            : dd;
                        createEl(
                            "option",
                            { text: item.name, value: item.field },
                            ele
                        );
                    }
                    for (const ds of this.dataStore.subgroupDataSources) {
                        dd.append(sgs[ds]);
                    }
                    createFilterElement(dd, holder);
                    this.paramSelects.push(dd);
                }
            }
        }
        const t = BaseChart.types[this.chartType.value];
        this.extraControls = {};
        if (t.extra_controls) {
            const controls = t.extra_controls(this.dataSource.dataStore);
            const parentDiv = this.paramDiv;
            for (const c of controls) {
                createEl(
                    "div",
                    {
                        text: c.label,
                        classes: ["ciview-title-div"],
                    },
                    parentDiv
                );
                // could this logic be shared with SettingsDialog? <<
                if (c.type === "dropdown") {
                    const sel = createEl(
                        "select",
                        {
                            styles: {
                                maxWidth: "200px",
                            },
                        },
                        parentDiv
                    );
                    createFilterElement(sel, parentDiv);
                    for (const item of c.values) {
                        const option = createEl(
                            "option",
                            { text: item.name, value: item.value },
                            sel
                        );
                        if (item.value === c.defaultVal) {
                            option.selected = true;
                        }
                    }
                    this.extraControls[c.name] = sel;
                } else if (c.type === "string") {
                    const el = createEl(
                        "input",
                        { value: c.defaultVal },
                        parentDiv
                    );
                    this.extraControls[c.name] = el;
                    //el.onchange // not using callback, value will be read on submit().
                } else if (c.type === "textbox") {
                    const el = createEl(
                        "textarea",
                        { value: c.defaultVal, styles: { height: "300px" } },
                        parentDiv
                    );
                    this.extraControls[c.name] = el;
                    //el.onchange // not using callback, value will be read on submit().
                } else if (c.type === "checkbox" || c.type === "check") {
                    const el = createEl(
                        "input",
                        { type: "checkbox" },
                        parentDiv
                    );
                    el.checked = c.defaultVal;
                    this.extraControls[c.name] = el;
                    //el.onchange // not using callback, value will be read on submit().
                }
            }
        }
    }
}
