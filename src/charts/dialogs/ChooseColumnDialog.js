import { createEl } from "../../utilities/Elements";
import { BaseDialog } from "../../utilities/Dialog.js";
import { getRandomString } from "../../utilities/Utilities";

/**
 * Creates a dialog for the user to choose multiple columns
 * @param {DataStore} dataStore - the dataStore the columns will be chosen from
 * @param {function} callback - A function called when the user has selected the columns
 * The callback is provided with a list of chosen column fields(ids)
 * @param {string} [filter=all] - The type of column the use can choose
 */
export class ChooseColumnDialog extends BaseDialog {
    constructor(dataStore, callback, filter = "all") {
        const config = {
            footer: true,
            width: 250,
            maxHeight: 500,
            title: "Select Columns",
            buttons: [{ text: "OK", method: "getColumns" }],
        };
        super(config, {
            dataStore: dataStore,
            callback: callback,
            filter: filter,
        });
    }
    init(content) {
        this.ds = content.dataStore;
        const gd = createEl("div", { styles: { padding: "8px" } });
        const rName = getRandomString();
        createEl("div", { text: "Groups" }, this.dialog);

        const cgs = Object.keys(this.ds.columnGroups);
        cgs.unshift("All");
        for (const group of cgs) {
            const d = createEl(
                "span",
                {
                    styles: {
                        display: "inline-block",
                        whiteSpace: "nowrap",
                        marginRight: "5px",
                    },
                },
                gd
            );
            createEl("span", { text: group }, d);
            createEl(
                "input",
                {
                    type: "radio",
                    value: group,
                    name: rName,
                },
                d
            ).addEventListener("click", (e) => {
                this.checkAllInGroup(e.target.value);
            });
        }
        this.dialog.append(gd);
        createEl("div", { text: "Select Individual Columns" }, this.dialog);
        const cd = createEl("div", { style: { padding: "8px" } });
        const cols = this.ds.getColumnList(content.filter);
        this.checks = [];
        this.callback = content.callback;
        //todo let this work with createFilterElement? <-- desperately needed in ad-car... hacking something together for now
        const filter = createEl(
            "input",
            {
                placeholder: "Search columns",
                type: "text",
            },
            cd
        );
        filter.addEventListener("input", (e) => {
            const val = e.target.value.toLowerCase().split(" ");
            for (const check of this.checks) {
                const f = val.some((v) => !check[1].toLowerCase().includes(v));
                check[0].parentElement.style.display = f ? "none" : "";
            }
        });
        for (const col of cols) {
            const d = createEl(
                "div",
                {
                    styles: {
                        //display:"inline-block",
                        whiteSpace: "nowrap",
                        // marginRight:"5px"
                    },
                },
                cd
            );
            //createEl("span",{text:col.name},d);
            const cb = createEl(
                "input",
                {
                    type: "checkbox",
                },
                d
            );
            this.checks.push([cb, col.field]);
            createEl("span", { text: col.name }, d);
        }
        this.dialog.append(cd);
    }
    checkAllInGroup(group) {
        if (group === "All") {
            for (const check of this.checks) {
                check[0].checked = true;
            }
        } else {
            const cols = this.ds.columnGroups[group].columns;
            for (const check of this.checks) {
                check[0].checked = cols.indexOf(check[1]) !== -1;
            }
        }
    }

    getColumns() {
        const cols = [];
        for (const check of this.checks) {
            if (check[0].checked) {
                cols.push(check[1]);
            }
        }
        this.callback(cols);
        this.close();
    }
}
