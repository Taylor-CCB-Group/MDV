import BaseChart from "./BaseChart";
import { createEl, makeSortable } from "../utilities/Elements.js";
import noUiSlider from "nouislider";
import { getRandomString } from "../utilities/Utilities";
import MultiSelectDropdown from "../utilities/MultiSelect";

class SelectionDialog extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        const c = this.config;
        const con = c.param;
        this.dims = {};
        this.textFilters = {};
        this.numberFilters = {};
        c.filters = c.filters || {};
        this.hasFiltred = false;
        this.contentDiv.style.overflowY = "auto";

        for (const c of con) {
            const col = dataStore.columnIndex[c];
            const div = createEl(
                "div",
                {
                    styles: { padding: "5px" },
                },
                this.contentDiv,
            );
            //marker to identify which div associated with which column
            div.__fcol__ = col.field;
            createEl(
                "div",
                {
                    text: col.name,
                    classes: ["i-title"],
                    styles: { fontWeight: "bold" },
                },
                div,
            );
            if (col.datatype.match(/text/)) {
                this.addTextFilter(col, div);
            } else if (
                col.datatype === "integer" ||
                col.datatype === "double" ||
                col.datatype === "int32"
            ) {
                this.addNumberFilter(col, div);
            }
        }
        makeSortable(this.contentDiv, {
            handle: "i-title",
            //re-arrange params in config
            sortEnded: (li) => {
                this.config.param = li.map((x) => x.__fcol__);
            },
        });
        if (this.hasFiltered) {
            setTimeout(() => this.dataStore.triggerFilter(), 0);
        }
    }

    removeFilter() {
        const f = this.config.filters;
        for (const c in this.textFilters) {
            const t = this.textFilters[c];
            t.unSelectAll();
            f[c].category = [];
        }
        for (const c in this.numberFilters) {
            const mm = this.dataStore.getMinMaxForColumn(c);
            const t = this.numberFilters[c];
            t.noUiSlider.set([mm[0], mm[1]]);
            f[c] = null;
        }
        this.resetButton.style.display = "none";
        for (const c in this.dims) {
            this.dims[c].removeFilter(false);
        }
        this.dataStore.triggerFilter();
    }

    onDataFiltered(dim) {}

    addNumberFilter(col, div) {
        const sl = createEl(
            "div",
            {
                styles: {
                    margin: "10px 10px",
                },
            },
            div,
        );
        const dd = createEl(
            "div",
            {
                styles: {
                    padding: "3px 10px",
                    overflow: "hidden",
                },
            },
            div,
        );
        createEl(
            "span",
            { text: ">", styles: { float: "left", fontWeight: "bold" } },
            dd,
        );
        const greaterThan = createEl(
            "input",
            {
                styles: {
                    width: "70px",
                    float: "left",
                },
            },
            dd,
        );
        createEl(
            "span",
            { text: "<", styles: { float: "right", fontWeight: "bold" } },
            dd,
        );
        const lessThan = createEl(
            "input",
            {
                styles: {
                    width: "70px",
                    float: "right",
                },
            },
            dd,
        );

        const c = this.config;
        const dim = this.dataStore.getDimension("range_dimension");
        dim.noClear = true;
        this.dims[col.field] = dim;
        const mm = this.dataStore.getMinMaxForColumn(col.field);

        const fil = c.filters[col.field];
        let cv = [mm[0], mm[1]];
        if (fil) {
            cv = [fil[0], fil[1]];
        }
        noUiSlider.create(sl, {
            start: [cv[0], cv[1]],
            range: {
                min: mm[0],
                max: mm[1],
            },
            //step:s.step || null,
            tooltips: true,
            documentElement: this.__doc__,
        });
        sl.noUiSlider.on("set", (values) => {
            const min = Number.parseFloat(values[0]);
            const max = Number.parseFloat(values[1]);
            if (min <= mm[0] && max >= mm[1]) {
                c.filters[col.field] = null;
            } else {
                c.filters[col.field] = [min, max];
            }
            lessThan.value = max;
            greaterThan.value = min;

            this.filterCategories(col.field);
        });

        lessThan.value = cv[1];
        greaterThan.value = cv[0];

        greaterThan.addEventListener("blur", (e) => {
            let n = Number.parseFloat(greaterThan.value);

            const cvi = sl.noUiSlider.get();
            if (Number.isNaN(n) || n > cvi[1]) {
                return;
            }
            n = n < mm[0] ? mm[0] : n;
            sl.noUiSlider.set([n, cvi[1]], true, true);
        });
        greaterThan.addEventListener("keypress", (e) => {
            if (e.key === "Enter") {
                greaterThan.blur();
            }
        });

        lessThan.addEventListener("blur", (e) => {
            let n = Number.parseFloat(lessThan.value);

            const cvi = sl.noUiSlider.get();
            if (Number.isNaN(n) || n < cvi[0]) {
                return;
            }
            n = n > mm[1] ? mm[1] : n;
            sl.noUiSlider.set([cvi[0], n], true, true);
        });
        lessThan.addEventListener("keypress", (e) => {
            if (e.key === "Enter") {
                lessThan.blur();
            }
        });

        if (fil) {
            this.filterCategories(col.field, false);
            this.hasFiltered = true;
        }
        this.numberFilters[col.field] = sl;
    }

    filterCategories(col, notify = true) {
        const f = this.config.filters[col];
        const type = this.dataStore.columnIndex[col].datatype;
        if (type === "text" || type === "text16") {
            if (!f.category || f.category.length === 0) {
                this.dims[col].removeFilter();
            } else {
                this.resetButton.style.display = "inline";
                this.dims[col].filter(
                    "filterCategories",
                    [col],
                    f.category,
                    notify,
                );
            }
        } else if (type === "multitext") {
            if (f.category.length === 0) {
                this.dims[col].removeFilter();
            } else {
                this.resetButton.style.display = "inline";
                const args = f.category.slice(0);
                args.operand = f.operand;
                this.dims[col].filter("filterCategories", [col], args, notify);
            }
        } else {
            if (!f) {
                this.dims[col].removeFilter();
            } else {
                this.resetButton.style.display = "inline";
                this.dims[col].filter(
                    "filterRange",
                    [col],
                    { min: f[0], max: f[1] },
                    notify,
                );
            }
        }
    }

    addTextFilter(col, div) {
        const ismulti = col.datatype === "multitext";
        //div.style.whiteSpace="nowrap";
        const c = this.config;
        const dim = this.dataStore.getDimension("category_dimension");
        dim.noClear = true;
        this.dims[col.field] = dim;
        let fil = c.filters[col.field];
        if (!fil) {
            fil = c.filters[col.field] = {
                category: [],
                operand: ismulti ? "or" : null,
            };
        } else {
            //convert legacy
            if (typeof fil.category === "string") {
                if (fil.category === "__none__") {
                    fil.category = [];
                } else {
                    fil.category = [fil.category];
                }
            }
            if (fil.exclude) {
                fil.category = this.dataStore
                    .getColumnValues(col.field)
                    .filter((x) => x !== fil.category[0]);
            }
            //end
            this.filterCategories(col.field, false);
            this.hasFiltered = true;
        }

        const hdiv = createEl("div", {}, div);
        if (ismulti) {
            const div1 = createEl("div", {}, hdiv);
            const rname = getRandomString();
            createEl("span", { text: "and" }, div1);
            const rb = createEl(
                "input",
                {
                    type: "radio",
                    name: rname,
                    value: "and",
                },
                div1,
            );
            rb.checked = fil.operand === "and";
            rb.addEventListener("click", (e) => {
                fil.operand = "and";
                this.filterCategories(col.field);
            });
            createEl("span", { text: "or" }, div1);
            const ra = createEl(
                "input",
                {
                    type: "radio",
                    name: rname,
                    value: "or",
                },
                div1,
            );
            ra.checked = fil.operand === "or";
            ra.addEventListener("click", (e) => {
                fil.operand = "or";
                this.filterCategories(col.field);
            });
        }
        const dd = createEl("div", {}, div);
        const options = col.values
            .slice(0)
            .sort((a, b) => a.localeCompare(b))
            .map((x) => {
                return {
                    text: x,
                    value: x,
                    selected: fil.category.indexOf(x) !== -1,
                };
            });
        const ms = new MultiSelectDropdown(dd, options, {
            onchange: (vals) => {
                fil.category = vals;
                this.filterCategories(col.field);
            },
        });
        this.textFilters[col.field] = ms;
    }

    remove(notify = true) {
        for (const c in this.dims) {
            this.dims[c].destroy(false);
        }
        Object.values(this.textFilters).forEach((x) => x.remove());
        if (notify) {
            this.dataStore.triggerFilter();
        }
        super.remove();
    }
}

BaseChart.types["selection_dialog"] = {
    name: "Selection Dialog",
    class: SelectionDialog,
    params: [
        {
            type: "_multi_column:all",
            name: "Columns To filter",
        },
    ],
};

export default SelectionDialog;
