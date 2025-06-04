import SVGChart from "./SVGChart.js";
import { action } from "mobx";

class CategoryChart extends SVGChart {
    constructor(dataStore, div, config, axisTypes) {
        if (!Array.isArray(config.param)) {
            // as much as it's nice to have a simple single value, we are now more formally consistently
            // making config.param an array, so we need to make sure it is one here.
            // Currently, the amount of extra faff using config.param[0] here 
            // is much less than would be needed elsewhere with more generic code.
            config.param = [config.param];
        }
        // would be nice to have a more generic way to handle this...
        // if we have an explicitly set title, we keep it, otherwise we use the column name & update it if it changes
        const colName = dataStore.getColumnName(config.param[0]);
        const hadDefaultTitle = !config.title || config.title === colName;
        config.title = config.title || colName;
        const originalTitle = config.title;
        super(dataStore, div, config, axisTypes);
        this.dim = this.dataStore.getDimension("category_dimension");
        this.colors = this.dataStore.getColumnColors(config.param[0]);
        this.filter = [];
        this.mobxAutorun(() => {
            // we need to react to config.param changes and updateData() or similar...
            console.log("config.param changed", this.config.param);
            setTimeout(action(() => {
                // colors could probably be `@computed` from dataStore.getColumnColors(this.config.param[0])...
                this.colors = this.dataStore.getColumnColors(this.config.param[0]);
                this.onDataFiltered();
                // this.updateData();
                if (hadDefaultTitle && this.config.title === originalTitle) {
                    this.config.title = dataStore.getColumnName(this.config.param[0]);
                }
            }), 0);
        });
    }

    remove(notify = true) {
        this.dim.destroy(notify);
        super.remove();
    }

    setSize(x, y) {
        super.setSize(x, y);
        this.drawChart();
    }

    removeFilter() {
        this.dim.removeFilter();
        this.filter = [];
        this.drawChart();
    }

    filterCategories(cat, append) {
        if (append) {
            if (this.filter.indexOf(cat) !== -1) {
                this.filter = this.filter.filter((x) => x !== cat);
            } else if (this.filter.length === 0) {
                const vs = this.dataStore.getColumnValues(this.config.param[0]);
                this.filter = vs.filter((x) => x !== cat);
            } else {
                this.filter.push(cat);
            }
        } else {
            this.filter = [cat];
        }
        this.resetButton.style.display = "inline";
        this.drawChart(100);
        this.dim.filter("filterCategories", [this.config.param[0]], this.filter);
    }

    getFilter() {
        const f = {};
        if (!this.filter || this.filter.length === 0) {
            return null;
        }
        f[this.config.param[0]] = this.filter.slice(0);
        return f;
    }

    pinChart() {
        this.isPinned = true;
    }

    unpinChart() {
        this.isPinned = false;
        this.onDataFiltered();
    }

    updateData() {
        const c = this.config;

        this.rowData = this.data.slice(0);
        if (c.filter_zeros) {
            this.rowData = this.rowData.filter((x) => x[0] !== 0);
        }
        if (c.sort === "size") {
            this.rowData.sort((a, b) => b[0] - a[0]);
        }
        if (c.sort === "name") {
            const v = this.dataStore.getColumnValues(c.param);
            this.rowData.sort((a, b) => v[a[1]].localeCompare(v[b[1]]));
        }
        if (this.rowData.length > 40 && !c.show_limit) {
            c.show_limit = 40;
        }
        if (c.show_limit) {
            if (c.show_limit < this.rowData.length)
                this.rowData.splice(c.show_limit);
        }
    }

    onDataFiltered(dim) {
        //no need to change anything
        if (this.dim === dim || this.isPinned) {
            return;
        }
        if (dim === "all_removed") {
            this.filter = [];
            this.resetButton.style.display = "none";
        }
        const config = {};
        this.dim.getCategories(
            (data) => {
                this.data = new Array(data.length);
                this.maxCount = 1;
                const c = this.config;
                for (let n = 0; n < data.length; n++) {
                    this.data[n] = [data[n], n];
                    this.maxCount = Math.max(data[n], this.maxCount);
                }

                this.updateData();
                this.drawChart();
            },
            this.config.param[0],
            config,
        );
    }
}

export default CategoryChart;
