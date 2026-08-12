import { BaseReactChart } from "../../react/components/BaseReactChart";
import { g } from "../../lib/utils";
import BaseChart from "../BaseChart";
import CategoryHeatmapComponent from "./CategoryHeatmapComponent";

class CategoryHeatMap extends BaseReactChart {
    constructor(dataStore, div, config) {
        if (!Array.isArray(config.param)) {
            config.param = [config.param];
        }
        const xName = config.param[0] ? dataStore.getColumnName(config.param[0]) : "Category";
        const yName = config.param[1] ? dataStore.getColumnName(config.param[1]) : "Category";
        config.title = config.title || `${xName} vs ${yName}`;
        config.x_display_categories = config.x_display_categories || [];
        config.y_display_categories = config.y_display_categories || [];
        super(dataStore, div, config, CategoryHeatmapComponent);
        this.dim = this.dataStore.getDimension("category_dimension");
        this.cellFilter = null;
    }

    remove(notify = true) {
        this.dim.destroy(notify);
        super.remove();
    }

    fetchAggregation() {
        const params = this.config.param || [];
        if (params.length < 2) {
            return Promise.resolve({
                xLabels: [],
                yLabels: [],
                counts: [],
                maxCount: 0,
                totalCells: 0,
            });
        }

        return new Promise((resolve) => {
            this.dim.getCategoryHeatmap(
                (aggregation) => {
                    resolve(aggregation);
                },
                [params[0], params[1]],
                {},
            );
        });
    }

    filterCell(xCategory, yCategory) {
        const [xColumn, yColumn] = this.config.param || [];
        if (!xColumn || !yColumn) {
            return;
        }

        if (
            this.cellFilter &&
            this.cellFilter[0] === xCategory &&
            this.cellFilter[1] === yCategory
        ) {
            this.removeFilter();
            return;
        }

        this.cellFilter = [xCategory, yCategory];
        this.resetButton.style.display = "inline";
        this.dim.filter(
            "filterCategoryPair",
            [xColumn, yColumn],
            [xCategory, yCategory],
        );
    }

    removeFilter() {
        this.cellFilter = null;
        this.dim.removeFilter();
        this.resetButton.style.display = "none";
    }

    getFilter() {
        if (!this.cellFilter) {
            return null;
        }
        const [xColumn, yColumn] = this.config.param || [];
        if (!xColumn || !yColumn) {
            return null;
        }
        return {
            [xColumn]: [this.cellFilter[0]],
            [yColumn]: [this.cellFilter[1]],
        };
    }

    getSettings() {
        const settings = super.getSettings();
        const [xColumn, yColumn] = this.config.param || [];
        const xValues = xColumn ? this.dataStore.getColumnValues(xColumn) : [];
        const yValues = yColumn ? this.dataStore.getColumnValues(yColumn) : [];

        return settings.concat([
            g({
                type: "multidropdown",
                current_value: this.config.x_display_categories || [],
                values: [xValues.map((x) => [x]), 0, 0],
                label: "Display X Categories",
                func: (selected) => {
                    this.config.x_display_categories = selected || [];
                },
            }),
            g({
                type: "multidropdown",
                current_value: this.config.y_display_categories || [],
                values: [yValues.map((y) => [y]), 0, 0],
                label: "Display Y Categories",
                func: (selected) => {
                    this.config.y_display_categories = selected || [];
                },
            }),
        ]);
    }
}

BaseChart.types["category_heatmap"] = {
    name: "Category Heatmap",
    class: CategoryHeatMap,
    params: [
        {
            type: ["text", "text16", "multitext"],
            name: "X-axis categories",
        },
        {
            type: ["text", "text16", "multitext"],
            name: "Y-axis categories",
        },
    ],
};

export default CategoryHeatMap;
