import "gridstack/dist/gridstack.min.css";
import { GridStack } from "gridstack";
import { debounce } from "../utilities/Utilities";
import type { ChartManager, DataSource } from "./charts";
import type BaseChart from "./BaseChart";
export type Chart = BaseChart<unknown>;
function clearPosition(div: HTMLElement) {
    div.style.position = "";
    div.style.left = "";
    div.style.right = "";
    div.style.top = "";
    div.style.bottom = "";
    div.style.width = "";
    div.style.height = "";
}

declare global {
    interface HTMLElement {
        gridstackNode?: any;
        gridstackPopoutCallback?: () => void;
    }
}
export type GridInstance = {
    grid: GridStack;
    charts: Set<Chart>;
    icon: HTMLElement;
};

export default class GridStackManager {
    cellHeight = 150;
    cellApproxWidth = 300;
    grids: Map<DataSource, GridInstance>;
    dragHandle = ".ciview-chart-title";
    chartManager: ChartManager;
    constructor(chartManager: ChartManager) {
        this.grids = new Map();
        this.chartManager = chartManager;
    }

    getGrid(ds: DataSource) {
        if (!this.grids.has(ds)) {
            const div = ds.contentDiv;
            div.classList.add("grid-stack");
            const rect = div.getBoundingClientRect();
            const column = Math.round(rect.width / this.cellApproxWidth);
            console.log(column);
            const rows = Math.round(rect.height / this.cellHeight);
            this.cellHeight = rect.height / rows; //hack, we might want different sizes for different views, not one 'global'
            console.log("gridstack cellHeight", this.cellHeight);
            const grid = GridStack.init(
                {
                    cellHeight: this.cellHeight,
                    handle: this.dragHandle,
                    float: true,
                    oneColumnSize: 400,
                    resizable: {
                        handles: "e,se,s,sw,w",
                    }
                    // these options not working as expected...
                    // margin: 10
                    // column
                },
                div,
            );
            grid.on("resizestop", (_: Event, el: HTMLElement) => {
                el.style.filter = "";
            });
            grid.on("resizestart", (_: Event, el: HTMLElement) => {
                el.style.filter = "blur(1px) opacity(0.5)";
            });

            const i = this.chartManager.addMenuIcon(
                ds.name,
                "fas fa-compress-arrows-alt",
                "compact layout",
                () => grid.compact(),
            );
            this.grids.set(ds, { grid: grid, charts: new Set(), icon: i });
        }

        return this.grids.get(ds);
    }

    destroy(ds: DataSource) {
        if (!this.grids.has(ds)) {
            console.error(`destroying non-existent gridstack - ${ds.name}`);
            return;
        }
        const gi = this.grids.get(ds);
        if (!gi) {
            console.error(`destroying non-existent gridstack - ${ds.name}`);
            return;
        }
        const div = ds.contentDiv;
        //store sizes/positions of div elements
        const sizes = new Map();
        for (const chart of gi.charts) {
            const d = chart.getDiv();
            sizes.set(d, [
                d.offsetWidth,
                d.offsetHeight,
                d.offsetLeft,
                d.offsetTop,
            ]);
            chart.removeLayout?.();
            chart.config.gssize = undefined;
            chart.config.gsposition = undefined;
        }

        gi.grid.destroy(false);

        //convert back to absolute positioning plus other clean up on the div
        for (const chart of gi.charts) {
            const d = chart.getDiv();
            const s = sizes.get(d);
            d.style.position = "absolute";
            d.style.width = `${s[0] - 5}px`;
            d.style.height = `${s[1] - 5}px`;
            d.style.left = `${s[2]}px`;
            d.style.top = `${s[3]}px`;
            chart.config.size = [s[0] - 5, s[1] - 5];
            chart.config.position = [s[2], s[3]];
            d.classList.remove("grid-stack-item");
        }
        div.classList.remove("grid-stack");
        gi.icon.remove();
        this.grids.delete(ds);
    }

    getCellDimensions(ds: DataSource) {
        const gridInstance = this.getGrid(ds);
        if (!gridInstance) {
            console.error("no grid instance for", ds.name);
            return [300, 150];
        }
        const grid = gridInstance.grid;
        return [grid.cellWidth(), this.cellHeight];
    }

    manageChart(chart: Chart, ds: DataSource, autoPosition?: boolean) {
        const gridInstance = this.getGrid(ds);
        if (!gridInstance) throw new Error(`no grid instance for ${ds.name}`);
        const grid = gridInstance.grid;
        gridInstance.charts.add(chart);
        const div = chart.getDiv();
        const rect = div.getBoundingClientRect();
        let w = Math.round(rect.width / grid.cellWidth());
        let h = Math.round(rect.height / this.cellHeight);
        let x = Math.round(div.offsetLeft / grid.cellWidth()); //relative location
        //cases where the old container was wider than the current one
        x = x > 11 ? 11 : x;
        let y = Math.floor(div.offsetTop / this.cellHeight); //relative location
        if (chart.config.gsposition) {
            [x, y] = chart.config.gsposition;
            [w, h] = chart.config.gssize;
            autoPosition = false;
        }
        //maybe better not to add to dom before getting here, and use grid.addWidget
        //current version seems to be working moderately ok.
        clearPosition(div);
        div.style.transition = "filter 0.5s ease";

        // using ResizeObserver rather than gridstack callbacks
        // means we have closure on 'chart' without needing to maintain another data structure
        // or attach more properties to the div.
        const ro = new ResizeObserver(
            debounce(() => {
                try {
                    chart.setSize();
                } catch {}
            }, 20),
        );
        ro.observe(div);

        // div.remove();
        // grid.addWidget(div, {w, h, x, y, autoPosition});

        grid.makeWidget(div);
        //nb, autoPosition property doesn't apply in update()?
        //passing options to makeWidget() or addWidget() does not evoke joy.
        if (!autoPosition) {
            console.log("gridstack:", chart.config.type, { w, h, x, y });
            grid.update(div, { w, h, x, y });
        } else {
            grid.update(div, { w, h });
        }

        const oldChangeBase = addPopOutHandler(); //revertModifications() used here refers to values stored by subsequent calls
        const lockButton = addPositionLock();
        const oldRemove = addRemoveHandler();
        function addPositionLock() {
            let locked = false;
            const lockButton = chart.addMenuIcon(
                "fas fa-unlock",
                "lock position",
            );
            const lockIcon = lockButton.children[0];
            lockButton.addEventListener("click", () => {
                locked = !locked;
                lockIcon.classList.toggle("fa-lock", locked);
                lockIcon.classList.toggle("fa-unlock", !locked);
                grid.update(div, { locked });
                div.classList.toggle("gridLock", locked);
            });
            return lockButton;
        }
        function addRemoveHandler() {
            const { remove } = chart;
            chart.remove = () => {
                grid.removeWidget(div, true);
                gridInstance?.charts.delete(chart);

                //div.gridstackNode = undefined; //doesn't help
                remove.apply(chart);
            };
            return remove;
        }
        function addPopOutHandler() {
            const { changeBaseDocument } = chart;
            div.gridstackPopoutCallback = () => {
                // console.log('removing from gridstack');
                grid.removeWidget(div, true);
            };
            chart.changeBaseDocument = (doc: Document) => {
                //by now, it's too late... the parent element is no longer the grid.
                //grid.removeWidget(div, true);
                changeBaseDocument.apply(chart, [doc]);
                if (doc === document)
                    console.error(
                        "changeBaseDocument should be restored to normal by the time we get here...",
                    );
                revertModifications();
            };
            return changeBaseDocument;
        }
        // this happens when we close PopOut... we expect to add them all again straight after
        // >>> this implementation might change if we change more of the ChartManager logic for when manageChart is called <<<
        function revertModifications() {
            grid.removeWidget(div, true); //doesn't remove listeners from handle...
            //// leaving listeners as a bug for now.
            // lockButton.remove();
            chart.remove = oldRemove;
            chart.changeBaseDocument = oldChangeBase;
        }
        chart.removeLayout = () => {
            grid.removeWidget(div, false); //doesn't remove listeners from handle...
            lockButton.remove();
            div.gridstackPopoutCallback = undefined;
            chart.remove = oldRemove;
            chart.changeBaseDocument = oldChangeBase;
            ro.disconnect();
        };
    }
}

/// some more layout stuff, subject to change / moving to a different file

export type P = [number, number];
export type Config = Partial<{ size: P; position: P; gssize: P; gsposition: P }>;
function getVisibleChartBounds(dataSource: DataSource) {
    const charts = Object.entries(
        window.mdv.chartManager.charts,
    ) as unknown as [[string, { dataSource: DataSource; chart: Chart }]];
    return charts
        .filter((c) => c[1].dataSource === dataSource)
        .map((c) => c[1].chart.getDiv().getBoundingClientRect());
}

function getGridInfo(dataSource: DataSource) {
    const [cellW, cellH] =
        window.mdv.chartManager.gridStack.getCellDimensions(dataSource);
    const rect = dataSource.contentDiv.getBoundingClientRect();
    const rows = Math.floor(rect.height / cellH) * 2; //hack for more rows in layout
    const cols = Math.floor(rect.width / cellW);
    return { cellW, cellH: cellH / 2, rows, cols };
}

function findFreeSpace(dataSource: DataSource) {
    const occupiedRects = getVisibleChartBounds(dataSource);
    const { cellW, cellH, rows, cols } = getGridInfo(dataSource);
    const mainRect = dataSource.contentDiv.getBoundingClientRect();
    const width = 2;
    const height = 1;

    // Initialize the grid
    const grid = new Array(rows);
    for (let i = 0; i < rows; i++) {
        grid[i] = new Array(cols).fill(false);
    }

    // Mark the occupied cells
    for (const rect of occupiedRects) {
        // get outer-bounds in terms of nominal grid... with appropriate offset relative to mainRect...
        const left = Math.floor((rect.left - mainRect.left) / cellW);
        const top = Math.floor((rect.top - mainRect.top) / cellH);
        const right = Math.min(
            cols,
            Math.ceil((rect.right - mainRect.left) / cellW),
        );
        const bottom = Math.min(
            rows,
            Math.ceil((rect.bottom - mainRect.top) / cellH),
        );
        for (let row = top; row < bottom; row++) {
            for (let col = left; col < right; col++) {
                grid[row][col] = true;
            }
        }
    }

    // Find a free space
    for (let row = 0; row <= rows - height; row++) {
        for (let col = 0; col <= cols - width; col++) {
            let isFree = true;
            for (let i = 0; i < height && isFree; i++) {
                for (let j = 0; j < width && isFree; j++) {
                    if (grid[row + i][col + j]) {
                        isFree = false;
                    }
                }
            }
            if (isFree) {
                return { left: col * cellW, top: row * cellH };
            }
        }
    }

    // If no free space was found, return [10, 10]
    return { left: 10, top: 10 };
}

export function positionChart(dataSource: DataSource, config: Config) {
    // some legacy code extracted to here and may be reviewed at some point.
    let width = 300;
    let height = 300; //consider having a preferred-size...
    let left = 10;
    let top = 10;
    if (config.size) {
        width = config.size[0];
        height = config.size[1];
    }
    const { chartManager } = window.mdv;
    if (config.position) {
        left = config.position[0];
        top = config.position[1];
    } else {
        try {
            const pos = findFreeSpace(dataSource);
            left = pos.left;
            top = pos.top;
        } catch (e) {
            console.error(e);
        }
    }
    //hack approx position of grid stack elements

    if (
        chartManager.viewData.dataSources[dataSource.name].layout ===
            "gridstack" &&
        config.gssize
    ) {
        const cellDim = chartManager.gridStack.getCellDimensions(dataSource);
        width = Math.round(config.gssize[0] * cellDim[0]);
        height = Math.round(config.gssize[1] * cellDim[1]);
        left = config.gsposition ? Math.round(config.gsposition[0] * (cellDim[0] + 5)) : 20;
        top = config.gsposition ? Math.floor(config.gsposition[1] * (cellDim[1] + 5)) : 20;
    }

    return { width, height, left, top };
}
