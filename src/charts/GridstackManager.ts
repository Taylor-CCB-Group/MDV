import "gridstack/dist/gridstack.min.css";
import { GridStack } from "gridstack";
import { debounce } from "../utilities/Utilities";
import type { ChartManager, DataSource } from "./charts";
import type BaseChart from "./BaseChart";
export type Chart = BaseChart<any>;
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
type ManagedChart = {
    outer: HTMLElement;
    inner: HTMLElement;
    lockButton: HTMLElement;
    _gsResizeObserver?: ResizeObserver;
    _gsMutationObserver?: MutationObserver;
} & Chart;
export type GridInstance = {
    grid: GridStack;
    charts: Set<ManagedChart>;
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
                        handles: "e,se,s,sw,w,n,ne,nw",
                    },
                    // now working - as a have inner and outer containers
                     margin: "3px"
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
                () => this.compact(ds),
            );
            this.grids.set(ds, { grid: grid, charts: new Set(), icon: i });
        }

        return this.grids.get(ds);
    }

    //when the grid calls compact the size and move observers are not called
    //therefore need to do it manually
    compact(ds: DataSource) {
        const gi = this.grids.get(ds);
        if (!gi) {
            console.error(`no gridstack for ${ds.name}`);
            return;
        }
        gi.grid.compact();
        //store gs position of each chart
        for (const chart of gi.charts) {
            const div = chart.getDiv();
            const d = div.parentElement?.parentElement;
            if (d){
                const { x, y, w, h } = d.gridstackNode;
                chart.config.gssize = [w, h];
                chart.config.gsposition = [x, y];
            }

        }
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
        //get an array - don't want to loop through set as removeLayout 
        //will remove the chart from the set
        const charts = Array.from(gi.charts);
        for (const chart of charts){
            const cdiv = chart.getDiv();
            //bit of a hack to get the outer gridstack container
            const d = cdiv.parentElement?.parentElement;
            if (d){
                //get the size and position of the chart in gridstack
                const w =  d.offsetWidth;
                const h = d.offsetHeight;
                const l = d.offsetLeft;
                const t= d .offsetTop;
                //apply these to the chart's div for absolute positioning
                cdiv.style.position = "absolute";
                cdiv.style.width = `${w - 6}px`;
                cdiv.style.height = `${h - 12}px`;
                cdiv.style.left = `${l}px`;
                cdiv.style.top = `${t}px`;
                chart.config.size = [w- 6, h - 12];
                chart.config.position = [l,t]
            }
            //remove the gridstack data from the chart's config
            chart.config.gssize = undefined;
            chart.config.gsposition = undefined;
            chart.removeLayout?.();
        }

        gi.grid.destroy(false);
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

    manageChart(originalChart: Chart, ds: DataSource, autoPosition?: boolean, remanageChart?: boolean) {
        try {
            const gridInstance = this.getGrid(ds);
            if (!gridInstance) throw new Error(`no grid instance for ${ds.name}`);
            const grid = gridInstance.grid;
            //note that we're being a bit type-unsafe here
            //everything we do with managedChart will be to the same chart instance
            //we're adding properties to it and using it as a managed chart
            //would just take a little bit more re-work to make it properly safe
            const chart = originalChart as ManagedChart;
            gridInstance.charts.add(chart);
            const div = chart.getDiv();
            //work out the size location of the chart in gridstack
            //doesn't work that well with complex layouts
            const rect = div.getBoundingClientRect();
            let w = Math.round((rect.width+6) / grid.cellWidth());
            let h = Math.round((rect.height+6) / this.cellHeight);
            let x = Math.round(div.offsetLeft / grid.cellWidth()); //relative location
            //cases where the old container was wider than the current one
            x = x > 11 ? 11 : x;
            let y = Math.floor(div.offsetTop / this.cellHeight); //relative location
            //if the chart has gridstack data i.e. in saved view or being re-added from a popout
            //then use these
            if (chart.config.gsposition) {
                [x, y] = chart.config.gsposition;
                [w, h] = chart.config.gssize;
                autoPosition = false;
            }
            //maybe better not to add to dom before getting here, and use grid.addWidget
            //current version seems to be working moderately ok.
            clearPosition(div);
            div.style.transition = "filter 0.5s ease";


            // Add the "move" cursor
            const handle = div.querySelector(this.dragHandle) as HTMLElement;
            if (handle) {
                handle.style.cursor = "move"; // Add the "move" cursor
            }

            // a bit long-winded but for gridstack to work correctly we 
            // the chart wrapped in an outer and inner container
            //! check the handling these variables if any issues occur
            let outer: HTMLElement;
            let inner: HTMLElement;
            // Create new inner and outer divs only during the creation
            if (remanageChart) {
                outer = chart.outer;
                inner = chart.inner;
                if (!inner || !outer) {
                    // If the containers don't exist call manageChart again with the remanage as false
                    this.manageChart(chart, ds, autoPosition, false);
                    return;
                }
            } else {
                outer = document.createElement("div");
                inner = document.createElement("div");
                chart.outer = outer;
                chart.inner = inner;

                inner.classList.add("grid-stack-item-content");
                inner.style.overflow = "visible"; // able to see shadow
                outer.appendChild(inner);
                inner.appendChild(div);
                grid.el.appendChild(outer);
            }
            const origHeight = div.style.height;
            div.style.height="100%";

            //saves thr chart's position and size to the chart's config
            const saveChartPosition = () => { 
                const { x, y, w, h } = outer.gridstackNode;
                chart.config.gssize = [w, h];
                chart.config.gsposition = [x, y];
            }

            // Disconnect existing observers
            //! If any problems occur try deleting the observers
            if (remanageChart) {
                if (chart._gsResizeObserver) {
                chart._gsResizeObserver.disconnect();
                }
                if (chart._gsMutationObserver) {
                chart._gsMutationObserver.disconnect();
                }
            }
            // using ResizeObserver rather than gridstack callbacks
            // means we have closure on 'chart' without needing to maintain another data structure
            // or attach more properties to the div.
            const ro = new ResizeObserver(
                debounce(() => {
                    try {
                        //this is causing 'No data for row chart' error
                        //when loading as saved view
                        //it redraws chart before the data is processed
                        //^ perhaps we could use something related to chart.deferredInit if this is a problem
                        //would possibly need to improve the design of that mechanism, 
                        //we could possibly move `ro` creation to after that promise resolves.

                        // also leads to a condition in which initial chart.config.size is not stable
                        // immediately after adding it... and also if the window is resized, 
                        // or the view is opened on a different screen... we need to think about how
                        // we serialize the layout information.
                        chart.setSize();
                        saveChartPosition();

                    } catch {}
                }, 20),
            );
            ro.observe(outer);
            // Storing the observer
            chart._gsResizeObserver = ro;

            //stores the position and size of the chart is repositioned
            //could also check if width, height has changed as well and have one observer,
            //but probably best to use the proper resize observer above
            const mo = new MutationObserver(
                debounce(() => {
                    try {
                        saveChartPosition();
                    } catch {}
                }, 20),
            );
            mo.observe(outer, { attributes: true, attributeFilter: ["style"] });
            // Storing the observer
            chart._gsMutationObserver = mo;
            // workaround for absolute positioning still being saved in config
            // this interferes with detecting unsaved changes
            // but if we want to compare state with state from server, it won't work.
            // it would still also have a false positive if window is resized.
            // only applies during chart init, so not when _alreadyCheckedDeferReady is true.
            if (!chart._alreadyCheckedDeferReady) {
                chart.deferredInit(() => chart.setSize(), 100);
            }

            // Create widget only during creation
            if (!remanageChart)
                grid.makeWidget(outer);

            //nb, autoPosition property doesn't apply in update()?
            //passing options to makeWidget() or addWidget() does not evoke joy.
            if (!autoPosition) {
                console.log("gridstack:", chart.config.type, { w, h, x, y });
                grid.update(outer, { w, h, x, y });
            } else {
                grid.update(outer, { w, h });
            }
        

            let lockButton: HTMLElement;
            // Create lock button only the first time during creation
            if (remanageChart) {
                lockButton = chart.lockButton as HTMLElement;
            } else {
                lockButton = addPositionLock();
                chart.lockButton = lockButton;
            }
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
                    grid.update(outer, { locked });
                    outer.classList.toggle("gridLock", locked);
                });
                return lockButton;
            }
            //could this be handled by a decorator?
            const oldRemove = addRemoveHandler();
            function addRemoveHandler() {
                const { remove } = chart;
                chart.remove = () => {
                    //grid.removeWidget(div, true);
                    chart.removeLayout?.();
                    remove.apply(chart);
                };
                return remove;
            }

            div.gridstackPopoutCallback = () => {
                chart.removeLayout?.();
            };

            //the below was leading subtle errors
            //not sure why you need to alter changeBaseDocument?
            /*
            function addPopOutHandler() {
                const { changeBaseDocument } = chart;
                div.gridstackPopoutCallback = () => {
                    // console.log('removing from gridstack');
                    grid.removeWidget(div, true);
                    if (chart.removeLayout)
                        chart.removeLayout();
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
            }*/

            //cleanup - this is called when the chart is removed from the gridstack
            //a) we transfer the chart to a popout
            //b) we remove the chart completely
            //c) we change the layout back to absolute positioning
            chart.removeLayout = () => {
                //remove the chart from the gridstack
                grid.removeWidget(outer, false);

                // assign the height of the chart only if origHeight exists
                if (origHeight) {
                    div.style.height = origHeight;
                }
                //take the chart out of the containers
                grid.el.append(div);
                //destroy the containers
                outer?.remove();
                inner?.remove();
                lockButton?.remove();
                div.gridstackPopoutCallback = undefined;
                chart.removeLayout = undefined;
                chart.remove = oldRemove;
                if (handle) {
                    handle.style.cursor = ""; // Reset to default cursor
                }
                //remove chart
                gridInstance.charts.delete(chart);
                //remove the observers
                if (chart._gsResizeObserver) {
                    chart._gsResizeObserver.disconnect();
                    // not using delete because of biome performance warning
                    chart._gsResizeObserver = undefined;
                }
                if (chart._gsMutationObserver) {
                    chart._gsMutationObserver.disconnect();
                    // not using delete because of biome performance warning
                    chart._gsMutationObserver = undefined;
                }
                ro.disconnect();
                mo.disconnect();
            };
        } catch (error) {
            console.error("Error occurred during manage chart:", error);
        }
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
