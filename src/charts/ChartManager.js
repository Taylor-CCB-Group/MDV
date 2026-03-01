import {
    createEl,
    makeDraggable,
    makeResizable,
    MDVProgress,
    removeDraggable,
    removeResizable,
    createMenuIcon,
    splitPane,
} from "../utilities/Elements";
import BaseChart from "./BaseChart";
import DataStore from "../datastore/DataStore.js";
import CustomDialog from "./dialogs/CustomDialog.js";
import { ContextMenu } from "../utilities/ContextMenu";
import { BaseDialog } from "../utilities/Dialog.js";
import { getRandomString } from "../utilities/Utilities";
import { csv, tsv, json } from "d3-fetch";
import ColorChooser from "./dialogs/ColorChooser";
import GridStackManager, { positionChart } from "./GridstackManager"; //nb, '.ts' unadvised in import paths... should be '.js' but not configured webpack well enough.
// this is added as a side-effect of import HmrHack elsewhere in the code, then we get the actual class from BaseDialog.experiment
import FileUploadDialogReact from "./dialogs/FileUploadDialogWrapper";

//default charts
import "./HistogramChart.js";
import "./RowChart.js";
import "./TableChart.js";
import "./WGL3DScatterPlot.js";
import "./WGLScatterPlot.js";
import "./RingChart.js";
import "./TextBoxChart.js";
import "./HeatMap.js";
import "./ViolinPlot.js";
import "./BoxPlot.js";
import "./SankeyChart.js";
import "./MultiLineChart.js";
import "./DensityScatterPlot";
// import "./SelectionDialog.js"; //now replaced with SelectionDialogReact (currently imported as a side-effect of import HmrHack)
import "./StackedRowChart";
import "./TreeDiagram";
import "./CellNetworkChart";
import "./FlexibleNetworkChart";
import "./SingleHeatMap";
import "./VivScatterPlot";
import "./DotPlot";
import "./ImageTableChart";
import "./CellRadialChart";
import "./RowSummaryBox";
import "./VivScatterPlot";
import "./ImageTableChart";
import "./ImageScatterChart";
// import "./WordCloudChart"; //todo: this only works in vite build, not webpack
import "./CustomBoxPlot";
import "./SingleSeriesChart";
import "./GenomeBrowser";
import "./DeepToolsHeatMap";
import connectIPC from "../utilities/InterProcessCommunication";
import { addChartLink } from "../links/link_utils";
import popoutChart from "@/utilities/Popout";
import { makeObservable, observable, action } from "mobx";
import { createMdvPortal } from "@/react/react_utils";
import ViewManager from "./ViewManager";
import ErrorComponentReactWrapper from "@/react/components/ErrorComponentReactWrapper";
import ViewDialogWrapper from "./dialogs/ViewDialogWrapper";
import { deserialiseParam, getConcreteFieldNames } from "./chartConfigUtils";
import AddChartDialogReact from "./dialogs/AddChartDialogReact";
import MenuBarWrapper from "@/react/components/MenuBarComponent";


//order of column data in an array buffer
//doubles and integers (both represented by float32) and int32 need to be first
// followed by multitext/text16 (uint16) then text/unique (uint8)
const column_orders = {
    double: 0,
    integer: 0,
    int32: 0,
    multitext: 1,
    text16: 1,
    text: 2,
    unique: 2,
};

const themes = {
    dark: {
        title_bar_color: "#222",
        main_panel_color: "black",
        text_color: "white",
        background_color: "#333",
    },
    light: {
        title_bar_color: "white",
        main_panel_color: "#f1f1f1",
        text_color: "black",
        background_color: "#bababa",
    },
};

//https://stackoverflow.com/questions/56393880/how-do-i-detect-dark-mode-using-javascript
function getPreferredColorScheme() {
    if (window.matchMedia) {
        if (window.matchMedia("(prefers-color-scheme: dark)").matches) {
            return "dark";
        }
        return "light";
    }
    return "light";
}
function listenPreferredColorScheme(callback) {
    if (window.matchMedia) {
        const colorSchemeQuery = window.matchMedia(
            "(prefers-color-scheme: dark)",
        );
        colorSchemeQuery.addEventListener("change", () =>
            callback(getPreferredColorScheme()),
        );
    }
}

/**
* The object to manage charts [Chart Manager](../../docs/extradocs/chartmanager.md)
* 
*/
export class ChartManager {
    /**
     * @param {string|HTMLElement} div - The DOM element or id of the element to house the app
     * @param {import("@/charts/charts").DataSource[]} dataSources - An array of datasource configs - see  [Data Source](../../docs/extradocs/datasource.md).
     * Each config must contain the size parameter, giving the number of rows in the DataStore.
     * @param {object} dataloader - An object containing the following
     *   - function - The [function](../../docs/extradocs/dataloader.md) to load the data
     *   (can be omitted if data loaded from a file)
     *   - viewLoader - The function that will load the each view  (not necessay if only one view)
     *   - rowDataLoader - (optional) an asunc function which is given the datasource name and row index
     *     returns unstructured data . A datasource's config requires row_data_loader:true to activate the loader
     *   - files - specifies the files to load the data 
     * @param {Object} config extra settings
     * @param {Object[]} [config.initialCharts] A list of chart configs to initially load if
     * no views are specified
     * @param {string[]} [config.all_views] A list of views that can be loaded (a suitable 
     * view loader is required to atually load the view)
     * @param {string} [config.current_view] the current view (only used with all views)
     * @param {string} [config.permission] the level of permission the user has. This just makes certain
     * options unavaliable. Any logic should be handled when a state_saved event is broadcast
     * @param {boolean} [config.gridstack] whether to arrange the charts in a grid
     * @param {boolean?} [config.chat_enabled] 
     * @param {string?} [config.mdv_api_root] 
     * @param {function} [listener] - A function to listen to events. `(eventType: string, cm: ChartManager, data: any) => void | Promise<void>`
     * beware: the way 'event listeners' are implemented is highly unorthodox and may be confusing.
     * 
     */
    constructor(div, dataSources, dataLoader, config = {}, listener = null) {
        if (!window.isSecureContext) {
            alert(
                "This application requires a secure context (https / localhost)",
            );
            throw new Error(
                "This application requires a secure context (https / localhost)",
            );
        }
        // manage global singleton
        if (!window.mdv) {
            // @ts-ignore this is unclean but somewhat ok
            window.mdv = {};
        }
        if (window.mdv.chartManager) {
            console.warn(
                "Making another ChartManager... had believed this to be a singleton",
            );
        }
        window.mdv.chartManager = this;

        this.listeners = {};
        this.infoAlerts = {};
        this.progressBars = {};
        this.setTheme(getPreferredColorScheme());
        //maybe better to stop listening once explicit option has been set
        //or to allow the user to explicitly say 'system default'
        listenPreferredColorScheme((t) => this.setTheme(t));
        /** !typed according to previous comments - but I was using it in a way that was working and doesn't match the comments...
         * each entry in dataSources will contain
         *  dataStore - the actual dataStore object (previously this comment erroneously referred to as 'dataSource')
         *  name - the name given to this data source
         *  menuBar the dom menu associated with this element
         *  contentDiv the div that the charts associated with the datastore will be added
         * @typedef {import("@/charts/charts").DataSource} DataSource
         * @type {DataSource[]}
         */
        this.dataSources = [];
        /** @type {{[k: string]: DataSource | undefined}} */
        this.dsIndex = {};
        this.columnsLoading = {};
        for (const d of dataSources) {
            const ds = {
                name: d.name,
                dataStore: new DataStore(d.size, d, dataLoader),
                link_to: d.link_to,
                index_link_to: d.index_link_to,
                color: d.color || themes[this.theme].background_color,
                column_link_to: d.column_link_to,
                links: d.links,
                custom: d.custom || {},
            };
            this.dataSources.push(ds);
            this.dsIndex[d.name] = ds;
            this._addDSListeners(ds);
            this.columnsLoading[d.name] = {};
        }
        if (listener) {
            this.addListener("_default", listener);
        }
        makeObservable(this, {
            // making these observable caused bugs - in particular, when we got to 'state_saved', ds.contentDiv was undefined
            // dataSources: observable,
            // dsIndex: observable,
            theme: observable,
            setTheme: action,
        });
        this.transactions = {};

        // View Manager
        this.viewManager = new ViewManager(config.initial_view, config.all_views);

        //set up container and top(main menu)
        /** @type {HTMLElement} */
        this.containerDiv =
            typeof div === "string" ? document.getElementById(div) : div;
        this.containerDiv.style.display = "flex";
        this.containerDiv.style.flexDirection = "column";

        /** @type {HTMLDivElement} */
        // classes: ["ciview-main-menu-bar"],
        this.menuBar = createEl("div", {}, this.containerDiv);

        createMdvPortal(MenuBarWrapper(), this.menuBar);

        this._setupThemeContextMenu();

        /** @type {HTMLDivElement} */
        this.contentDiv = createEl(
            "div",
            {
                styles: {
                    flex: "1 1 auto",
                    position: "relative",
                },
            },
            this.containerDiv,
        );
        this.contentDiv.classList.add("ciview-contentDiv");
        this.config = config;

        //!!! ChatMDV specific, but we really should be using websocket for other things
        //let's have a think. we may have websocket, but not chat... 
        //we *do* anticipate finding out about chat with a similar mechanism (flag as part of state.json)
        //but websocket should be ubiquitous
        if (config.websocket) {
            console.log('websocket is enabled');
            const fn = async () => {
                // previously, we were always calling connectIPC - but it was only relevant to earlier experiment with Unity
                // and we had disabled websocket on server.
                // started experimenting with socketio for chatMDV - mechanism is working, to an extent... 
                // but actually, REST is probably best for this (maybe a protocol agnostic abstraction).
                // try/catch doesn't help when it gets stuck in await...
                // console.warn('websocket is not currently supported but used as flag for chat experiment - will be fixed very soon')
                try {
                    const { socket, sendMessage } = await connectIPC(this);
                    console.log('connected to socketio');
                    this.ipc = { socket, sendMessage };
                } catch (error) {
                    console.error('Failed to connect to websocket', error);
                }
            };
            fn();
        }

        /** @type {GridStackManager} */
        this.gridStack = new GridStackManager(this);

        //each entry in charts will contain
        //  chart - the actual chart
        //  win - the popout window it is in (or null)
        //  dataSource - the data source associated with it
        /** @type {{[id: string]: {chart: import("./charts").Chart, win: Window | null, dataSource: import("./charts").DataSource}}} */
        this.charts = {};

        const c = this.config;
        c.chart_color = c.chart_color || "white";

        //load any files first

        this.dataLoader = dataLoader.function; // || async function defaultDataLoaderFunction() { console.warn(`ceci n'est pas une dataLoader`) };
        this.viewLoader = dataLoader.viewLoader;

        this.layoutMenus = {};
        this.isFullscreen = false;

        if (dataLoader.files) {
            this.filesToLoad = dataLoader.files.length;
            for (const item of dataLoader.files) {
                this.loadFile(item, () => {
                    this.filesToLoad--;
                    if (this.filesToLoad === 0) {
                        this._loadView(config, dataLoader, true);
                    }
                });
            }
        } else {
            this._loadView(config, dataLoader, true);
        }

        //add links
        for (const ds of this.dataSources) {
            const links = ds.dataStore.links;
            if (links) {
                for (const ods in links) {
                    const _ods = this.dsIndex[ods].dataStore;
                    if (!_ods) {
                        console.warn(
                            `datasource ${ds.name} has link to non-existent datasource ${ods}`,
                        );
                        return;
                    }
                    // pjt: could there be cases where we want more than one of a given type of link between two data sources?
                    // if so, we can support that by making these functions undertand how to interpret the links object appropriately.
                    if (links[ods].interactions) {
                        this._addInteractionLinks(
                            ds.dataStore,
                            _ods,
                            links[ods].interactions,
                        );
                    }
                    if (links[ods].valueset) {
                        this._addValuesetLink(
                            ds.dataStore,
                            _ods,
                            links[ods].valueset,
                        );
                    }
                    if (links[ods].column_pointer) {
                        console.warn(
                            "experimental column_pointer links not currently supported",
                        );
                        // this._addColumnPointerLink(ds.dataStore, _ods,links[ods].column_pointer);
                    }
                }
            }
        }
    }

    /**
     * Adds a link that will filter a column in one datasource based on the values in another:
     * when the values in the source column are filtered, the values in the destination column
     * will be filtered so that only those appearing in the selected values of the source column
     * will be displayed.
     */
    _addValuesetLink(ds, ods, link) {
        const valuesetFilter = ods.getDimension("valueset_dimension");
        const srcCol = ds.columnIndex[link.source_column];
        const destCol = ods.columnIndex[link.dest_column];
        const isTextLike = srcCol.values && destCol.values;
        Promise.all([
            this._getColumnsAsync(ds.name, [link.source_column]),
            this._getColumnsAsync(ods.name, [link.dest_column]),
        ]).then(() => {
            ds.addListener(
                `${ds.name}-${ods.name}_valueset`,
                async (type, data) => {
                    // if the current view doesn't use the ods, don't bother filtering
                    if (!this.viewData.dataSources[ods.name]) return;
                    if (type === "filtered") {
                        //`data` may be a (Range)Dimension (which has a `filterArray` property),
                        //or 'undefined', or a string like "all_removed", or who knows what else...
                        if (
                            !data ||
                            typeof data === "string" ||
                            !data.filterArray
                        ) {
                            valuesetFilter.removeFilter();
                            return;
                        }
                        await this._getColumnsAsync(ds.name, [
                            link.source_column,
                        ]);
                        const resultSet = await data.getValueSet(
                            link.source_column,
                        );
                        await this._getColumnsAsync(ods.name, [
                            link.dest_column,
                        ]);
                        valuesetFilter.filter(
                            "filterValueset",
                            [link.dest_column],
                            resultSet,
                        );
                    } else if (type === "data_highlighted") {
                        await this._getColumnsAsync(ds.name, [
                            link.source_column,
                        ]);
                        await this._getColumnsAsync(ods.name, [
                            link.dest_column,
                        ]);
                        const { indexes } = data;
                        if (isTextLike) {
                            // given a set of indexes into the source column, get the corresponding values
                            // the result should be a set of strings, and "filterValueset" will be responsible
                            // for filtering the destination column based on those strings
                            const resultSet = new Set(
                                indexes.map(
                                    (i) => srcCol.values[srcCol.data[i]],
                                ),
                            );
                            valuesetFilter.filter(
                                "filterValueset",
                                [link.dest_column],
                                resultSet,
                            );
                        } else {
                            valuesetFilter.filter(
                                "filterValueset",
                                [link.dest_column],
                                new Set(indexes),
                            );
                        }
                    }
                },
            );
        });
    }

    /**
     * @param {DataStore} ds
     * @param {DataStore} ods other data store
     * @param {object} link - information about what properties to target in ods
     * For instance, `{ source_column: <column_name>, target_chart: <chart_id>, target_property: <property_name> }`
     * more concretely: `{ source_column: "feature", target_chart: "quadrat_image_chart", target_property: "color_by" }`
     */
    _addColumnPointerLink(ds, ods, link) {
        const srcCol = ds.columnIndex[link.source_column];
        // do we expect the chart to be there when it might belong to a different view?
        // look it up in the listener for now, and just ignore/log if it's not there
        // const targetChart = this.charts[link.target_chart];
        // if (!targetChart) throw new Error(`No chart with id ${link.target_chart} found`);

        const isTextLike = srcCol.values; // perhaps 'unique' should also be considered text-like / valid?
        if (!isTextLike)
            throw new Error(
                "Only text-like columns are supported for column pointer links",
            );
        ds.addListener(
            `${ds.name}-${ods.name}_column_pointer`,
            async (type, data) => {
                //don't think we need mobx action as autoObservable is used
                if (type === "data_highlighted") {
                    const targetChart = this.charts[link.target_chart]?.chart;
                    if (!targetChart) {
                        console.log(
                            `No chart with id ${link.target_chart} found - bypassing column pointer link`,
                        );
                        return;
                    }
                    await this._getColumnsAsync(ds.name, [link.source_column]);
                    // data probably looks something like `{indexes: Array(1), source: undefined, data: null, dataStore: DataStore}`
                    const { indexes } = data;
                    const newValue = srcCol.values[srcCol.data[indexes[0]]];
                    targetChart.config[link.target_property] = newValue;
                    targetChart.setTitle(newValue); //could be optional - or use some kind of template?
                    if (link.target_property === "color_by") {
                        targetChart.colorByColumn(newValue);
                    }
                }
            },
        );
    }

    _addInteractionLinks(ds, ods, links) {
        const interactionFilter = ods.getDimension("category_dimension");
        const icols = links.interaction_columns;
        const c1 = icols[0];
        const c2 = icols[1];
        const pc = links.pivot_column;
        //sync the colors of the matching columns
        ds.syncColumnColors.push({
            dataSource: ods.name,
            columns: [
                { link_to: icols[2], col: icols[0] },
                { link_to: icols[2], col: icols[1] },
                { link_to: pc, col: pc },
            ],
        });
        ds.addListener(
            `${ds.name}-${ods.name}_interaction`,
            async (type, data) => {
                if (type === "data_highlighted") {
                    //get the two interacting items plus pivot
                    await this._getColumnsAsync(ds.name, [c1, c2, pc]);
                    const info = ds.getRowAsObject(data.indexes[0], [
                        c1,
                        c2,
                        pc,
                    ]);
                    //get pivot from the other datasource
                    await this._getColumnsAsync(ods.name, [icols[2]]);
                    //filter the two interacting items
                    //would be nice if this could be maybe async / lazy / different ways of composing filters?

                    const args = [info[c1], info[c2]];
                    // @ts-ignore should revisit this at some point
                    args.operand = "and"; //force multitext to use "and"

                    interactionFilter.filter(
                        "filterCategories",
                        [icols[2]],
                        args,
                    );
                    //show the region if not already displayed
                    if (
                        links.is_single_region &&
                        ods.regions &&
                        this.dsPanes[ods.name]
                    ) {
                        if (
                            !this.getAllCharts(ods.name).find(
                                (x) => x.config.region === info[pc],
                            )
                        ) {
                            const conf = {
                                type: "image_scatter_plot",
                            };
                            //add the default parameters
                            BaseChart.types["image_scatter_plot"].init(
                                conf,
                                ods,
                                { region: info[pc] },
                            );
                            //add the chart
                            this.addChart(ods.name, conf);
                        }
                    }
                }
            },
        );
    }

    _setUpChangeLayoutMenu(ds) {
        this.layoutMenus[ds.name] = new ContextMenu(() => {
            const lo = this.viewData.dataSources[ds.name].layout || "absolute";
            return [
                {
                    text: "Absolute",
                    ghosted: lo === "absolute",
                    func: () => this.changeLayout("absolute", ds),
                },
                {
                    text: "Grid Stack",
                    ghosted: lo === "gridstack",
                    func: () => this.changeLayout("gridstack", ds),
                },
            ];
        });
    }

    changeLayout(type, ds) {
        const view = this.viewData.dataSources[ds.name];
        const current = view.layout || "absolute";
        if (type === current) {
            return;
        }
        //remove existing layouts on charts
        if (current === "gridstack") {
            this.gridStack.destroy(ds);
        } else if (current === "absolute") {
            this.getAllCharts(ds.name).forEach((x) => {
                const div = x.getDiv();
                removeResizable(div);
                removeDraggable(div);
            });
        }
        
        view.layout = type;
        //will add the appropriate layout depending on the current layout type
        this.getAllCharts(ds.name).forEach((x) => {
            //the chart is popped out so not subject to the layout manager
            if (x.__doc__ !== document){
                return;
            }
            this._makeChartRD(x, ds)
        });
    }

    _setupThemeContextMenu() {
        this.themeMenu = new ContextMenu(() => {
            const mItems = [];
            for (const t in themes) {
                mItems.push(this.__getMenuItem(t));
            }
            return mItems;
        });
    }

    __getMenuItem(theme) {
        return {
            text: theme,
            ghosted: this.theme === theme,
            func: () => this.setTheme(theme),
        };
    }

    setTheme(theme) {
        //thinking about doing everything with css
        // there could be graphics rendering of other sorts as well...
        // nothing I can see at the moment that responds to theme.
        this.theme = theme;
        document.getElementsByTagName("html")[0].className = theme;
        //only chart this is required for is the genome browser
        //it uses canvas and thus has to redraw the canavas which is
        //just a png so won't be effected by css changes
        if (!this.charts) {
            return;
        }
        for (const ch in this.charts) {
            const c = this.charts[ch];
            if (c.chart.themeChanged) {
                c.chart.themeChanged();
            }
        }
    }

    //sync color columns
    _sync_colors(columns, from, to) {
        for (const item of columns) {
            const from_col = from.columnIndex[item.link_to];
            const to_col = to.columnIndex[item.col];
            const newColors = new Array(to_col.values);
            const colors = from.getColumnColors(item.link_to);
            //"cannot read properties of undefined (reading 'length')"
            //seems to be related to from_col being type 'integer' when it should be 'text'
            if (!from_col.values || !to_col.values) {
                //todo display a 'toast' error...
                console.error(
                    `failed to _sync_colors '${item.link_to}', '${item.col}' - expected 'text' or similar columns`,
                );
                continue;
            }
            for (let i = 0; i < from_col.values.length; i++) {
                const val = from_col.values[i];
                const index = to_col.values.indexOf(val);
                if (index !== -1) {
                    newColors[index] = colors[i];
                }
            }
            to_col.colors = newColors;
        }
    }

    _initiateOffsets(dataSource) {
        const ds = dataSource.dataStore;
        const o = ds.offsets;
        const p = o.param;
        //need to make sure all columns are loaded
        const cols = [p[0], p[1], o.groups];
        if (o.background_filter) {
            cols.push(o.background_filter);
        }
        this._getColumnsThen(dataSource.name, cols, () => {
            ds.initiateOffsets();
            //update x,y offsets
            ds.updateColumnOffsets();
            //update rotation and update
            ds.updateColumnOffsets(null, true, true);
        });
    }

    //load the view metadata or use initialCharts then call _init to load the view
    _loadView(config, dataLoader, firstTime = false) {
        //load view, then initialize
        if (config.all_views) {
            const currentView = config.initial_view || config.all_views[0];
            this.viewManager.setView(currentView);
            dataLoader.viewLoader(currentView).then(async (data) => {
                try {
                    await this._init(data, firstTime);
                    if (currentView) {
                        const state = this.getState();
                        if (!state.view?.viewImage) {
                            // todo: update the error handler function in future
                            await this.viewManager.saveView((state) => {
                                console.log("Error occurred: ", state.chartErrors);
                                return false;
                            });
                        }
                    }
                } catch (error) {
                    console.error("Error during view initialization:", error);
                    // Consider adding user-facing error handling here
                }
            });
        }
        //only one view hard coded in config
        //! This else block is not called, but if it is called at some point, make sure the state save works properly
        else {
            this._init(config.only_view, firstTime)
            .then(async () => {
                    const state = this.getState();
                    if (!state.view?.viewImage) {
                        // todo: update the error handler function in future
                        await this.viewManager.saveView((state) => {
                            console.log("Error occurred: ", state.chartErrors);
                            return false;
                        });
                    }
            });
        }
    }

    /**
     * Caution: doesn't return a 'DataSource', but a 'DataStore' (which is a property of a 'DataSource')
     * @param {string} name
     * @returns {DataStore}
     */
    getDataSource(name) {
        return this.dsIndex[name].dataStore;
    }

    async _init(view, firstTime = false) {
        //no initial view just make one with all available
        //DataSources but no charts
        if (!view) {
            const dts = {};
            for (const ds in this.dsIndex) {
                dts[ds] = { layout: "gridstack" };
            }
            this.viewData = { dataSources: dts, initialCharts: {} };
        } else {
            //legacy data (which only has initialCharts)
            //need to add which DataSources to display
            if (!view.dataSources) {
                view.dataSources = {};
                for (const ds in view.initialCharts) {
                    view.dataSources[ds] = {};
                }
            }
            this.viewData = view;
        }

        for (const ds of this.dataSources) {
            ds.contentDiv = undefined;
            ds.menuBar = undefined;
        }
        const dsToView = Object.keys(this.viewData.dataSources);

        let widths = [];
        for (const ds of dsToView) {
            const w = this.viewData.dataSources[ds].panelWidth;
            if (!w) {
                widths = null;
                break;
            }
            widths.push(w);
        }
        this.dsPanes = {};

        //add all the appropriate panes (one per datasource)
        const panes = splitPane(this.contentDiv, {
            number: dsToView.length,
            sizes: widths,
        });
        //thinking about adding an optional index to the view so we can control the order
        for (let n = 0; n < dsToView.length; n++) {
            const p = panes[n];

            p.style.display = "flex";
            p.style.flexDirection = "column";
            const ds = this.dsIndex[dsToView[n]];
            this.dsPanes[ds.name] = p;
            this.columnsLoading[ds.name] = {};
            ds.charts = [];
            ds.menuBar = createEl(
                "div",
                {
                    classes: ["ciview-menu-bar"],
                },
                p,
            );
            const d = createEl(
                "div",
                {
                    styles: {
                        flex: "1 1 auto",
                        position: "relative",
                        overflow: "auto",
                        height: "100px",
                    },
                },
                p,
            );
            this._setUpMenu(ds);
            // might move styles from here into .css
            ds.contentDiv = createEl(
                "div",
                {
                    styles: {
                        //flex:"1 1 auto",
                        position: "relative",
                        height: "100%",

                        //overflow:"auto",
                        // background:col
                    },
                },
                d,
            );
            ds.contentDiv.classList.add("ciview-contentDiv");
            this._setUpChangeLayoutMenu(ds);
            //need to add
        }

        //any first time initiation
        if (firstTime) {
            for (const d of this.dataSources) {
                const ds = d.dataStore;
                //initiate offsets if any
                if (ds.offsets) {
                    this._initiateOffsets(d);
                }
                //sync any columns
                //phasing out
                //todo ^^ PJT this 'phasing out' comment was written a year ago as of 24-01-15...
                //and this code can cause problems.
                if (d.column_link_to) {
                    this._sync_colors(
                        d.column_link_to.columns,
                        this.dsIndex[d.column_link_to.dataSource].dataStore,
                        ds,
                    );
                }
                for (const scc of ds.syncColumnColors) {
                    try {
                        this._sync_colors(
                            scc.columns,
                            this.dsIndex[scc.dataSource].dataStore,
                            ds,
                        );
                    } catch (error) {
                        console.warn(
                            `error syncing colors for '${ds.name}' and '${scc.dataSource}'`,
                        );
                    }
                }
            }
        }

        //need to create a set to create track of
        //charts loaded
        const charts = view ? view.initialCharts || {} : {};
        this._toLoadCharts = new Set();
        for (const ds in charts) {
            for (const ch of charts[ds]) {
                this._toLoadCharts.add(ch);
            }
        }
        //nothing to load - call any listeners
        if (this._toLoadCharts.size === 0) {
            this._toLoadCharts = undefined;
            if (this.viewManager.current_view === undefined) {
                if (this.dataSources.length === 0) new FileUploadDialogReact();
                else {
                    // todo - separate out view code, with a solid model and start making some nice ui etc...
                    this.showAddViewDialog();
                }
            } else {
                this._callListeners("view_loaded", this.viewManager.current_view);
            }
        }
        //add charts - any columns will be added dynamically
        this._inInit = true;
        const chartPromises = [];
        for (const ds in charts) {
            for (const ch of charts[ds]) {
                chartPromises.push(this.addChart(ds, ch));
            }
        }
        console.log('before await Promise.all(chartPromises);');
        await Promise.all(chartPromises);
        // it prints some "decorated loadColumnData method called" after this line...
        // is it ever with a value that would cause the chart.getConfig() to have a different result?
        // We are relying on the config being stable/settled when these promises resolve.
        console.log('after await Promise.all(chartPromises);');

        // Restore highlighted indices from view data (e.g. after loading or switching view)
        for (const dsName of Object.keys(this.viewData.dataSources || {})) {
            const spec = this.viewData.dataSources[dsName];
            const highlight = spec?.highlight;
            const dstore = this.dsIndex[dsName]?.dataStore;
            if (!dstore?.dataHighlighted) continue;
            if (Array.isArray(highlight) && highlight.length > 0) {
                dstore.dataHighlighted(highlight, null);
            } else {
                dstore.dataHighlighted([], null);
            }
        }

        //this could be a time to _callListeners("view_loaded",this.currentView)
        //but I'm not going to interfere with the current logic
        this._inInit = false;
    }

    _addDSListeners(ds) {
        ds.dataStore.addListener("l1", (type, data) => {
            if (type === "column_removed") {
                this._columnRemoved(ds, data);
            } else if (type === "data_highlighted") {
                data.dataStore = ds.dataStore;
                this._callListeners(type, data);
            } else if (type === "filtered") {
                if (!this.progressBars[ds.name]) {
                    return;
                }
                const n1 = ds.dataStore.size;
                const n2 = ds.dataStore.filterSize;
                this.progressBars[ds.name].setValue(n2);
                this.progressBars[ds.name].setText(n2);
                this._callListeners(type, data);
            }
        });
    }

    showAddViewDialog() {
        new ViewDialogWrapper("add");
    }

    /**
     * @param {()=>void=} action
     * @param {string=} content
     */
    showSaveViewDialog(action, content) {
        new ViewDialogWrapper("save", action, content);
    }

    showDeleteViewDialog() {
        new ViewDialogWrapper("delete");
    }

    changeView(view) {
        this.viewManager.checkAndChangeView(view);
    }

    _columnRemoved(ds, col) {
        const ids_to_delete = [];
        for (const id in this.charts) {
            const info = this.charts[id];
            if (info.dataSource === ds) {
                const ch = info.chart;
                const div = ch.getDiv();
                const del = ch.onColumnRemoved(col);
                if (del) {
                    div.remove(false);
                    ids_to_delete.push(id);
                    this._removeLinks(ch);
                    this._callListeners("chart_removed", ch);
                }
            }
        }
        //onColumnRemoved will remove the chart if it contains
        //data from the column, it will also remove the filter,
        //but not call any listeners
        if (ids_to_delete.length > 0) {
            ds.dataStore._callListeners("filtered");
        }
        for (const id of ids_to_delete) {
            delete this.charts[id];
        }
    }

    _getColumnsRequiredForChart(config) {
        const set = new Set();
        let p = config.param;

        if (!p) {
            //no 'parameters',
            //*but there could be other config entries / methods that refer to columns*
            // return [];
        } else if (typeof p === "string") {
            // pretty sure there's nothing in BaseChart.types that would get here - single param is ["string"]
            // but the LLM might still generate a config with a single string, or an old config might have one
            console.error(`Unexpected param string '${config.param}' for ${config.name} - expected array`);
            set.add(p);
            p = [p];
            config.param = p;
        } else {
            for (const i of p) {
                set.add(i);
            }
        }
        if (config.color_by) {
            // LegacyColorBy - can we nip it in the bud so we can use narrower types elsehwere?
            if (config.color_by.column) {
                set.add(config.color_by.column.field);
            } else {
                set.add(config.color_by);
            }
        }
        if (config.tooltip) {
            if (config.tooltip.column) {
                if (Array.isArray(config.tooltip.column)) {
                    for (const i of config.tooltip.column) {
                        set.add(i);
                    }
                } else {
                    set.add(config.tooltip.column);
                }
            }
        }
        if (config.background_filter) {
            set.add(config.background_filter.column);
        }

        //are there any config entries that refer to column(s)
        const t = BaseChart.types[config.type];
        // this will probably be handled differently in the (near) future
        // maybe later we could use a decorator to add this to the config?
        // (not sure it'd be necessary - as long as they have a consistent interface)
        if (t.configEntriesUsingColumns) {
            t.configEntriesUsingColumns.forEach((x) => {
                let e = config[x];
                if (e) {
                    e = Array.isArray(e) ? e : [e];
                    for (const i of e) {
                        set.add(i);
                    }
                }
            });
        }
        return [...set];
    }

    /**
     * Loads data from a remote file -the file must have headers (keys in the
     * case of json) which which match a columns field/id
     * @param {object} info A config describing the file -
     * @param {string} info.type - either csv,tsv ot json
     * @param {string} info.dataSource - the name of the datasource to load the data into
     * @param {string} info.url  - the url of the file
     * @param {function} [callback]  - a function to run once the data has loaded
     */
    loadFile(info, callback) {
        const meths = { csv: csv, json: json, tsv: tsv };
        const iid = this.createInfoAlert("loading file", { spinner: true });
        meths[info.type](info.url).then((data) => {
            const cols = {};
            const dataSource = info.dataSource;
            const ds = this.dsIndex[dataSource].dataStore;
            const all_cols = ds.getAllColumns();
            //which columns are present in the datastore
            for (const c of data.columns) {
                if (ds.columnIndex[c]) {
                    cols[c] = [];
                }
            }
            for (let i = 0; i < data.length; i++) {
                const row = data[i];
                if (i + (1 % 100) === 0) {
                    this.updateInfoAlert(
                        iid,
                        `processed ${i}/${data.length} rows`,
                    );
                }
                for (const col in cols) {
                    cols[col].push(row[col]);
                }
            }
            let proc = 0;
            for (const col in cols) {
                ds.setColumnData(col, cols[col]);
                proc++;
                this.updateInfoAlert(
                    iid,
                    `processed ${proc}/${all_cols.length} columns`,
                );
            }
            this.updateInfoAlert(iid, "complete", { duration: 2000 });
            if (callback) {
                callback();
            }
        });
    }

    _getUpdatedColumns(dataStore) {
        const dc = dataStore.dirtyColumns;
        const rv = {
            columns: [],
            added: [],
            removed: [],
            colors_changed: [],
        };
        for (const c in dc.added) {
            const td = getMd(c);
            rv.columns.push(td);
            rv.added.push(c);
        }
        for (const r in dc.removed) {
            rv.removed.push(r);
        }

        for (const c in dc.data_changed) {
            if (!rv.columns[c]) {
                const td = getMd(c);
                rv.columns.push(td);
            }
        }

        for (const cc in dc.colors_changed) {
            rv.colors_changed.push({
                column: cc,
                colors: dataStore.columnIndex[cc].colors,
            });
        }

        return rv;

        function getMd(c) {
            const cl = dataStore.columnIndex[c];
            const md = {
                values: cl.values,
                datatype: cl.datatype,
                name: cl.name,
                editable: true,
                field: cl.field,
            };
            const numRows = dataStore.size;
            
            // Add stringLength to metadata for unique columns (required by server)
            if (cl.datatype === "unique" && cl.stringLength) {
                md.stringLength = cl.stringLength;
            }
            
            let arr;
            if (cl.datatype === "unique") {
                // For unique columns, convert Uint8Array to array of strings (expected by server)
                const textDecoder = new TextDecoder();
                const stringLength = cl.stringLength;

                if (!stringLength || typeof stringLength !== "number" || stringLength <= 0) {
                    console.error(
                        `Column ${c} has invalid or missing stringLength: ${stringLength}.`
                    );
                    // Fallback as empty values to keep it consistent with other columns
                    return { metadata: md, data: new Array(numRows).fill("") };
                }

                arr = new Array(numRows);
                
                for (let i = 0; i < numRows; i++) {
                    const baseIndex = i * stringLength;

                    if (!cl.data || baseIndex + stringLength > cl.data.length) {
                        throw new Error(
                            `Invalid unique-column buffer for '${c}' at row ${i} (stringLength=${stringLength}).`
                        );
                    }

                    const rowBytes = cl.data.subarray(baseIndex, baseIndex + stringLength);
                    const decoded = textDecoder.decode(rowBytes);
                    // Remove null padding characters
                    arr[i] = decoded.replace(/\0+$/, "");
                }
            } else {
                // For other datatypes, get the values from data array
                arr = new Array(cl.data.length);
                for (let i = 0; i < cl.data.length; i++) {
                    arr[i] = cl.data[i];
                }
            }
            
            return { metadata: md, data: arr };
        }
    }

    saveState() {
        const state = this.getState();
        this._callListeners("state_saved", state);
    }

    getState() {
        const initialCharts = {};
        const updatedColumns = {};
        const metadata = {};
        const chartErrors = [];
        // const twidth = this.contentDiv.offsetWidth;
        for (const ds of this.dataSources) {
            if (ds.contentDiv) {
                initialCharts[ds.name] = [];
                const w = this.dsPanes[ds.name].style.width;
                const re2 = /calc\((.+)\%.+/;
                this.viewData.dataSources[ds.name].panelWidth =
                    Number.parseFloat(w.match(re2)[1]);
            }

            updatedColumns[ds.name] = this._getUpdatedColumns(ds.dataStore);
            const dstore = ds.dataStore;

            if (dstore.dirtyMetadata.size !== 0) {
                metadata[ds.name] = {};
                for (const param of dstore.dirtyMetadata) {
                    metadata[ds.name][param] = dstore[param];
                }
            }
        }
        for (const chid in this.charts) {
            const chInfo = this.charts[chid];
            const chart = chInfo.chart;
            try {
                // as of now, some charts throw errors when calling getConfig()
                // in particular, if they haven't finished loading properly
                const config = chart.getConfig();
                const div = chart.getDiv();
                const d = this.viewData.dataSources[chInfo.dataSource.name];
                if (d.layout === "gridstack") {
                   //this is handled by gridstack now
                } else {
                    config.position = [div.offsetLeft, div.offsetTop];
                }
                initialCharts[chInfo.dataSource.name].push(config);
            } catch (error) {
                console.error(
                    `chart ${chid} failed to getConfig - ${error.message}`,
                );
                chartErrors.push(error);
            }
        }

        const view = JSON.parse(JSON.stringify(this.viewData));
        view.initialCharts = initialCharts;
        for (const ds of this.dataSources) {
            const h = ds.dataStore.getHighlightedData?.();
            if (!view.dataSources[ds.name]) view.dataSources[ds.name] = {};
            if (Array.isArray(h) && h.length > 0) {
                view.dataSources[ds.name].highlight = [...h];
            } else {
                // biome-ignore lint/performance/noDelete: setting it to undefined would mean there would still be the key, not a perf issue here.
                delete view.dataSources[ds.name].highlight;
            }
        }
        const all_views = this.viewManager.all_views ? this.viewManager.all_views : null;
        return {
            view: view,
            currentView: this.viewManager.current_view,
            all_views: all_views,
            updatedColumns: updatedColumns,
            metadata: metadata,
            //untested - what happens if we include this (including in what we save to server?)
            //could be useful... but I'm leaving it out now in case of unexpected issues.
            chartErrors, 
        };
    }

    setAllColumnsClean() {
        for (const ds of this.dataSources) {
            ds.dataStore.setAllColumnsClean();
        }
    }

    /** Displays a dialog
     * @param {Object} config extra settings
     */

    showCustomDialog(config) {
        new CustomDialog(config);
    }

    /**Adds a menu icon to either the main menubar or a datasource menubar
     * @param {string} dataSource The name of data source or _main if adding
     * an icon to the main (top) toolbar
     * @param {string} icon The class name(s) of the icon
     * @param {string} text Text that will be displayed in a tooltip
     * @param {function} func The function that will be called when the icon is pressed
     * @param {boolean} right If true (and `dataSource = "_main"`) the icon will be added to the right of the menu bar
     */
    addMenuIcon(dataSource, icon, text, func, right = false) {
        const pos = dataSource === "_main" ? "bottom-right" : "bottom";
        const el =
            dataSource === "_main"
                ? right ? this.rightMenuBar : this.leftMenuBar
                : this.dsIndex[dataSource].menuBar;
        return createMenuIcon(
            icon,
            {
                tooltip: {
                    text: text,
                    position: pos,
                },
                func: func,
            },
            el,
        );
    }

    createInfoAlert(msg, config = {}) {
        const id = getRandomString();
        const len = Object.keys(this.infoAlerts).length;
        config.type = config.type || "info";
        const div = createEl(
            "div",
            {
                classes: ["ciview-info-alert", `ciview-alert-${config.type}`],
                styles: {
                    right: "10px",
                    top: `${50 + len * 40}px`,
                },
            },
            this.containerDiv,
        );
        let spinner = null;
        const text = createEl("span", { text: msg }, div);
        if (config.spinner) {
            spinner = createEl(
                "i",
                {
                    classes: [
                        "fas",
                        "fa-spinner",
                        "fa-spin",
                        "ciview-info-alert-spin",
                    ],
                },
                div,
            );
        }
        this.infoAlerts[id] = {
            div: div,
            text: text,
            spinner: spinner,
            type: config.type,
        };
        if (config.duration) {
            this.removeInfoAlert(id, config.duration);
        }
        return id;
    }

    updateInfoAlert(id, msg, config = {}) {
        const al = this.infoAlerts[id];
        if (al) {
            if (config.type && al.type !== config.type) {
                al.div.classList.remove(`ciview-alert-${al.type}`);
                al.div.classList.add(`ciview-alert-${config.type}`);
                al.type = config.type;
            }
            al.text.textContent = msg;
            if (config.duration) {
                this.removeInfoAlert(id, config.duration);
            }
        }
    }

    removeInfoAlert(id, delay = 2000) {
        const spinner = this.infoAlerts[id].spinner;
        if (spinner) {
            spinner.remove();
        }
        setTimeout(() => {
            if (!this.infoAlerts[id]) return; //PJT allow for clearing list.
            this.infoAlerts[id].div.remove();
            delete this.infoAlerts[id];
            let top = 50;
            for (const i in this.infoAlerts) {
                this.infoAlerts[i].div.style.top = `${top}px`;
                top += 40;
            }
        }, delay);
    }
    clearInfoAlerts() {
        for (const i in this.infoAlerts) {
            this.infoAlerts[i].div.remove();
        }
        this.infoAlerts = {};
    }

    /**
     * Loads data for specified columns into the appropriate dataStore
     * @param {string[]} columns An array of column fields/ids
     * @param {string} dataSource The name of the dataSource
     * @param {function} callback A function which will be run once all the
     * columns are loaded, with any failed columns as an argument (although it's not clear that the underlying code actually does this,
     * or that any code that calls this function actually uses the argument)
     * @param {number} [split=10]  the number of columns to send with each request
     * @param {number} [threads=2]  the number of concurrent requests
     */
    loadColumnSet(columns, dataSource, callback, split = 10, threads = 2) {
        const ds = this.getDataSource(dataSource);
        // nb, if any items are in columnsLoading, we don't filter them here
        // because we'd have to think of how to listen for their completion before calling the callback
        // !! could be an issue if a column has been edited by another user and this method is supposed to reload it?
        columns = columns.filter((x) => !ds.columnsWithData.includes(x));
        if (columns.length === 0) {
            callback(); //should there be any args? hard to tell without types
            return;
        }
        const id = getRandomString();
        const lc = this.config.dataloading || {};
        split = lc.split || 10;
        threads = lc.threads || 2;
        this.transactions[id] = {
            callback: callback,
            columns: [],
            totalColumns: columns.length,
            failedColumns: [],
            nextColumn: 0,
            columnsLoaded: 0,
            id: id,
        };
        let col_list = [];
        const t = this.transactions[id];
        for (const col of columns) {
            this.columnsLoading[dataSource][col] = true;
            col_list.push(col);
            if (col_list.length === split) {
                t.columns.push(col_list);
                col_list = [];
            }
        }
        if (col_list.length !== 0) {
            t.columns.push(col_list);
            col_list = [];
        }
        t.alertID = this.createInfoAlert(
            `Loading Columns:0/${columns.length}`,
            { spinner: true },
        );
        const max = Math.min(t.columns.length, threads);

        for (let n = 0; n < max; n++) {
            this._loadColumnData(t, dataSource);
        }
    }
    loadColumnSetAsync(columns, dataSource, split = 10, threads = 2) {
        return new Promise((resolve, reject) => {
            this.loadColumnSet(columns, dataSource, resolve, split, threads);
        });
    }

    _loadColumnData(trans, dataSource) {
        const dataStore = this.dsIndex[dataSource].dataStore;

        const col_list = trans.columns[trans.nextColumn++];
        const columns = [];
        for (const col of col_list) {
            columns.push(dataStore.getColumnInfo(col));
        }
        //float32 columns need to be at the beginning of the byte stream
        //as you can't create an array from  an arry buffer starting at
        //a byte position not divisible by 4
        columns.sort((a, b) => {
            return column_orders[a.datatype] - column_orders[b.datatype];
        });

        //"this.dataLoader is not a function" with e.g. "cell_types"
        this.dataLoader(columns, dataSource, dataStore.size)
            .then((resp) => {
                for (const col of resp) {
                    dataStore.setColumnData(col.field, col.data);
                }
                trans.columnsLoaded++;
            })
            .catch((error) => {
                //! this is an error that is not being handled properly... what happens to trans.failedColumns?
                //! they do get passed to a callback... what does that callback do?
                console.log(error);
                trans.columnsLoaded++;
                trans.failedColumns.push(columns);
            })
            .finally(() => {
                const total = trans.columns.length;
                const loaded = trans.columnsLoaded;
                let all_loaded = loaded * col_list.length;
                for (const col of col_list) {
                    delete this.columnsLoading[dataSource][col];
                }
                all_loaded =
                    all_loaded > trans.totalColumns
                        ? trans.totalColumns
                        : all_loaded;
                if (trans.failedColumns.length > 0) {
                    this.updateInfoAlert(
                        trans.alertID,
                        `Failed to load ${trans.failedColumns.length} columns`,
                        { type: "danger" },
                    );
                    // return;
                } else {
                    this.updateInfoAlert(
                        trans.alertID,
                        `Loading Columns:${all_loaded}/${trans.totalColumns}`,
                    );
                }
                if (loaded >= total) {
                    this.updateInfoAlert(
                        trans.alertID,
                        `Loaded ${total} column${total === 1 ? "" : "s"}`,
                        { duration: 2000 },
                    );
                    trans.callback(trans.failedColumns);
                    delete this.transactions[trans.id];
                }
                if (trans.nextColumn < total) {
                    this._loadColumnData(trans, dataSource);
                }
            });
    }

    /**
     * @param {{dataStore: DataStore}} ds
     */
    _setUpMenu(ds) {
        const dataStore = ds.dataStore;
        createMenuIcon(
            "fas fa-chart-bar",
            {
                tooltip: {
                    text: "Add Chart",
                    position: "bottom-right",
                },
                func: () => {
                    // new AddChartDialog(ds, (config) =>
                    //     this.addChart(ds.name, config, true),
                    // );
                    // if (import.meta.env.DEV) new BaseDialog.experiment["AddChartDialogReact"](dataStore);
                    new AddChartDialogReact(dataStore);
                },
            },
            ds.menuBar,
        );

        createMenuIcon(
            "fas fa-sync-alt",
            {
                tooltip: {
                    text: "Reset All Filters",
                    position: "bottom-right",
                },
                func: () => {
                    dataStore.removeAllFilters();
                },
            },
            ds.menuBar,
        );
        createMenuIcon(
            "fas fa-palette",
            {
                tooltip: {
                    text: "Change Color Scheme",
                    position: "bottom-right",
                },
                func: () => {
                    try {
                        new ColorChooser(this, ds);
                    } catch (error) {
                        console.error("error making ColorChooser", error);
                        this.createInfoAlert("Error making color chooser", {
                            type: "warning",
                            duration: 2000,
                        });
                    }
                },
            },
            ds.menuBar,
        );
        createMenuIcon(
            "fas fa-th",
            {
                tooltip: {
                    text: "Change layout",
                    position: "bottom-right",
                },
                func: (e) => {
                    this.layoutMenus[ds.name].show(e);
                },
            },
            ds.menuBar,
        );
        //previously we shared a TagModel between invocations of the annotation dialog
        //but they should be able to change columns - not sure if the sharing had any benefits
        this.addMenuIcon(ds.name, "fas fa-tags", "Tag Annotation", () => {
            //todo - check whether we have a reason for hacky import here
            new BaseDialog.experiment["AnnotationDialogReact"](ds.dataStore);
        });

        const idiv = createEl(
            "div",
            {
                styles: {
                    float: "right",
                    lineHeight: "1.0",
                },
            },
            ds.menuBar,
        );
        createEl(
            "span",
            {
                text: ds.name,
                styles: {
                    verticalAlign: "top",
                    fontSize: "16px",
                    marginRight: "4px",
                },
            },
            idiv,
        );
        const size = ds.dataStore.size;
        ds.filterBar = createEl(
            "progress",
            {
                value: size,
            },
            ds.menBar,
        );
        const pb = createEl(
            "div",
            {
                styles: {
                    width: "100px",
                    display: "inline-block",
                    marginTop: "2px",
                },
            },
            idiv,
        );
        const pbConf = {
            max: size,
            value: size,
            text: `${size}`,
        };
        this.progressBars[ds.name] = new MDVProgress(pb, pbConf);
        this._addFullscreenIcon(ds);
    }

    _addFullscreenIcon(ds) {
        const iconElement = createMenuIcon(
            "fas fa-expand",
            {
                tooltip: {
                    text: "Full Screen",
                    position: "bottom-right",
                },
                func: async () => {
                    //nb, not sure best way to access the actual div I want here
                    //this could easily break if the layout structure changes
                    // ds.contentDiv.parentElement.requestFullscreen();
                    try {
                        if (!this.isFullscreen) {
                            await ds.menuBar.parentElement.requestFullscreen();
                        } else {
                            await document.exitFullscreen();
                        }
                    } catch (error) {
                        console.error("fullscreen error caused: ", error);
                    }
                },
            },
            ds.menuBar,
        );

        document.addEventListener("fullscreenchange", () => {
            const iconEl = iconElement.querySelector("i");
            if (document.fullscreenElement) {
                if (ds.menuBar.parentElement === document.fullscreenElement) {
                    if (iconEl) {
                        iconEl.classList.remove("fa-expand");
                        iconEl.classList.add("fa-compress");
                    }
                    iconElement.setAttribute("aria-label", "Exit Full Screen");
                    this.isFullscreen = true;
                }
             } else {
                if (iconEl) {
                    iconEl.classList.remove("fa-compress");
                    iconEl.classList.add("fa-expand");
                }
                iconElement.setAttribute("aria-label", "Full Screen");
                this.isFullscreen = false;
            }
        });

    }

    /**
     * Adds a listener to the ChartManager.
     *
     * Note that the signature is different from what might be expected (this is an MDV pattern that may be reviewed in future):
     * The first argument is a string that identifies the listener (not an event type to listen for).
     * Subsequent calls to addListener with the same id will overwrite the previous listener (which could lead to
     * unexpected behaviour, for example if a corresponding `removeListener()` call is made).
     * The second argument is a function that will be called when the listener is triggered - all registered listeners will be called
     * for any event.
     * @param {string} id
     * @param {function} func Callback function.
     * First argument is the type of event (`"state_saved" | "view_loaded" | "chart_added" | "chart_removed"`),
     * second is the `ChartManager` instance, third is the data.
     */
    addListener(id, func) {
        this.listeners[id] = func;
    }

    /**
     * Removes a listener from the ChartManager.
     * @param {string} id
     */
    removeListener(id) {
        delete this.listeners[id];
    }
    _callListeners(type, data) {
        for (const id in this.listeners) {
            this.listeners[id](type, this, data);
        }
    }

    /**
     * Adds a chart to the app. Returns asyncronously when needed columns have loaded and the chart has been added
     * (will reject if there is an error)
     * @param {string} dataSource The name of the chart's data source
     * @param {any} config The chart's config
     * @param {boolean} [notify=false] If true any listeners will be informed that
     * a chart has been loaded
     */
    async addChart(dataSource, config, notify = false) {
        if (!BaseChart.types[config.type]) {
            this.createInfoAlert(
                `Tried to add unknown chart type '${config.type}'`,
                { type: "danger", duration: 2000 },
            );
            throw `Unknown chart type ${config.type}`;
        }
        if (typeof config.param === "string") {
            console.error(`Unexpected param string '${config.param}' for ${config.name} - expected array`);
            config.param = [config.param];
        }
        //check if columns need loading
        const neededCols = this._getColumnsRequiredForChart(config);
        //check which columns need loading
        if (config.location) {
            const l = config.location;
            const b = 5;
            config.size = [
                l.width * 90 + l.width * b - b,
                l.height * 40 + l.height * b - b,
            ];
            config.position = [
                (l.x + 1) * b + l.x * 90,
                (l.y + 1) * b + l.y * 40,
            ];
        }
        //**convert legacy data***********
        const ds = this.dsIndex[dataSource];
        const { width, height, left, top } = positionChart(ds, config);

        const t = themes[this.theme];
        // PJT may want different behaviour for gridstack
        //MJS this is very messy - create divs in (hopefully) the right location and add chart when data loaded
        //ideally create the chart with a waiting icon and  update it when the data has loaded
        //However, no way of creating charts at the moment without data - charts need separate init method?
        const div = createEl(
            "div",
            {
                styles: {
                    position: "absolute",
                    width: `${width}px`,
                    height: `${height}px`,
                    left: `${left}px`,
                    top: `${top}px`,
                    background: "var(--main_panel_color)",
                    zIndex: "2",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                },
            },
            ds.contentDiv,
        );
        const spinner = createEl(
            "i",
            {
                classes: ["fas", "fa-circle-notch", "fa-spin"],

                styles: {
                    fontSize: "30px",
                    color: t.text_color,
                },
            },
            div,
        );
        const ellipsis = createEl(
            "div",
            {
                styles: {
                    position: "absolute",
                    overflow: "hide",
                    textAlign: "center",
                    top: "3px",
                    color: t.text_color,
                    textOverflow: "ellipsis",
                    wordBreak: "break-all",
                    fontSize: "16px",
                },
                text: config.title,
            },
            div,
        );
        try {
            // this can go wrong if the dataSource doesn't have data or a dynamic dataLoader.
            // when it goes wrong, it can cause problems outside the creation of this chart
            // - other charts wanting to use similar neededCols end up not having data?
            // links should be initialised here as well...
            await this._getColumnsAsync(dataSource, neededCols);
            await this._addChart(dataSource, config, div, notify);
        } catch (error) {
            this.clearInfoAlerts();
            spinner.remove();
            ellipsis.remove();
            div.style.border = "1px solid var(--border_menu_bar_color)";
            const debugNode = createEl(
                "div",
                {
                    styles: {
                        position: "absolute",
                        width: "100%",
                        height: "100%",
                        display: "flex",
                        justifyContent: "center",
                        alignItems: "center",
                        backdropFilter: "blur(10px)",
                    },
                },
                div,
            );
            const errorObj = error instanceof Error ? 
                error 
                :
                typeof(error) === "string" ?
                    { message: error }
                    :
                    { message: "An error occurred while creating the chart" };
            createMdvPortal(ErrorComponentReactWrapper({ error: errorObj, height, width, extraMetaData: { config } }), debugNode);
            const closeButtonContainer = createEl(
                "div",
                {
                    styles: {
                        position: "absolute",
                        top: "10px",
                        right: "10px",
                    }
                },
                div
            );

            createMenuIcon(
                "fas fa-times",
                {
                    tooltip: {
                        text: "Remove Chart",
                        position: "bottom-right",
                    },
                    func: () => {
                        if (this.charts[config.id] && this.charts[config.id].chart) {
                            const chart = this.charts[config.id].chart;
                            chart.remove();
                            delete this.charts[chart.config.id];
                            this._removeLinks(chart);
                            this._callListeners("chart_removed", chart);
                        }
                        div.remove();
                    },
                },
                closeButtonContainer,
            )
            //not rethrowing doesn't help recovering from missing data in other charts.
            //throw new Error(error); //probably not a great way to handle this
        }
    }

    _getColumnsFromOtherSource(
        dataSource,
        otherDataSource,
        columns,
        indexCol,
        func,
    ) {
        this._getColumnsThen(otherDataSource, columns.concat(indexCol), () => {
            const ds = this.dsIndex[dataSource].dataStore;
            const ods = this.dsIndex[otherDataSource].dataStore;
            const oindex = ods.getColumnIndex(indexCol);
            const ic = ds.columnIndex[indexCol];
            const index = ic.values.map((x) => oindex[x]);
            const colInfo = columns.map((x) => {
                const c1 = ds.columnIndex[x];
                const c2 = ods.columnIndex[x];
                if (c2.values) {
                    c1.values = c2.values;
                }
                if (c2.minMax) {
                    c1.minMax = c2.minMax;
                }
                if (c2.quantiles) {
                    c1.quantiles = c2.quantiles;
                }

                const buf = new SharedArrayBuffer(
                    ds.size * (c1.datatype === "text" ? 1 : 4),
                );
                const arrType =
                    c1.datatype === "text" ? Uint8Array : Float32Array;
                return {
                    col: x,
                    data: buf,
                    arr: new arrType(buf),
                    odata: c2.data,
                };
            });
            for (let n = 0; n < ds.size; n++) {
                const i = index[ic.data[n]];
                for (const c of colInfo) {
                    c.arr[n] = c.odata[i];
                }
            }

            for (const c of colInfo) {
                ds.setColumnData(c.col, c.data);
            }
            func();
        });
    }

    async _getColumnsAsync(dataSource, columns) {
        // there could be links as well as column `fields` in columns.
        // let's make a mechanism for awaiting the link - and initial linked cols.
        const ds = this.dsIndex[dataSource].dataStore;
        const links = columns.map((c) => deserialiseParam(ds, c)).filter(c => typeof c !== 'string');
        const linkPromises = links.map(link => link.initialize());
        await Promise.all(linkPromises);
        const linkedFields = links.flatMap(link => link.fields);
        const allColumns = columns.concat(linkedFields);
        

        // let's add some error handling here...
        const result = await new Promise((resolve) => {
            this._getColumnsThen(dataSource, allColumns, resolve);
        });
        return result;
    }

    _getColumnsThen(dataSource, columns, func) {
        const ds = this.dsIndex[dataSource];
        const dStore = ds.dataStore;
        //check if need to load column data from linked data set
        if (dStore.links) {
            for (const ods in dStore.links) {
                const link = dStore.links[ods];
                if (link.columns) {
                    const otherCols = [];
                    const thisCols = [];
                    for (const c of columns) {
                        if (link.columns.indexOf(c) === -1) {
                            thisCols.push(c);
                        } else if (!dStore.columnIndex[c].data) {
                            otherCols.push(c);
                        }
                    }
                    //get the other datasource's columns first
                    if (otherCols.length > 0) {
                        //get index column
                        this._getColumnsThen(dataSource, [link.index], () => {
                            //then get all the other columns
                            this._getColumnsFromOtherSource(
                                dataSource,
                                ods,
                                otherCols,
                                link.index,
                                () => {
                                    this._getColumnsThen(
                                        dataSource,
                                        thisCols,
                                        func,
                                    );
                                },
                            );
                        });
                        return;
                    }
                }
            }
        }
        const reqCols = columns.filter((x) => {
            //column already loading - but what if something went wrong earlier?
            if (this.columnsLoading[dataSource][x]) {
                return false;
            }
            const col = dStore.columnIndex[x];
            //no record of column- need to load it (plus metadata)
            if (!col) {
                // what if x is something like a MulticolumnQuery?
                if (typeof x !== "string") {
                    // we could make dataStore understand it as a 'field'...
                    // or if we return false to filter it out, chart deserialise can handle it?
                    return false;
                } else {
                    dStore.addColumnFromField(x);
                    return true;
                }
            }
            //only load if has no data
            return !col.data;
        });

        //No columns needed
        //but columns requested by other actions may still be loading
        if (reqCols.length === 0) {
            this._haveColumnsLoaded(columns, dataSource, func);
        }
        //load required columns, then check all requested are loaded
        else {
            this.loadColumnSet(reqCols, dataSource, (failedColumns) => {
                if (failedColumns.length) {
                    console.warn(
                        "got columns with some failures",
                        failedColumns,
                    );
                }
                this._haveColumnsLoaded(columns, dataSource, func);
            });
        }
    }

    /*getIndexedData(dataSource,columns,indexColumn,callback,config={}){
        const col = this.dsIndex[dataSource].dataStore;
        this._getColumnsThen(dataSource,column,indexColumn],()=>{
            const index = ds.getColumnIndex(column);
            const cf = ds.getColorFunction(column,config);
            callback((val)=>{
                cf(index[val])
            })
        })

    }*/

    //check all columns have loaded - if not recursive call after
    //time out, otherwise add the chart
    _haveColumnsLoaded(neededCols, dataSource, func) {
        for (const col of neededCols) {
            if (this.columnsLoading[dataSource][col]) {
                setTimeout(() => {
                    this._haveColumnsLoaded(neededCols, dataSource, func);
                }, 500);
                return;
            }
        }
        func();
    }

    async _addChart(dataSource, config, div, notify = false) {
        //**convert legacy data***********
        const ds = this.dsIndex[dataSource];
        div.innerHTML = "";
        div.style.display = "";
        div.style.alignItems = "";
        div.style.justifyContent = "";
        const chartType = BaseChart.types[config.type];
        // sometimes the constructor doesn't understand active link, but chart.setParams() does...
        // it somewhat appears to work if we first manifest concrete field names, 
        // then set the real params in a setTimeout()...
        // but we're resorting to setTimeout for a few things and it's going to be hard to manage...
        // are we sure that this timeout will happen after the one that instruments decorators etc?
        // seems dodgy to rely on setTimeout order. 
        // *We now have chart.deferredInit() for this purpose*
        // As long as the order in which these methods are called, this should allow us to at least
        // assert that the initialisation is complete before we consider the chart (and ultimately the view) to be loaded.
        // Also setParams() itself *should* always work - but depends on the chart implementation,
        // and not all charts implement that.
        const originalParam = config.param.map(p => deserialiseParam(ds, p));
        config.param = originalParam.flatMap(getConcreteFieldNames);
        const chart = new chartType.class(ds.dataStore, div, config);
        if (originalParam.some(p => typeof p !== 'string')) {
            console.log('setting param with active queries in timeout');
            chart.deferredInit(() => {
                chart.setParams(originalParam);
            });
        }
        this.charts[chart.config.id] = {
            chart: chart,
            dataSource: ds,
        };
        this._makeChartRD(chart, ds);
        // @ts-ignore
        chart.popoutIcon = chart.addMenuIcon(
            "fas fa-external-link-alt",
            "popout",
            {
                func: () => {
                    this._popOutChart(chart);
                },
            },
        );
        chart
            .addMenuIcon("fas fa-times", "remove chart")
            .addEventListener("click", () => {
                chart.remove();
                div.remove();
                delete this.charts[chart.config.id];
                this._removeLinks(chart);
                this._callListeners("chart_removed", chart);
            });

        // @ts-ignore
        if (chart.setupLinks) {
            //phasing out
            if (ds.index_link_to) {
                this._giveChartAccess(
                    chart,
                    this.dsIndex[ds.index_link_to.dataSource].dataStore,
                    ds.index_link_to.index,
                );
            }
            for (const lnk of ds.dataStore.accessOtherDataStore) {
                this._giveChartAccess(
                    chart,
                    this.dsIndex[lnk.dataSource].dataStore,
                    lnk.index,
                );
            }
        }

        //I think this is obsolete now
        //(the above comment is itself very old)
        const cll = ds.column_link_to;
        // @ts-ignore
        if (cll && chart.createColumnLinks) {
            const func = (columns, callback) => {
                //make sure index is loaded before use
                this._getColumnsThen(cll.dataSource, columns, callback);
            };
            // @ts-ignore
            chart.createColumnLinks(
                this.dsIndex[cll.dataSource].dataStore,
                cll.columns,
                func,
            );
        }

        await chart.deferredInitsReady();

        if (notify) {
            this._callListeners("chart_added", chart);
        }
        //check to see if all inital charted loaded , then can call any listeners
        if (this._toLoadCharts) {
            this._toLoadCharts.delete(config);
            if (this._toLoadCharts.size === 0) {
                this._toLoadCharts = undefined;
                if (this.viewData.links) {
                    for (const l of this.viewData.links) {
                        this._setUpLink(l);
                    }
                }
                this._callListeners("view_loaded", this.viewManager.current_view);
            }
        }
        return chart;
    }

    //gives a chart access to another datasource
    _giveChartAccess(chart, dataSource, index) {
        const getDataFunction = async (columns, callback) => {
            //make sure index is loaded before use
            columns.push(index); //pjt: might just push `undefined`?
            await this._getColumnsAsync(dataSource.name, columns);
            callback();
        };
        chart.setupLinks(dataSource, index, getDataFunction);
    }

    //sets up a link between charts
    _setUpLink(link) {
        if (!link.id) {
            link.id = getRandomString();
        }
        switch (link.type) {
            // some other types of link may be interpreted within chart code (e.g. react effect does "view_state")
            case "color_by_column":
                link.set_color = true;
                console.warn(
                    "legacy color_by_column link type - use chart_columnval_link with set_color=true instead",
                );
                addChartLink(link, this);
                break;
            case "chart_columnval_link":
                addChartLink(link, this);
                break;
        }
    }

    //if a chart has been removed, work out which links need to be removed
    _removeLinks(chart) {
        const linksToRemove = [];
        const cid = chart.config.id;
        const links = this.viewData.links;
        if (!links) {
            return;
        }
        for (let i = 0; i < links.length; i++) {
            const link = links[i];
            if (link.source_chart === cid) {
                linksToRemove.push(i);
            }
            try {
                // nb, for example viewState link doesn't have target_charts, it has linked_charts...
                //pjt: slight hack, maybe ok, want to have some more typing on links
                const targets = link.target_charts || link.linked_charts;
                const index = targets.indexOf(cid);
                if (index !== -1) {
                    targets.splice(index, 1);
                    if (targets.length === 0) {
                        linksToRemove.push(i);
                    }
                }
            } catch (error) {
                console.error(
                    `Error removing links for chart '${cid}'`,
                    error,
                    link,
                );
            }
        }
        for (const i of linksToRemove) {
            this.removeLink(i);
        }
    }

    removeLink(linkIndex) {
        const link = this.viewData.links[linkIndex];
        switch (link.type) {
            case "color_by_column":
            case "chart_columnval_link": {
                const chart = this.charts[link.source_chart].chart;
                chart.removeListener(link.id);
            }
        }
        this.viewData.links.splice(linkIndex, 1);
    }

    /**
     * Removes all charts from the ChartManager.
     * If `dataSources` is provided, only charts associated with those data sources will be removed.
     * @param {string[]?} [dataSources] - An array of data source names. If provided, only charts associated with these data sources will be removed.
     */
    removeAllCharts(dataSources) {
        const allCharts = [];
        for (const cn in this.charts) {
            const ch = this.charts[cn];
            if (dataSources && dataSources.indexOf(ch.dataSource.name) === -1) {
                continue;
            }
            allCharts.push([ch.chart, ch.window]);
        }
        for (const ci of allCharts) {
            try {
                if (ci[1]) {
                    ci[1].close();
                }
                ci[0].remove();
                ci[0].div.remove();
            } catch (error) {
                console.error("Error occurred while removing the chart: ", error);
            }
        }
        this.charts = {};
    }

    getAllFilters(dataSorce) {
        const charts = this.getAllCharts(dataSorce);
        const fs = [];
        for (const c of charts) {
            const filter = c.getFilter();
            if (filter) {
                fs.push(filter);
            }
        }
        return fs;
    }

    getChart(id) {
        const cinfo = this.charts[id];
        if (!cinfo) {
            return null;
        }
        return cinfo.chart;
    }

    /**
     * Get all the charts for a data sorce
     * @param {string} dataSource - The name of the data source
     * @returns {Array} - An array of chart objects
     */
    getAllCharts(dataSource) {
        const charts = [];
        for (const id in this.charts) {
            const ch = this.charts[id];
            if (ch.dataSource.name === dataSource) {
                charts.push(ch.chart);
            }
        }
        return charts;
    }

    setChartsAsGrid(rowLength = 5, size = [300, 300], margin = 10) {
        let top = margin;
        let left = margin;
        let rowSize = 0;
        for (const id in this.charts) {
            const info = this.charts[id];
            const d = info.chart.getDiv();
            d.style.left = `${left}px`;
            d.style.top = `${top}px`;
            //info.chart.setSize(size[0],size[1]);
            left += size[0] + margin;
            rowSize++;
            if (rowSize === rowLength) {
                rowSize = 0;
                left = margin;
                top += size[1] + margin;
            }
        }
    }

    addButton(text, callback, tooltip) {
        createEl(
            "button",
            {
                classes: ["ciview-button"],
                text: text,
                styles: {
                    position: "fixed",
                    bottom: "40px",
                    right: "40px",
                    fontSize: "18px",
                    zIndex: "100",
                },
            },
            this.contentDiv,
        ).addEventListener("click", () => callback());
    }

    _popOutChart(chart) {
        popoutChart(chart);
    }

    _sendAllChartsToBack(ds) {
        for (const id in this.charts) {
            const c = this.charts[id];
            if (ds === c.dataSource) {
                c.chart.div.style.zIndex = "";
            }
        }
    }

    /**
     * @param {BaseChart} chart 
     * @param {import("@/charts/charts").DataSource} ds
     */
    _makeChartRD(chart, ds) {
        //if (!ds) console.error(`_makeChartRD called without ds - resize / drag etc may not work properly`);
        //^^ actually doesn't make much difference to non-gridStack in practice.
        if (
            ds &&
            this.gridStack &&
            this.viewData.dataSources[ds.name].layout === "gridstack"
        ) {
            this.gridStack.manageChart(chart, ds, this._inInit);
            return;
        }
        const div = chart.getDiv();
        makeDraggable(div, {
            handle: ".ciview-chart-title",
            contain: "topleft",
            ondragstart: (e) => {
                this._sendAllChartsToBack(ds);
                div.style.zIndex = 2;
            },
        });
        makeResizable(div, {
            onResizeStart: () => {
                this._sendAllChartsToBack(ds);
                div.style.zIndex = 2;
            },
            onresizeend: (width, height) => chart.setSize(width, height),
        });
    }
}

export default ChartManager;
