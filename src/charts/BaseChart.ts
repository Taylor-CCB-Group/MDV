import { getRandomString } from "../utilities/Utilities";
import { ContextMenu } from "../utilities/ContextMenu.js";
import { createEl } from "../utilities/Elements.js";
import { type ChartTypeMap, chartTypes } from "./ChartTypes";
import DebugJsonDialogReactWrapper from "../react/components/DebugJsonDialogReactWrapper";
import SettingsDialogReactWrapper from "../react/components/SettingsDialogReactWrapper";
import { makeAutoObservable, action, autorun, IReactionDisposer, IAutorunOptions } from "mobx";
import type DataStore from "@/datastore/DataStore";
import type { BaseDialog } from "@/utilities/Dialog";
import type { DataColumn, FieldName, GuiSpec, GuiSpecs } from "./charts";
import type Dimension from "@/datastore/Dimension";
import { g } from "@/lib/utils";
import { serialiseConfig, initialiseConfig } from "./chartConfigUtils";
type ChartEventType = string;
type Listener = (type: ChartEventType, data: any) => void;
type LegacyColorBy = { column: DataColumn<any> }
export type BaseConfig = {
    id: string;
    size: [x: number, y: number];
    title: string;
    legend: string;
    type: string;
    param: FieldName[]; // | string,
    //nb we might want to represent columns as something other than string soon, in which case this will need to be updated
    color_by?: FieldName | LegacyColorBy;
    title_color?: string;
};
type ColumnChangeEvent = { columns: FieldName[], hasFiltered: boolean };
type ColorOptions = any;
class BaseChart<T> {
    config: any;
    __doc__: Document;
    dataStore: DataStore;
    listeners: Record<string, Listener>;
    title: HTMLDivElement;
    div: HTMLElement;
    titleBar: HTMLDivElement;
    menuSpace: HTMLDivElement;
    contentDiv: HTMLDivElement;
    resetButton: HTMLButtonElement | HTMLSpanElement;
    contextMenu: ContextMenu;
    dialogs: BaseDialog[] = [];
    legendIcon: HTMLElement;
    observable: { container: HTMLElement };
    width = 0;
    height = 0;
    legend: any;
    activeQueries: any;
    /**
     * The base constructor
     * @param {import("./charts.js").DataStore} dataStore - The datastore object that contains the data for this chart
     * @param {string | HTMLDivElement} div - The id of the div element or the element itself to house the chart
     * @param {Object} config - The config describing the chart
     */
    constructor(dataStore: DataStore, div: HTMLDivElement | string, config: T & BaseConfig) {
        //******adapt legacy configs
        if (config.color_by) {
            const { color_by } = config;
            //nb we might want to represent columns as something other than string soon, in which case this will need to be updated
            if (typeof color_by !== "string" && color_by.column) {
                config.color_by = color_by.column.field;
            }
        }
        //**********

        //copy the config, 
        // this.config = JSON.parse(JSON.stringify(config));
        //^^ previously, only react charts had observable config
        //and make it observable etc... may apply some other processing e.g. in case of 'virtual columns'
        //   we've been experimenting with applying this to all charts
        //   but there are issues with mobx actions that result in mutations to the config
        //... so perhaps the idea of keeping that property to react charts is a good one?
        this.config = initialiseConfig(config, this);

        //required in case added to separate browser window
        this.__doc__ = document;

        this.dataStore = dataStore;
        this.listeners = {};

        //create the DOM elements
        const divElement = typeof div === "string" ? document.getElementById(div) : div;
        if (!divElement) throw new Error("failed to get div element");
        this.div = divElement;
        this.div.classList.add("ciview-chart-panel");
        this.titleBar = createEl("div", {
            classes: ["ciview-chart-titlebar"],
        });

        this.title = createEl(
            "div",
            {
                classes: ["ciview-chart-title"],
                text: config.title,
            },
            this.titleBar,
        );

        this.menuSpace = createEl(
            "div",
            {
                classes: ["ciview-chart-menuspace"],
            },
            this.titleBar,
        );

        this.contentDiv = createEl("div", {
            classes: ["ciview-chart-content"],
        });

        if (config.title_color) {
            this.titleBar.style.backgroundColor = config.title_color;
        }

        //reset button
        this.resetButton = this.addMenuIcon("fas fa-sync", "remove filter");
        this.resetButton.style.display = "none";
        this.resetButton.addEventListener("click", () => {
            this.removeFilter();
            this.resetButton.style.display = "none";
        });

        //register with datastore to listen to filter events
        this.dataStore.addListener(this.config.id, (type, data) => {
            if (type === "filtered") {
                this.onDataFiltered(data);
            } else if (type === "data_changed") {
                this.onDataChanged(data);
            } else if (type === "data_added") {
                this.onDataAdded(data);
            } else if (type === "data_highlighted") {
                if (data.source === this) {
                    this._callListeners("data_highlighted", data);
                }
                if (this.onDataHighlighted) {
                    this.onDataHighlighted(data);
                }
            }
        });

        //set up context menu and icon which opens it
        this.contextMenu = new ContextMenu((data) => {
            const menu = this.getContextMenu(data);
            menu.push({
                text: "debug chart",
                icon: "fas fa-bug",
                func: () => {
                    window.mdv.debugChart = this;
                    this.dialogs.push(
                        new DebugJsonDialogReactWrapper(this.config, this),
                    );
                },
            });
            menu.push({
                text: "copy config JSON to clipboard",
                icon: "fas fa-copy",
                func: () =>
                    navigator.clipboard.writeText(
                        JSON.stringify(this.config, null, 2),
                    ),
            });
            return menu;
        });

        this.addMenuIcon("fas fa-bars", "more").addEventListener("click", (e) =>
            this.contextMenu.show(e),
        );

        this.addMenuIcon("fas fa-cog", "settings").addEventListener(
            "click",
            (e) => this._openSettingsDialog(e),
        );

        this.addMenuIcon("fas fa-copy", "duplicate chart").addEventListener(
            "click",
            () => {
                const newConfig = JSON.parse(JSON.stringify(this.getConfig()));
                newConfig.id = getRandomString();
                newConfig.title = `${newConfig.title} copy`;
                window.mdv.chartManager.addChart(
                    this.dataStore.name,
                    newConfig,
                );
            },
        );

        //info icon
        this.legendIcon = this.addMenuIcon(
            "fas fa-info",
            config.legend || "No description",
            { size: "medium" },
        );

        let oldSize = config.size;
        this.contentDiv.addEventListener(
            "fullscreenchange",
            action(() => {
                //nb, debounced version of setSize also being called by gridstack - doesn't seem to cause any problems
                if (document.fullscreenElement) {
                    if (this.contentDiv !== document.fullscreenElement)
                        console.error("unexpected fullscreen element");
                    this.observable.container = this.contentDiv;
                    const rect = window.screen;
                    this.setSize(rect.width, rect.height);
                    for (const d of this.dialogs) {
                        d.setParent(this.contentDiv);
                    }
                } else {
                    this.observable.container = this.__doc__.body;
                    this.setSize(...oldSize);
                    for (const d of this.dialogs) {
                        d.setParent(null);
                    }
                }
            }),
        );
        this.addMenuIcon("fas fa-expand", "fullscreen", {
            func: async () => {
                oldSize = this.config.size;
                await this.contentDiv.requestFullscreen();
            },
        });

        this.div.append(this.titleBar);
        this.div.append(this.contentDiv);
        this.listeners = {};

        //work out width and height based on container
        this._setDimensions();
        // provide some observable mobx state for useOuterContainer()
        this.observable = makeAutoObservable({
            container: this.__doc__.body,
        });
    }
    reactionDisposers: IReactionDisposer[] = [];
    mobxAutorun(fn: ()=>void, opts?: IAutorunOptions) {
        const disposer = autorun(fn, opts);
        this.reactionDisposers.push(disposer);
        return disposer;
    }
    _getContentDimensions() {
        return {
            //PJT to review re. gridstack.
            top: 5,
            left: 5,
            height: this.height,
            width: this.width - 5,
        };
    }
    get dataSource() {
        //this fails if we're still in chart construction
        //return window.mdv.chartManager.charts[this.config.id].dataSource;
        const name = this.dataStore.name;
        return window.mdv.chartManager.dataSources.find((ds) => ds.name === name);
    }

    getFilter(): any {}

    /** appears to be an appropriate type signature from glancing at uses in the code */
    setFilter(active: boolean) {}

    /**
     * Adds a listener to the datastore that will be called when an event occurs,
     * passing the event type and any data. There are the following different types
     * of event:-
     * <ul>
     * <li> removed - called when the chart has been removed </li>
     * </ul>
     * @param id - a unique id indetifying the listener
     * @param listener - a function that accepts two paramaters, the type
     * of event and the data associated with it
     */
    addListener(id: string, listener: Listener) {
        this.listeners[id] = listener;
    }

    /**
     * Removes the specified listener from the chart
     * @param id The id of the listener to remove
     */
    removeListener(id: string) {
        if (!Object.hasOwn(this.listeners, id)) console.warn("listener not found");
        delete this.listeners[id];
    }

    _callListeners(type: ChartEventType, data: any) {
        for (const id in this.listeners) {
            this.listeners[id](type, data);
        }
    }

    /**
     * Adds a menu icon with tooltip to the title bar
     * @param icon - the css classs of the icon (space delimited).
     * @param tooltip - the tooltip text
     * @param config - extra inormation about the icon/tooltip
     * @param [config.size=small] - the size of the tooltip
     * @param [config.position=bottom] - the position of the tooltip.
     * @param [config.func] - a function that is called when the icon is clicked
     * @returns - the icon
     */
    addMenuIcon(icon: string, tooltip: string, config:
        {size?: string, position?: string, func?: (e: MouseEvent)=>void} = {}
    ) {
        const sp = createEl(
            "span",
            {
                "aria-label": tooltip,
                "data-microtip-color": "red",
                role: "tooltip",
                "data-microtip-size": config.size || "small",
                "data-microtip-position": config.position || "bottom-left",
                styles: {
                    margin: "0px 1px",
                },
            },
            this.menuSpace,
        );

        createEl(
            "i",
            {
                classes: ["ciview-chart-icon"].concat(icon.split(" ")), //a11y - we could use an actual button
            },
            sp,
        );
        if (config.func) {
            sp.addEventListener("click", (e) => config.func?.(e));
        }
        return sp;
    }

    /**
     * needs to be implemented by subclasses
     */
    removeFilter() {}

    /**
     * Called by the datastore when the data is filtered. Needs to
     * be implemented on any subclasses.
     * @param dim - the dimension that has been filtered
     */
    onDataFiltered(dim?: Dimension) {}

    colorByColumn?(c: FieldName): void;
    colorByDefault?(): void;
    /**Check if chart is composed of any columns whose data has
     * changed. if so re-calculate and re-draw chart (call onDataFiltered)
     * @param data - a list of column/fields whose data has been modified
     * @param data.columns a list of column ids whose data has changed
     * @param data.hasFiltered Whether a 'filtered' callback has already been
     * issued
     */
    onDataChanged(data: ColumnChangeEvent) {
        const columns = data.columns;
        //update any charts which use data from the columns
        //(if they haven't already been updated by the filter changing)
        if (!data.hasFiltered) {
            let cols = this.config.param;
            let isDirty = false;
            if (typeof this.config.param === "string") {
                cols = [this.config.param];
            }
            for (const p of cols) {
                if (columns.indexOf(p) !== -1) {
                    isDirty = true;
                    break;
                }
            }
            if (isDirty) {
                this.onDataFiltered();
            }
        }
        //recolor any charts coloured by the column
        if (columns.indexOf(this.config.color_by) !== -1) {
            this.colorByColumn?.(this.config.color_by);
        }
    }

    getColorLegend() {
        const conf = {
            overideValues: {
                colorLogScale: this.config.log_color_scale,
            },
        };
        this._addTrimmedColor(this.config.color_by, conf);

        return this.dataStore.getColorLegend(this.config.color_by, conf);
    }

    // getQunatile;

    _addTrimmedColor(column: FieldName, conf: any) {
        const tr = this.config.trim_color_scale;
        const col = this.dataStore.columnIndex[column];
        if (tr && tr !== "none") {
            if (col.quantiles && col.quantiles !== "NA") {
                conf.overideValues.min = col.quantiles[tr][0];
                conf.overideValues.max = col.quantiles[tr][1];
            }
        }
    }

    /**
     * adds (or removes) the color legend depending on the chart's
     * config color_legend.display value - assumes chart has a
     * get colorLegend method
     */
    setColorLegend() {
        if (!this.config.color_legend.display) {
            if (this.legend) {
                this.config.color_legend.pos = [
                    this.legend.offsetLeft,
                    this.legend.offsetTop,
                ];
                this.legend.remove();
                this.legend = undefined;
            }
            return;
        }
        const box = this._getContentDimensions();
        let lt: string | number = 0;
        let ll: string | number = 0;
        if (this.legend) {
            ll = this.legend.style.left;
            lt = this.legend.style.top;
            this.legend.remove();
        } else {
            const cl = this.config.color_legend;
            if (!cl.pos) {
                cl.pos = [box.left, box.top];
            }
            ll = `${cl.pos[0]}px`;
            lt = `${cl.pos[1]}px`;
        }
        this.legend = this.getColorLegend();
        if (!this.legend) {
            console.warn("no color legend");
            return;
        }
        this.contentDiv.append(this.legend);

        this.legend.style.left = ll;
        this.legend.style.top = lt;
        this.legend.__doc__ = this.__doc__;
    }

    getColorFunction(column: FieldName, asArray?: boolean) {
        this.config.color_by = column;
        const conf = {
            asArray: asArray,
            overideValues: {
                colorLogScale: this.config.log_color_scale,
                fallbackOnZero: this.config.fallbackOnZero,
            },
        };
        this._addTrimmedColor(column, conf);

        const colorFunc = this.dataStore.getColorFunction(column, conf);

        if (!this.config.color_legend) {
            this.config.color_legend = {
                display: true,
            };
        }

        this.setColorLegend();
        return colorFunc;
    }

    /**Checks to see if the column is used in the chart
     * If so, the chart will be removed but no callbacks will be involved
     * @returns `true` if the chart has been removed
     */
    onColumnRemoved(column: FieldName) {
        let cols = this.config.param;
        let isDirty = false;
        if (typeof this.config.param === "string") {
            cols = [this.config.param];
        }
        for (const p of cols) {
            if (column === p) {
                isDirty = true;
                break;
            }
        }
        if (isDirty) {
            this.remove(false);
            return true;
        }
        if (this.colorByColumn) {
            if (this.config.color_by === column) {
                this.config.color_by = undefined;
                this.colorByDefault?.();
            }
        }
        return false;
    }

    onDataAdded(newSize: number) {
        this.onDataFiltered();
    }

    onDataHighlighted(data: any) {}
    _tooltip?: HTMLDivElement;
    addToolTip() {
        this._tooltip = createEl(
            "div",
            {
                classes: ["ciview-tooltip"],
                stlyles: {
                    display: "none",
                    position: "absolute",
                },
            },
            this.__doc__.body,
        );
    }

    showToolTip(e: MouseEvent, msg: string) {
        if (!this._tooltip) this.addToolTip();
        if (!this._tooltip) throw 'expected tooltip to be added by now';
        this._tooltip.innerHTML = msg;
        this._tooltip.style.display = "inline-block";
        this._tooltip.style.left = `${3 + e.clientX}px`;
        this._tooltip.style.top = `${3 + e.clientY}px`;
    }

    hideToolTip() {
        if (!this._tooltip) return;
        this._tooltip.style.display = "none";
    }

    /**
     * Just removes the DOM elements, subclasses should do their own cleanup
     */
    remove(notify?: boolean) {
        this.titleBar.remove();
        this.contentDiv.remove();
        this.dataStore.removeListener(this.config.id);
        if (this._tooltip) {
            this._tooltip.remove();
        }
        for (const d of this.dialogs) {
            d.close();
        }
        for (const disposer of this.reactionDisposers) {
            disposer();
        }
        // dynamic props?
    }
    removeLayout?(): void;
    /**
     * Returns information about which controls to add to the settings dialog.
     * Subclasses should call this method and then add their own controls e.g.
     * <pre style='background:lightgray;padding:10px'>
     * getSettings(){
     *     let settings = super.getSettings();
     *     return settings.concat([{
     *       label:"A value"
     *       type:"slider",
     *       default_value:this.config.value,
     *       max:10,
     *       min:10,
     *       func:x=>{
     *           this.config.value=x;
     *           //update chart
     *       }
     *     }]);
     * }
     * </pre>
     *
     * wrapping controls with a call to @{link g} will perform type checking
     * todo- specifiy the link to `g` above properly/ document better
     * @returns an array of objects describing tweakable parameters
     */
    getSettings() {
        const c = this.config;
        const settings: GuiSpecs = [
            {
                type: "text",
                label: "Chart Name",
                current_value: c.title,
                func: (v) => {
                    this.setTitle(v);
                },
            },
        ];
        const colorOptions = this.getColorOptions();

        if (colorOptions.colorby) {
            //cannot color by unique (at the moment)
            const filter = (
                colorOptions.colorby === "all"
                    ? [
                          "int32",
                          "text",
                          "integer",
                          "double",
                          "text16",
                          "multitext",
                      ]
                    : colorOptions.colorby
            );

            const colorSettings: GuiSpecs = [];
            settings.push(g({
                type: "folder",
                label: "Color Settings",
                current_value: colorSettings,
            }));

            // change this to use setting with type: "column" - so it should understand "filter" in a similar way
            colorSettings.push(g({
                label: "Color By",
                type: "column",
                current_value: c.color_by,
                columnSelection: { filter },
                func: (x) => {
                    if (x === "_none") {
                        c.color_by = undefined;
                        this.colorByDefault?.();
                    } else {
                        c.color_by = x;
                        this.colorByColumn?.(x);
                    }
                },
            }));
            if (colorOptions.color_overlay !== undefined) {
                colorSettings.push({
                    label: "Color Overlay",
                    type: "slider",
                    current_value: c.color_overlay,
                    func: (x) => {
                        c.color_overlay = x;
                        this.colorByColumn?.(c.color_by);
                    },
                });
            }
            colorSettings.push({
                label: "Show Color Legend",
                type: "check",

                current_value: c.color_legend ? c.color_legend.display : true,
                func: (x) => {
                    if (!c.color_by) {
                        return;
                    }
                    c.color_legend.display = x;
                    this.setColorLegend();
                },
            });
            colorSettings.push({
                label: "SymLog Color Scale",
                type: "check",

                current_value: c.log_color_scale,
                func: (x) => {
                    c.log_color_scale = x;
                    if (c.color_by) {
                        this.colorByColumn?.(c.color_by);
                    }
                },
            });
            colorSettings.push({
                label: "Treat zero as missing",
                type: "check",

                current_value: c.fallbackOnZero,
                func: (x) => {
                    c.fallbackOnZero = x;
                    if (c.color_by) {
                        this.colorByColumn?.(c.color_by);
                    }
                },
            });
            colorSettings.push({
                type: "radiobuttons",
                label: "Trim Color Scale to Percentile",
                current_value: c.trim_color_scale || "none",
                choices: [
                    ["No Trim", "none"],
                    ["0.001", "0.001"],
                    ["0.01", "0.01"],
                    ["0.05", "0.05"],
                ],
                func: (v) => {
                    c.trim_color_scale = v;
                    if (c.color_by) {
                        this.colorByColumn?.(c.color_by);
                    }
                },
            });
        }

        return settings;
    }

    unpinIcon?: HTMLElement;
    _addUnpinIcon() {
        this.unpinIcon = this.addMenuIcon("fas fa-thumbtack", "unpin chart");
        this.unpinIcon.addEventListener("click", (e) => {
            this.unpinIcon?.remove();
            this.unpinChart?.();
        });
    }
    settingsDialog?: BaseDialog;
    _openSettingsDialog(e: MouseEvent) {
        if (!this.settingsDialog) {
            // the dialog will set `this.settingsDialog = null` (and remove itself from this.dialogs) when it closes.
            this.settingsDialog = new SettingsDialogReactWrapper(this, [
                e.pageX,
                e.pageY,
            ]);
            this.dialogs.push(this.settingsDialog);
        }
    }
    isPinned = false;
    pinChart?(): void;
    unpinChart?(): void;
    getChartData?(): any;
    /**
     * @returns an array of `ContextMenuItems`
    */
    addToContextMenu?(): any;
    getContextMenu(data?: any) {
        type ContextMenuItem = { text: string, icon: string, func: () => void };
        let menu: ContextMenuItem[] = [];
        if (this.getImage) {
            menu = menu.concat([
                {
                    text: "create svg image",
                    icon: "far fa-image",
                    func: () => this.downloadImage("svg"),
                },
                {
                    text: "create png image",
                    icon: "far fa-image",
                    func: () => this.downloadImage("png"),
                },
            ]);
        }
        if (this.pinChart && !this.isPinned) {
            menu.push({
                text: "pin chart",
                icon: "fas fa-thumbtack",
                func: () => {
                    this.pinChart?.();
                    this._addUnpinIcon();
                },
            });
        }

        if (this.getChartData) {
            menu.push({
                text: "download data",
                icon: "fas fa-download",
                func: () => this.downloadData(),
            });
        }
        if (this.addToContextMenu) {
            menu = menu.concat(this.addToContextMenu());
        }
        return menu;
    }

    /**
     * Returns the dom element that the chart is attached to
     */
    getDiv() {
        return this.div;
    }
    extra_legends?: any[];
    /**
     * Instructs the chart to use a different document. This is only required if you are
     * going to add the chart to a different browser window
     * @param doc - the document that the chart will use
     */
    changeBaseDocument(doc: Document) {
        //this needs to be reviewed for popout windows
        // - mouse events need to be on the right window ✅
        // - dialogs need to be on the right window ✅ / transferred
        // - fullscreen... some permutations of dialog behaviour / mouse events on deck are a bit odd
        //   ^^ that's not a changeBaseDocument() thing, it's a fullscreen thing...
        this.contextMenu.__doc__ = doc;
        action(() => (this.observable.container = doc.body))();
        this.__doc__ = doc;
        for (const d of this.dialogs) {
            // d.close();
            // how about changing the parent of the dialog?
            // like we do with fullscreen (make sure drag works after this)
            d.setParent(doc.body);
        }
        if (this.legend) {
            this.legend.__doc__ = doc;
        }
        if (this._tooltip) {
            this._tooltip.remove();
            this.addToolTip();
        }
        if (this.extra_legends) {
            for (const l of this.extra_legends) {
                //@ts-expect-error
                if (this[l]) {
                    //@ts-expect-error
                    this[l].__doc__ = doc;
                }
            }
        }
    }

    _setConfigValue(conf: T, value: keyof T, def: any) {
        if (conf[value] === undefined) {
            conf[value] = def;
        }
        return conf[value];
    }

    drawChart() {}

    /**
     * Returns a copy of the chart's config
     */
    getConfig() {
        if (this.legend) {
            this.config.color_legend = {
                display: true,
                pos: [this.legend.offsetLeft, this.legend.offsetTop],
            };
        }
        return serialiseConfig(this);
    }

    getColorOptions(): ColorOptions {
        return {};
    }

    /**
     * Sets the size of the graph. If no parameters are supplied
     * then the graph will be resized based on it container.
     * Subclasses should overide this, but call the super method
     * first which will calculate width and height of the content div
     * @param x - The new width
     * @param y The new height;
     */
    setSize(x?: number, y?: number) {
        //if supplied change the div dimensions
        if (x) {
            this.div.style.height = `${y}px`;
            this.div.style.width = `${x}px`;
        }
        //calculate width and height based on outer div
        this._setDimensions();
    }

    _setDimensions() {
        const rect = this.div.getBoundingClientRect();
        const y = Math.round(rect.height);
        const x = Math.round(rect.width);
        this.config.size = [x, y];
        const rect2 = this.contentDiv.getBoundingClientRect();
        this.height = Math.round(rect2.height);
        this.width = x;
    }

    getImage?(callback: (resp: any) => void, im_type: "png" | "svg"): void;
    /**
     * Downloads an image of the chart
     * @param im_type - either svg or png
     */
    downloadImage(im_type: "svg" | "png") {
        const originalColor = this.contentDiv.style.color;
        this.contentDiv.style.color = "black";
        this.getImage?.((resp) => {
            const link = document.createElement("a");
            const name = this.config.title || "image";
            link.download = `${name}.${im_type}`;
            if (im_type === "svg") {
                link.href = `data:image/svg+xml,${encodeURIComponent(resp)}`;
            } else {
                let url = resp.toDataURL("image/png");
                url = url.replace(
                    /^data:image\/png/,
                    "data:application/octet-stream",
                );
                link.href = url;
            }
            link.click();
            link.remove();
            this.contentDiv.style.color = originalColor;
        }, im_type);
    }

    downloadData() {
        const blob = this.getChartData?.();
        if (!blob) throw new Error("No data to download");
        const save = createEl(
            "a",
            {
                download: this.config.title,
                target: "_blank",
                href: window.URL.createObjectURL(blob),
            },
            this.__doc__.body,
        );
        save.click();
        save.remove();
    }

    setTitle(title: string) {
        this.title.textContent = title;
        // avoid issue with mobx autorun mutating this when reacting to config.title change
        if (this.config.title !== title) {
            this.config.title = title;
        }
    }

    getImageFromSVG(svg: any|SVGElement, callback: (canvas: HTMLCanvasElement, ctx: CanvasRenderingContext2D) => void) {
        const copy = svg.cloneNode(true);
        copyStylesInline(copy, svg);
        const canvas = document.createElement("canvas");
        //const bbox = svg.getBBox();
        copy.style.top = "0px";
        canvas.width = svg.width.baseVal.value;
        canvas.height = svg.height.baseVal.value;
        const ctx = canvas.getContext("2d");
        if (!ctx) throw "Could not get 2d context from canvas";
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        const data = new XMLSerializer().serializeToString(copy);
        const DOMURL = window.URL || window.webkitURL || window;
        const img = new Image();
        const svgBlob = new Blob([data], {
            type: "image/svg+xml;charset=utf-8",
        });
        const url = DOMURL.createObjectURL(svgBlob);
        img.src = url;
        img.onload = () => {
            ctx.drawImage(img, 0, 0);
            callback(canvas, ctx);
        };
    }
    static types: ChartTypeMap;
}

/**
 * A dictionary of all the chart types
 */
BaseChart.types = chartTypes;

function copyStylesInline(
    destinationNode: HTMLElement,
    sourceNode: HTMLElement
): void {
    const containerElements = ["svg", "g"];
    for (let cd = 0; cd < destinationNode.childNodes.length; cd++) {
        const child = destinationNode.childNodes[cd] as HTMLElement;
        const sourceChild = sourceNode.childNodes[cd] as HTMLElement;

        if (!child || !sourceChild) continue;

        if (containerElements.includes(child.tagName.toLowerCase())) {
            copyStylesInline(child, sourceChild);
            continue;
        }

        const style =
            (sourceChild as any).currentStyle ||
            window.getComputedStyle(sourceChild);

        if (!style) continue;

        for (let st = 0; st < style.length; st++) {
            const property = style[st];
            child.style.setProperty(
                property,
                style.getPropertyValue(property)
            );
        }
    }
}


export default BaseChart;
