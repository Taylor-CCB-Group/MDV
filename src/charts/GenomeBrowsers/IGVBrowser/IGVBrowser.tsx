import BaseChart, { type BaseConfig } from "../../BaseChart";
import { BaseReactChart } from "../../../react/components/BaseReactChart";
import type DataStore from "../../../datastore/DataStore";
import IGVBrowserComponent from "./IGVBrowserComponent";
import { g } from "../../../lib/utils";
import { action,toJS } from "mobx";
import {
    type IGVBrowserLocation,
    type IGVBrowserViewMargins,
    type IGVBaseFeature,
    applyViewMargins,
    buildBaseFeatures,
    getLocationFieldsFromGenome,
    locationFromFieldValues,
} from "./igvUtils";

interface IGVBrowserConfig extends BaseConfig {
    view_margins?: IGVBrowserViewMargins;
    location?: IGVBrowserLocation;
    feature_label?: string;
    max_initial_features?: number;
    igv_config?: Record<string, unknown>;
    highlight_selected_region?: boolean;
    base_track_visibility_window?: number;
}

const BASE_TRACK_ID = "_mdv_base_track";
const DEFAULT_MAX_INITIAL_FEATURES = 500000;
const DEFAULT_BASE_TRACK_VISIBILITY_WINDOW = 300000000;

function parseCoordinate(value: unknown): number | null {
    if (typeof value === "number") {
        return Number.isFinite(value) ? value : null;
    }
    if (typeof value === "string") {
        const normalized = value.replaceAll(",", "").trim();
        if (!normalized) return null;
        const parsed = Number(normalized);
        return Number.isFinite(parsed) ? parsed : null;
    }
    return null;
}

function toSafeLocusString(location: IGVBrowserLocation | undefined): string | null {
    const locations = Array.isArray(location) ? location : location ? [location] : [];
    const loci: string[] = [];

    for (const loc of locations) {
        if (!loc || typeof loc.chr !== "string" || !loc.chr.trim()) {
            continue;
        }
        const start = parseCoordinate(loc.start);
        const end = parseCoordinate(loc.end);
        if (start === null || end === null) {
            continue;
        }

        const s = Math.max(1, Math.floor(Math.min(start, end)));
        const e = Math.max(s, Math.floor(Math.max(start, end)));
        if (!Number.isFinite(s) || !Number.isFinite(e) || e < s) {
            continue;
        }
        loci.push(`${loc.chr}:${s}-${e}`);
    }

    return loci.length ? loci.join(" ") : null;
}

function getPathWithoutQueryOrHash(url: string): string {
    return url.split(/[?#]/, 1)[0] ?? url;
}

function appendSuffixBeforeQueryOrHash(url: string, suffix: string): string {
    const match = /^(?<path>[^?#]*)(?<query>\?[^#]*)?(?<hash>#.*)?$/.exec(url);
    if (!match || !match.groups) {
        return `${url}${suffix}`;
    }

    const path = match.groups.path;
    if (path.toLowerCase().endsWith(suffix.toLowerCase())) {
        return url;
    }

    const query = match.groups.query ?? "";
    const hash = match.groups.hash ?? "";
    return `${path}${suffix}${query}${hash}`;
}




import { observable, runInAction } from "mobx";

class IGVBrowser extends BaseReactChart<IGVBrowserConfig> {
    browser: any;
    baseTrack: any;

    igvContainer:{div: HTMLDivElement, css: string} | null = null;
    baseFeatures: IGVBaseFeature[] = [];
    colorFunction: ((rowIndex: number) => unknown) | null = null;
    locationFields: string[] = [];
    isSvs = false;
    guardrailMessage = "";
  
    // Stable observable spec objects for Track Name/URL inputs and status; set in getSettings()
    private _pendingNameSpec: { current_value: string } | null = null;
    private _pendingUrlSpec: { current_value: string } | null = null;
    private _trackStatusSpec: { current_value: string } | null = null;
    private _setSearchPending: ((pending: boolean) => void) | null = null;

    constructor(dataStore: DataStore, div: string | HTMLDivElement, config: IGVBrowserConfig) {
        config.view_margins = config.view_margins || { type: "percentage", value: 20 };
        config.max_initial_features =
            config.max_initial_features && config.max_initial_features > 0
                ? config.max_initial_features
                : DEFAULT_MAX_INITIAL_FEATURES;
        config.highlight_selected_region = config.highlight_selected_region ?? false;
        config.base_track_visibility_window =
            typeof config.base_track_visibility_window === "number" &&
            config.base_track_visibility_window > 0
                ? config.base_track_visibility_window
                : DEFAULT_BASE_TRACK_VISIBILITY_WINDOW;
        super(dataStore, div, config, IGVBrowserComponent);
        console.log("IGVBrowser config:", toJS(this.config));
        this.locationFields = getLocationFieldsFromGenome(this.dataStore.genome) || [];
        this.isSvs = Boolean(this.dataStore.genome?.svs?.sv_columns);
       
      
    }

    keepShadowDom(igvContainer: HTMLDivElement) {
         //keep the style sheet for popout windows as it needs to be added
        //to the shadow dom //keep the style sheet for popout windows as it needs to be added
        const sr = igvContainer.shadowRoot?.adoptedStyleSheets[0];
        let css = "";
        if (sr) {
             for (let rule of sr.cssRules) {
                css += rule.cssText + ' ';
            }
        }
        this.igvContainer = { div: igvContainer, css };
    }

     //keep the style sheet for popout windows as it needs to be added
        //to the shadow dom
    changeBaseDocument(doc: Document) {
        super.changeBaseDocument(doc);
        //create a style element in the right context
        if (!this.igvContainer) return;
        const styleElement = doc.createElement('style');
        styleElement.textContent = this.igvContainer.css;
        //replace it in the shadow root
        const  shadowRoot = this.igvContainer.div.shadowRoot;
        if (shadowRoot) {
            shadowRoot.adoptedStyleSheets = [];
            shadowRoot.appendChild(styleElement);
        }
    }

    clearPendingTrackFields() {
        runInAction(() => {
            if (this._pendingNameSpec) this._pendingNameSpec.current_value = "";
            if (this._pendingUrlSpec) this._pendingUrlSpec.current_value = "";
            if (this._trackStatusSpec) this._trackStatusSpec.current_value = "";
        });
    }

    setSearchPendingHandler(handler: ((pending: boolean) => void) | null) {
        this._setSearchPending = handler;
    }

    private setSearchPending(pending: boolean) {
        this._setSearchPending?.(pending);
    }
    //gets all the features in format igv understands chr,start and end
    //Also id which is the datastore index
    getAllFeatures(){
        const rows: Record<string, unknown>[] = new Array(this.dataStore.size);
        const columns = [...this.locationFields];
   
        for (let i = 0; i < this.dataStore.size; i++) {
            rows[i] = this.dataStore.getRowAsObject(i, columns) as Record<string, unknown>;
        }
        return buildBaseFeatures(rows, this.locationFields, this.isSvs, this.config.max_initial_features || DEFAULT_MAX_INITIAL_FEATURES).features;
    }


    getInitialBrowserConfig() {
        if (this.baseFeatures.length === 0) {
            this.baseFeatures = this.getAllFeatures();
        }
        const initialConfig = {
            genome: this.dataStore.genome?.assembly  || "hg38",
            ...(this.config.igv_config || {}),
            tracks: [
                {
                    id: BASE_TRACK_ID,
                    name: "MDV Features",
                    type: "mdv_feature_track",
                    format: "bed",
                    isSvs: this.isSvs,
                    features: this.baseFeatures,
                    displayMode: "EXPANDED",
                    height: 220,
                    searchable: false,
                    supportsWholeGenome: true,
                    visibilityWindow: -1,
                    __mdvGetFeatureStyle: (rowIndex: number) => this.getFeatureStyle(rowIndex),
                    __mdvIsFeatureVisible: (rowIndex: number) => this.isFeatureVisible(rowIndex),
                },
                ...((((this.config.igv_config as any)?.tracks as any[]) || []).filter((t) => t.id !== BASE_TRACK_ID)),
            ],
        } as any;

        const safeLocus = toSafeLocusString(this.config.location);
        if (safeLocus) {
            initialConfig.locus = safeLocus;
        }
        initialConfig.__mdvDataStore = this.dataStore;
        return initialConfig;
    }

    attachBrowser(browser: any, baseTrack: any) {
        this.browser = browser;
        this.baseTrack = baseTrack;
       /* if (this.config.feature_label) {
            this.setLabelFunction(this.config.feature_label);
        }
        if (this.config.color_by && typeof this.config.color_by === "string") {
            this.colorByColumn(this.config.color_by);
        }
        this.onDataFiltered();
        */
    }

    getColorOptions() {
        return {
            colorby: "all",
            has_default_color: true,
        };
    }

    colorByDefault() {
        this.config.color_by = undefined;
        if (!this.baseTrack) return;
        this.updateMDVFeatures();
    }

    colorByColumn(column: string) {
        if (!column) return;
        this.config.color_by = column;
        this.colorFunction= this.getColorFunction(column);
        if (!this.baseTrack) return;
        this.updateMDVFeatures();
    }

    setLabelFunction(column: string | null) {
        if (!this.baseTrack) return;
        if (!column || column === "_none") {
            this.config.feature_label = undefined;
        } else {
            this.config.feature_label = column;
        }
        this.updateMDVFeatures();
    }

    getFeatureStyle(rowIndex: number) {
        const labelColumn = this.config.feature_label;
        const name = labelColumn ? `${this.dataStore.getRowText(rowIndex, labelColumn)}` : undefined;
        let color: string | unknown;
        color = this.colorFunction ? this.colorFunction(rowIndex) : "#888888";
        return { name, color };
    }

    isFeatureVisible(rowIndex: number) {
        if (!this.dataStore.isFiltered()) {
            return true;
        }
        return this.dataStore.isRowFiltered(rowIndex);
    }

    async updateMDVFeatures() {
        if (!this.browser || !this.baseTrack) return;
        this.baseTrack.updateMDVFeatures();
    }

    async addTrackFromUrl(name: string, url: string) {
        if (!this.browser) {
            console.warn("IGV browser is not ready; cannot add track yet.");
            return false;
        }
        const trimmedName = name.trim();
        const trimmedUrl = url.trim();
        if (!trimmedName || !trimmedUrl) {
            return false;
        }
        try {
            const trackConfig: Record<string, unknown> = {
                name: trimmedName,
                url: trimmedUrl,
            };

            const loweredPath = getPathWithoutQueryOrHash(trimmedUrl).toLowerCase();
            if (loweredPath.endsWith(".cram")) {
                trackConfig.type = "alignment";
                trackConfig.format = "cram";
                trackConfig.indexURL = appendSuffixBeforeQueryOrHash(trimmedUrl, ".crai");
            } else if (loweredPath.endsWith(".bam")) {
                trackConfig.type = "alignment";
                trackConfig.format = "bam";
                trackConfig.indexURL = appendSuffixBeforeQueryOrHash(trimmedUrl, ".bai");
            }

            await this.browser.loadTrack(trackConfig);
            return true;
        } catch (e) {
            console.error("Failed to add IGV track", e);
            return false;
        }
    }

    onDataFiltered() {
        void this.updateMDVFeatures();
    }

    async onDataHighlighted(data: any) {
        if (data.source === this || !this.locationFields.length || !data?.indexes?.length) {
            return;
        }
        const row = this.dataStore.getRowAsObject(data.indexes[0], this.locationFields) as Record<string, unknown>;
        const chrValue = row[this.locationFields[0]] as string;
        const startValue = (row[this.locationFields[1]]) as number;
        const endValue = (row[this.locationFields[2]]) as number;
      
        const rawLocation = locationFromFieldValues(
            {
                chr: chrValue,
                start: startValue,
                end: endValue,
                chr2: this.isSvs ? row[this.locationFields[3]] as string: undefined,
            },
            this.isSvs,
        );
        if (!rawLocation) return;
        const locations = Array.isArray(rawLocation) ? rawLocation : [rawLocation];
        const nextLocations = locations.map((loc) => applyViewMargins(loc, this.config.view_margins!));
        let nextLocation: IGVBrowserLocation;
        if (nextLocations.length === 1) {
            const first = nextLocations[0];
            if (!first) return;
            nextLocation = first;
        } else {
            const first = nextLocations[0];
            const second = nextLocations[1];
            if (!first || !second) return;
            nextLocation = [first, second];
        }
        const searchLocus = toSafeLocusString(nextLocation);
        const roiSource = locations[0];
        action(() => {
            this.config.location = nextLocation;
        })();

        const features = [{...roiSource,name:"highlight1",color:"blue"}];
        if (nextLocations.length === 2) {
            features.push({chr:nextLocations[1].chr,start:endValue,end:endValue,name:"highlight",color:"blue"});
        }
        if (this.browser?.search && searchLocus) {
            this.setSearchPending(true);
            try {
                await this.browser.search(searchLocus);
                this.browser.clearROIs();
                const highlight ={
                    name:"highlight",
                    type:"annotation",
                    color: "rgba(68, 134, 247, 0.25)",
                    features: features
                };
                this.browser.loadROI(highlight);
            } finally {
                this.setSearchPending(false);
            }
        }
    }

    setLocationFromLocus(chr: string, start: number, end: number) {
        action(() => {
            this.config.location = { chr, start: Math.floor(start), end: Math.floor(end) };
        })();
    }

    getConfig() {
        const c = super.getConfig();
        if (this.browser?.toJSON) {
            const browserJson = this.browser.toJSON();
            const tracks = Array.isArray(browserJson?.tracks)
                ? browserJson.tracks.filter((t: any) => t.id !== BASE_TRACK_ID)
                : [];
            c.igv_config = {
                ...browserJson,
                tracks,
            };
        }
        return c;
    }

    setSize(x: number, y: number) {
        super.setSize(x, y);
        window.dispatchEvent(new Event("resize"));
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        const vm = c.view_margins || { type: "percentage", value: 20 };
        const cols = this.dataStore.getColumnList();
        cols.push({ name: "None", field: "_none" });
        settings.push(
            g({
                label: "Feature Label",
                type: "dropdown",
                values: [cols, "name", "field"],
                current_value: c.feature_label || "_none",
                func: (x: string) => {
                    this.setLabelFunction(x === "_none" ? null : x);
                },
            }),
        );
        settings.push(
            g({
                label: "View margin length",
                type: "text",
                current_value: vm.value.toString(),
                only_update_on_enter: true,
                func: (x: string) => {
                    let n = Number.parseInt(x, 10);
                    n = Number.isNaN(n) ? (vm.type === "percentage" ? 20 : 1000) : n;
                    c.view_margins = { type: vm.type, value: n };
                    const d = this.dataStore.getHighlightedData();
                    if (d != null) {
                        this.onDataHighlighted({ indexes: [d] });
                    }
                },
            }),
        );
        settings.push(
            g({
                label: "View margin type",
                type: "radiobuttons",
                choices: [
                    ["% of feature length", "percentage"],
                    ["absolute (bp)", "absolute"],
                    ["Fixed Length (bp)", "fixed_length"],
                ],
                current_value: vm.type,
                func: (x: any) => {
                    c.view_margins = { type: x, value: vm.value };
                    const d = this.dataStore.getHighlightedData();
                    if (d != null) {
                        this.onDataHighlighted({ indexes: [d] });
                    }
                },
            }),
        );
        settings.push(
            g({
                label: "Max initial features",
                type: "text",
                current_value: String(c.max_initial_features || DEFAULT_MAX_INITIAL_FEATURES),
                only_update_on_enter: true,
                func: (x: string) => {
                    const n = Number.parseInt(x, 10);
                    c.max_initial_features = Number.isNaN(n) || n < 1000 ? DEFAULT_MAX_INITIAL_FEATURES : n;
                },
            }),
        );
        settings.push(
            g({
                type: "folder",
                label: "Track Management",
                current_value: [],
            }),
        );
        const trackSection = settings[settings.length - 1].current_value as any[];
        // Use observable objects so TextComponent (observer) re-renders when current_value is cleared
        const nameSpec = observable.object({ label: "Track Name", type: "text" as const, current_value: "" });
        const urlSpec = observable.object({ label: "Track URL", type: "text" as const, current_value: "" });
        const statusSpec = observable.object({ label: "Status", type: "text" as const, current_value: "" });
        this._pendingNameSpec = nameSpec;
        this._pendingUrlSpec = urlSpec;
        this._trackStatusSpec = statusSpec;
        trackSection.push(nameSpec as any);
        trackSection.push(urlSpec as any);
        trackSection.push(
            g({
                label: "Add Track",
                type: "button",
                current_value: "Add",
                func: async () => {
                    // Use the local specs from this rendered settings panel.
                    // Class-level refs can become stale when settings are rebuilt.
                    const name = nameSpec.current_value.trim();
                    const url = urlSpec.current_value.trim();
                    if (!name || !url) {
                        runInAction(() => { statusSpec.current_value = "Please provide both a Track Name and URL."; });
                        return;
                    }
                    runInAction(() => { statusSpec.current_value = "Adding track..."; });
                    const added = await this.addTrackFromUrl(name, url);
                    if (added) {
                        this.clearPendingTrackFields();
                    } else {
                        runInAction(() => { statusSpec.current_value = `Failed to add track: ${name}`; });
                    }
                },
            }),
        );
        trackSection.push(statusSpec as any);
        return settings;
    }
}

BaseChart.types["igv_browser"] = {
    class: IGVBrowser,
    name: "IGV Browser",
    params: [],
    methodsUsingColumns: ["setLabelFunction"],
    configEntriesUsingColumns: ["feature_label"],
    required: (ds) => {
        return ds.genome?.genomic_location || ds.genome?.svs;
    },
    init(config, ds) {
        const locationFields = getLocationFieldsFromGenome(ds.genome);
        if (!locationFields) return;
        config.param = locationFields;
        if (ds.genome?.genomic_location) {
            const cols = ds.genome.genomic_location.columns;
            config.param = [cols.chr, cols.start, cols.end];
        } else if (ds.genome?.svs) {
            const cols = ds.genome.svs.sv_columns;
            config.param = [cols.chr1, cols.pos1, cols.pos2, cols.chr2,cols.svtype,cols.length];
        }
    },
};

export { BASE_TRACK_ID };
export type { IGVBrowserConfig };
export default IGVBrowser;
