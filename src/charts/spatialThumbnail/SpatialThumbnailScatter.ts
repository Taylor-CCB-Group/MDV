import { createElement } from "react";
import { Deck, OrthographicView } from "@deck.gl/core";
import { ScatterplotLayer } from "@deck.gl/layers";
import { getPostData, getProjectRoot, getProjectURL } from "@/dataloaders/DataLoaderUtil";
import BaseChart from "@/charts/BaseChart";
import { createEl } from "@/utilities/Elements";
import { ImageArray } from "@/webgl/ImageArray";
import { ImageArrayDeckExtension } from "@/webgl/ImageArrayDeckExtension";
import { COLOR_PALLETE } from "@/react/components/avivatorish/constants";
import type { GuiSpecs } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";
import { loadRegionPyramid } from "./loadRegionPyramid";
import { compositePlanes, planeDomain, type ChannelComposite } from "./pyramidLevel";
import { layoutCellPositions } from "./thumbnailPoints";
import {
    N_CELLS_FIELD,
    N_TRANSCRIPTS_FIELD,
    REGION_ID_FIELD,
    bestIdColumn,
    buildSampleCsv,
    summaryDatasourceName,
    summaryFieldName,
    TRANSCRIPT_COUNTS_FIELD,
} from "./sampleSummary";
import { ThumbnailChannelPanel } from "./ChannelPanel";

const DEFAULT_TONE = 0.5;
/** Same cap as the SpatialData viewer's `buildDefaultSelection`. */
const VIEWER_CHANNEL_COUNT = 4;
const LOAD_CONCURRENCY = 2;

export type ThumbnailSelection = { z: number; c: number; t: number };

export type ThumbnailChannels = {
    ids: string[];
    colors: [number, number, number][];
    channelsVisible: boolean[];
    selections: ThumbnailSelection[];
    contrastLimits: [number, number][];
    domains: [number, number][];
};

export type ThumbnailTone = {
    brightness: number[];
    contrast: number[];
};

type RegionEntry = {
    id: string;
    url: string | null;
    /** False when there is no OME-Zarr viv image. Shown as a label, not a second GL context. */
    zarr: boolean;
};

type Point = {
    index: number;
    x: number;
    y: number;
    label: string;
};

function parseStoredColor(value: unknown): [number, number, number] | null {
    if (Array.isArray(value) && value.length >= 3) {
        const rgb = [Number(value[0]), Number(value[1]), Number(value[2])];
        if (rgb.every((part) => Number.isFinite(part))) return [rgb[0], rgb[1], rgb[2]];
    }
    if (typeof value !== "string") return null;
    const hex = value.trim().replace(/^#/, "");
    if (!/^[0-9a-fA-F]{6}$/.test(hex)) return null;
    return [Number.parseInt(hex.slice(0, 2), 16), Number.parseInt(hex.slice(2, 4), 16), Number.parseInt(hex.slice(4, 6), 16)];
}

function paletteColor(index: number): [number, number, number] {
    const color = COLOR_PALLETE[index % COLOR_PALLETE.length];
    return [color[0], color[1], color[2]];
}

/** First channels, viewer palette. One mix for every thumbnail. */
function defaultChannels(count = VIEWER_CHANNEL_COUNT): ThumbnailChannels {
    const n = Math.max(1, count);
    return {
        ids: Array.from({ length: n }, (_, index) => `thumb-ch-${index}`),
        colors: Array.from({ length: n }, (_, index) => paletteColor(index)),
        channelsVisible: Array.from({ length: n }, () => true),
        selections: Array.from({ length: n }, (_, index) => ({ z: 0, c: index, t: 0 })),
        contrastLimits: Array.from({ length: n }, () => [0, 255] as [number, number]),
        domains: Array.from({ length: n }, () => [0, 1] as [number, number]),
    };
}

function defaultTone(count: number): ThumbnailTone {
    return {
        brightness: Array.from({ length: count }, () => DEFAULT_TONE),
        contrast: Array.from({ length: count }, () => DEFAULT_TONE),
    };
}

function isOmeTiff(file: string) {
    return /\.(ome\.)?tiff?$/i.test(file);
}

function regionImageUrl(base: string, viv: { file?: string; url?: string } | undefined): string | null {
    if (!viv) return null;
    if (viv.url) return viv.url;
    if (!viv.file || isOmeTiff(viv.file)) return null;
    const root = getProjectURL(base || "spatial/");
    return root + viv.file.replace(/^\//, "");
}

function labelThumbnail(text: string, size: number): Uint8ClampedArray {
    const canvas = document.createElement("canvas");
    canvas.width = size;
    canvas.height = size;
    const ctx = canvas.getContext("2d");
    if (!ctx) return new Uint8ClampedArray(size * size * 4);
    ctx.fillStyle = "#1c1c1c";
    ctx.fillRect(0, 0, size, size);
    ctx.fillStyle = "#e8e8e8";
    ctx.font = "12px sans-serif";
    const words = text.split(/[_-]/);
    let line = "";
    let y = 16;
    for (const word of words) {
        const next = line ? `${line} ${word}` : word;
        if (ctx.measureText(next).width > size - 8) {
            ctx.fillText(line, 4, y);
            line = word;
            y += 14;
            if (y > size - 4) break;
        } else {
            line = next;
        }
    }
    if (y <= size - 4) ctx.fillText(line, 4, y);
    return ctx.getImageData(0, 0, size, size).data;
}

async function mapPool<T>(items: T[], limit: number, fn: (item: T) => Promise<void>) {
    if (items.length === 0) return;
    let cursor = 0;
    const workers = Array.from({ length: Math.min(limit, items.length) }, async () => {
        while (cursor < items.length) {
            const item = items[cursor];
            cursor += 1;
            await fn(item);
        }
    });
    await Promise.all(workers);
}

class SpatialThumbnailScatter extends BaseChart<any> {
    canvas: HTMLCanvasElement;
    imageArray: ImageArray;
    deck: Deck<OrthographicView>;
    points: Point[] = [];
    regions: RegionEntry[] = [];
    /** region id -> channel index -> scaled plane */
    planes = new Map<string, Map<number, Float32Array>>();
    channelNames: string[] = [];
    channelColors: Array<[number, number, number] | null> = [];
    viewerMixApplied = false;
    textureSize: number;
    gridColumns: number;
    gap: number;
    channelListeners = new Set<() => void>();
    size = 1;
    opacity = 255;
    saturation = 1;
    summaryBusy = false;
    /** Built only while Show points is on. */
    cellPositions: Float32Array | null = null;
    cellPositionsDirty = true;
    disposed = false;
    resizeObserver: ResizeObserver | null = null;
    loadingNote: HTMLDivElement;

    constructor(dataStore: DataStore, div: HTMLDivElement, config: any) {
        super(dataStore, div, config);
        // `initialiseChartConfig` copies config before this runs, so defaults belong on `this.config`.
        this.textureSize = Number(this.config.texture_size) || 128;
        this.gridColumns = Number(this.config.grid_columns) || 8;
        this.gap = this.config.gap === undefined || this.config.gap === null ? 8 : Number(this.config.gap);
        if (!this.config.channels) this.config.channels = defaultChannels();
        if (!this.config.vivLayerProps) this.config.vivLayerProps = defaultTone(this.getChannels().ids.length);

        this.contentDiv.style.position = "relative";
        this.contentDiv.style.height = "100%";
        this.contentDiv.style.minHeight = "0";

        const canvas = (this.canvas = createEl("canvas", {}, this.contentDiv));
        canvas.style.display = "block";
        canvas.style.width = "100%";
        canvas.style.height = "100%";

        this.loadingNote = createEl(
            "div",
            { text: "Loading pyramid levels…", styles: { position: "absolute", color: "white", padding: "4px" } },
            this.contentDiv,
        );

        this.regions = this.collectRegions();
        const maxLayers = this.regions.length;
        this.imageArray = new ImageArray(null, canvas, null, {
            width: this.textureSize,
            height: this.textureSize,
            count: Math.max(1, maxLayers),
        });
        const shown = Math.min(this.regions.length, this.imageArray.depth);
        if (this.regions.length > shown) {
            this.regions = this.regions.slice(0, shown);
        }
        this.points = this.layout(this.regions);
        this.deck = this.createDeck();
        this.resizeObserver = new ResizeObserver(() => this.fitGrid());
        this.resizeObserver.observe(this.canvas);
        if (this.layoutMode() === "xy") void this.ensureSummaryColumns();
        void this.loadInitial();
    }

    setSize(x?: number, y?: number) {
        super.setSize(x, y);
        this.fitGrid();
    }

    collectRegions(): RegionEntry[] {
        const all = this.dataStore.regions?.all_regions ?? {};
        const base = this.dataStore.regions?.avivator?.base_url ?? "spatial/";
        const wanted: string[] | undefined = this.config.regions;
        const ids = wanted?.length ? wanted : Object.keys(all);
        const entries: RegionEntry[] = [];
        for (const id of ids) {
            const region = all[id];
            if (!region?.spatial) continue;
            const url = regionImageUrl(base, region.viv_image);
            entries.push({ id, url, zarr: Boolean(url) });
        }
        return entries;
    }

    layoutMode(): "grid" | "xy" {
        return this.config.layout === "xy" ? "xy" : "grid";
    }

    canEdit() {
        return window.mdv.chartManager.config?.permission === "edit";
    }

    summaryStoreName(): string | null {
        const index = window.mdv.chartManager.dsIndex;
        for (const name of ["samples", "sample_summary"]) {
            if (index[name]?.dataStore?.columnIndex?.[N_CELLS_FIELD]) return name;
        }
        return null;
    }

    summaryDataStore() {
        const name = this.summaryStoreName();
        if (!name) return null;
        return window.mdv.chartManager.dsIndex[name]?.dataStore ?? null;
    }

    applyLayout() {
        this.points = this.layoutMode() === "xy" ? this.layoutFromSummary(this.regions) : this.layout(this.regions);
        this.cellPositionsDirty = true;
        this.syncLayer();
        this.fitGrid();
    }

    layout(regions: RegionEntry[]): Point[] {
        const cell = this.textureSize + this.gap;
        return regions.map((region, index) => {
            const col = index % this.gridColumns;
            const row = Math.floor(index / this.gridColumns);
            return {
                index,
                x: col * cell + this.textureSize / 2,
                y: -(row * cell + this.textureSize / 2),
                label: region.id,
            };
        });
    }

    /** One image per region, placed by a row in the sample summary. */
    layoutFromSummary(regions: RegionEntry[]): Point[] {
        const store = this.summaryDataStore();
        const axes = this.sampleAxes();
        const idField = typeof this.config.sample_id_field === "string" ? this.config.sample_id_field : REGION_ID_FIELD;
        const idColumn = store?.columnIndex[idField];
        const xColumn = axes ? store?.columnIndex[axes[0]] : undefined;
        const yColumn = axes ? store?.columnIndex[axes[1]] : undefined;
        if (!store || !idColumn?.data || !xColumn?.data || !yColumn?.data) return this.layout(regions);
        const byId = new Map<string, { x: number; y: number }>();
        for (let i = 0; i < store.size; i++) {
            const x = Number(xColumn.data[i]);
            const y = Number(yColumn.data[i]);
            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
            byId.set(String(idColumn.getValue(i)), { x, y });
        }
        const finite = [...byId.values()];
        if (finite.length === 0) return this.layout(regions);
        let minX = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;
        for (const value of finite) {
            minX = Math.min(minX, value.x);
            maxX = Math.max(maxX, value.x);
            minY = Math.min(minY, value.y);
            maxY = Math.max(maxY, value.y);
        }
        const spanX = maxX - minX || 1;
        const spanY = maxY - minY || 1;
        const plot = this.textureSize * Math.max(4, Math.ceil(Math.sqrt(regions.length)));
        return regions.map((region, index) => {
            const value = byId.get(region.id);
            const u = value ? (value.x - minX) / spanX : 0.5;
            const v = value ? (value.y - minY) / spanY : 0.5;
            return { index, x: u * plot, y: v * plot, label: region.id };
        });
    }

    regionMembership(): { regionCode: ArrayLike<number> | null; slotCodes: number[] } | null {
        const regionField = this.dataStore.regions?.region_field;
        const regionCol = regionField ? this.dataStore.columnIndex[regionField] : undefined;
        const codes = regionCol?.data;
        const names = regionCol?.values;
        if (codes && names && (regionCol.datatype === "text" || regionCol.datatype === "text16")) {
            return {
                regionCode: codes,
                slotCodes: this.regions.map((region) => names.indexOf(region.id)),
            };
        }
        if (this.regions.length === 1) return { regionCode: null, slotCodes: [] };
        return null;
    }

    createDeck() {
        const cell = this.textureSize + this.gap;
        const cols = Math.min(this.gridColumns, Math.max(1, this.points.length));
        const rows = Math.ceil(this.points.length / cols) || 1;
        return new Deck<OrthographicView>({
            canvas: this.canvas,
            layers: this.layers(),
            views: new OrthographicView({}),
            controller: true,
            initialViewState: {
                target: [(cols * cell) / 2, -(rows * cell) / 2, 0],
                zoom: 0,
            },
            onViewStateChange: ({ viewState }) => {
                this.deck.setProps({ viewState });
            },
            getTooltip: (info) => {
                const point = info.object as Point | undefined;
                if (!info.picked || !point) return null;
                return point.label;
            },
        });
    }

    makeLayer() {
        return new ScatterplotLayer<Point>({
            id: "spatial-thumbnail-scatter",
            data: this.points,
            pickable: true,
            radiusUnits: "pixels",
            getRadius: (this.textureSize / 2) * this.size,
            getPosition: (d) => [d.x, d.y, 0],
            getFillColor: [255, 255, 255, 255],
            opacity: this.opacity / 255,
            // Extension attributes are not on ScatterplotLayer's own props.
            getImageIndex: (d: Point) => d.index,
            getImageAspect: () => 1,
            imageArray: this.imageArray,
            saturation: this.saturation,
            extensions: [new ImageArrayDeckExtension()],
        } as ConstructorParameters<typeof ScatterplotLayer<Point>>[0] & Record<string, unknown>);
    }

    showPoints() {
        return this.config.show_points === true;
    }

    onDataFiltered() {
        this.cellPositionsDirty = true;
        if (this.showPoints() && this.deck) this.syncLayer();
    }

    /** One scatter of cell centroids in the same deck. Absent until the setting is on. */
    layers() {
        const image = this.makeLayer();
        if (!this.showPoints()) return [image];
        const positions = this.ensureCellPositions();
        if (!positions || positions.length === 0) return [image];
        return [image, this.makePointsLayer(positions)];
    }

    ensureCellPositions() {
        if (!this.cellPositionsDirty && this.cellPositions) return this.cellPositions;
        const positions = this.collectCellPositions();
        if (!positions) return null;
        this.cellPositions = positions;
        this.cellPositionsDirty = false;
        return positions;
    }

    collectCellPositions(): Float32Array | null {
        const fields = this.dataStore.regions?.position_fields;
        if (!fields || fields.length < 2) return new Float32Array(0);
        const x = this.dataStore.columnIndex[fields[0]]?.data;
        const y = this.dataStore.columnIndex[fields[1]]?.data;
        const membership = this.regionMembership();
        if (!x || !y || !membership) return x && y ? new Float32Array(0) : null;
        const { regionCode, slotCodes } = membership;
        return layoutCellPositions({
            x,
            y,
            regionCode,
            slotCodes,
            filter: this.dataStore.filterArray,
            centers: this.points,
            size: this.textureSize,
        });
    }

    makePointsLayer(positions: Float32Array) {
        const count = positions.length / 2;
        return new ScatterplotLayer({
            id: "spatial-thumbnail-points",
            data: { length: count },
            pickable: false,
            radiusUnits: "pixels",
            getRadius: 1.5,
            getPosition: (_d: unknown, info: { index: number }) => {
                const i = info.index * 2;
                const position: [number, number, number] = [positions[i], positions[i + 1], 0];
                return position;
            },
            getFillColor: [255, 255, 255, 230],
            parameters: { depthTest: false },
        });
    }

    syncLayer() {
        this.deck?.setProps({ layers: this.layers() });
    }

    /** Fit the images in the canvas. Channel edits do not call this. */
    fitGrid() {
        if (!this.deck || !this.points.length) return;
        if (this.layoutMode() === "xy") {
            this.fitPoints();
            return;
        }
        const cell = this.textureSize + this.gap;
        const cols = Math.min(this.gridColumns, this.points.length);
        const rows = Math.ceil(this.points.length / cols);
        const gridW = cols * cell;
        const gridH = rows * cell;
        const canvasW = this.canvas.clientWidth || this.width;
        const canvasH = this.canvas.clientHeight || this.height;
        if (!canvasW || !canvasH || !gridW || !gridH) return;
        const zoom = Math.log2(Math.min(canvasW / gridW, canvasH / gridH) * 0.92);
        this.deck.setProps({
            viewState: {
                target: [gridW / 2, -gridH / 2, 0],
                zoom: Number.isFinite(zoom) ? zoom : 0,
            },
        });
    }

    fitPoints() {
        let minX = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;
        for (const point of this.points) {
            minX = Math.min(minX, point.x);
            maxX = Math.max(maxX, point.x);
            minY = Math.min(minY, point.y);
            maxY = Math.max(maxY, point.y);
        }
        const pad = this.textureSize / 2 + this.gap;
        const gridW = Math.max(this.textureSize, maxX - minX + pad * 2);
        const gridH = Math.max(this.textureSize, maxY - minY + pad * 2);
        const canvasW = this.canvas.clientWidth || this.width;
        const canvasH = this.canvas.clientHeight || this.height;
        if (!canvasW || !canvasH) return;
        const zoom = Math.log2(Math.min(canvasW / gridW, canvasH / gridH) * 0.92);
        this.deck.setProps({
            viewState: {
                target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
                zoom: Number.isFinite(zoom) ? zoom : 0,
            },
        });
    }

    markChannelsCustom() {
        this.config.channel_mix = "custom";
        this.viewerMixApplied = true;
    }

    /**
     * Match the SpatialData viewer: project `default_channels` when set, otherwise the first
     * four channels in the viewer palette (or the OME color when the store has one).
     */
    adoptViewerChannels() {
        if (this.viewerMixApplied || this.config.channel_mix === "custom") return;
        const names = this.channelNames;
        const available = names.length || this.channelColors.length;
        if (!available) return;
        const previous = new Set(this.getChannels().selections.map((selection) => selection.c));
        this.viewerMixApplied = true;
        const fromProject = this.channelsFromProjectDefaults(names);
        const next = fromProject ?? this.channelsFromImage(Math.min(VIEWER_CHANNEL_COUNT, available));
        this.config.channels = next;
        this.config.vivLayerProps = defaultTone(next.ids.length);
        this.notify();
        const wanted = next.selections.map((selection) => selection.c);
        if (wanted.some((c) => !previous.has(c))) void this.ensureChannels(wanted);
    }

    channelsFromImage(count: number): ThumbnailChannels {
        const channels = defaultChannels(count);
        channels.colors = channels.colors.map((fallback, index) => this.channelColors[index] ?? fallback);
        return channels;
    }

    channelsFromProjectDefaults(names: string[]): ThumbnailChannels | null {
        const defaults = this.dataStore.regions?.avivator?.default_channels;
        if (!Array.isArray(defaults) || defaults.length === 0) return null;
        const slots: { c: number; color: [number, number, number] }[] = [];
        defaults.forEach((entry: { name?: string; color?: unknown }, index: number) => {
            const named = entry?.name ? names.indexOf(entry.name) : -1;
            const c = named >= 0 ? named : index;
            if (names.length && c >= names.length) return;
            const parsed = parseStoredColor(entry?.color);
            slots.push({ c, color: parsed ?? paletteColor(index) });
        });
        if (!slots.length) return null;
        return {
            ids: slots.map((_, index) => `thumb-ch-${index}`),
            colors: slots.map((slot) => slot.color),
            channelsVisible: slots.map(() => true),
            selections: slots.map((slot) => ({ z: 0, c: slot.c, t: 0 })),
            contrastLimits: slots.map(() => [0, 255] as [number, number]),
            domains: slots.map(() => [0, 1] as [number, number]),
        };
    }

    /** Rendered at the top of the chart settings dialog. One mix for every thumbnail. */
    settingsHeader() {
        return createElement(ThumbnailChannelPanel, { chart: this });
    }

    sampleAxes(): [string, string] | null {
        const store = this.summaryDataStore();
        if (!store) return null;
        const choices = this.numericFields(store);
        const preferred = [
            ...choices.filter((field) => field !== N_CELLS_FIELD && field !== N_TRANSCRIPTS_FIELD),
            ...choices.filter((field) => field === N_CELLS_FIELD || field === N_TRANSCRIPTS_FIELD),
        ];
        const savedX = typeof this.config.sample_x === "string" ? this.config.sample_x : "";
        const savedY = typeof this.config.sample_y === "string" ? this.config.sample_y : "";
        const x = choices.includes(savedX) ? savedX : preferred[0];
        const y = choices.includes(savedY) ? savedY : (preferred.find((field) => field !== x) ?? preferred[0]);
        if (!x || !y) return null;
        return [x, y];
    }

    numericFields(store: { getColumnList: (filter: string) => { field?: string; name?: string }[] }) {
        const fields: string[] = [];
        for (const column of store.getColumnList("number")) {
            if (typeof column.field === "string") fields.push(column.field);
        }
        return fields;
    }

    ensureSummaryColumns() {
        const name = this.summaryStoreName();
        const axes = this.sampleAxes();
        if (!name || !axes) return;
        this.config.sample_summary = name;
        this.config.sample_id_field = REGION_ID_FIELD;
        this.config.sample_x = axes[0];
        this.config.sample_y = axes[1];
        const fields = [REGION_ID_FIELD, axes[0], axes[1]];
        const store = this.summaryDataStore();
        const missing = fields.filter((field) => !store?.columnIndex[field]?.data);
        if (missing.length === 0) {
            if (this.layoutMode() === "xy") this.applyLayout();
            return;
        }
        window.mdv.chartManager.loadColumnSet(missing, name, () => {
            if (this.disposed || this.layoutMode() !== "xy") return;
            this.applyLayout();
        });
    }

    summaryPaneMissing(name: string) {
        return !window.mdv.chartManager.viewData?.dataSources?.[name];
    }

    /** Save the view with a pane for the summary, and keep any summary charts getState would drop. */
    async persistSummaryPane(name: string) {
        const cm = window.mdv.chartManager;
        const state = cm.getState();
        const sources = state.view.dataSources;
        if (!sources[name]) {
            const share = 25;
            const scale = (100 - share) / 100;
            for (const entry of Object.values(sources)) {
                if (!entry || typeof entry !== "object" || !("panelWidth" in entry)) continue;
                if (typeof entry.panelWidth === "number") entry.panelWidth *= scale;
            }
            sources[name] = { layout: "absolute", panelWidth: share };
        }
        const kept: Record<string, unknown>[] = [];
        for (const entry of Object.values(cm.charts)) {
            if (entry.dataSource?.name !== name) continue;
            const title = entry.chart?.config?.title;
            if (Array.isArray(title)) entry.chart.config.title = name;
            const config = entry.chart.getConfig();
            const size = config.size;
            if (!Array.isArray(size) || size[0] < 40 || size[1] < 40) {
                config.size = [170, 400];
                config.position = [8, 8];
            }
            kept.push(config);
        }
        if (kept.length) state.view.initialCharts[name] = kept;
        await getPostData(`${getProjectRoot()}/save_state`, state);
    }

    async openSampleSummary(switchToXy: boolean) {
        if (!this.canEdit() || this.summaryBusy) return;
        const existing = this.summaryStoreName();
        if (existing) {
            if (switchToXy) {
                this.config.layout = "xy";
                this.config.sample_summary = existing;
                this.ensureSummaryColumns();
            }
            if (this.summaryPaneMissing(existing)) {
                await this.persistSummaryPane(existing);
                window.location.reload();
                return;
            }
            if (!switchToXy) {
                window.mdv.chartManager.createInfoAlert("Sample summary is in this view", { duration: 2500 });
            }
            return;
        }
        this.summaryBusy = true;
        try {
            if (switchToXy) this.config.layout = "xy";
            const csv = await this.buildSampleCsv();
            const names = window.mdv.chartManager.dataSources.map((source: { name: string }) => source.name);
            const samplesIsSummary = Boolean(window.mdv.chartManager.dsIndex.samples?.dataStore?.columnIndex?.[N_CELLS_FIELD]);
            const name = summaryDatasourceName(names, samplesIsSummary);
            this.config.sample_summary = name;
            this.config.sample_id_field = REGION_ID_FIELD;
            const cm = window.mdv.chartManager;
            try {
                await this.persistSummaryPane(name);
            } catch (error) {
                console.error(error);
            }
            const body = new FormData();
            body.append("name", name);
            body.append("view", cm.viewManager.current_view);
            body.append("file", new Blob([csv], { type: "text/csv" }), `${name}.csv`);
            const response = await fetch(getProjectURL("add_datasource", false), { method: "POST", body });
            if (!response.ok) {
                cm.createInfoAlert((await response.text()) || "Could not create the sample summary", {
                    type: "danger",
                    duration: 4000,
                });
                return;
            }
            window.location.reload();
        } finally {
            this.summaryBusy = false;
        }
    }

    async buildSampleCsv() {
        await this.ensureCountColumns();
        const counts = this.cellCounts();
        const extras = await this.userColumnsByRegion();
        const rows = this.regions.map((region) => {
            const row: Record<string, string | number> = {
                [REGION_ID_FIELD]: region.id,
                [N_CELLS_FIELD]: counts.nCells.get(region.id) ?? 0,
            };
            if (counts.nTranscripts) row[N_TRANSCRIPTS_FIELD] = counts.nTranscripts.get(region.id) ?? 0;
            const extra = extras.get(region.id);
            if (extra) Object.assign(row, extra);
            return row;
        });
        const headers: string[] = [];
        for (const row of rows) {
            for (const key of Object.keys(row)) if (!headers.includes(key)) headers.push(key);
        }
        const numeric = headers.filter(
            (header) => header !== REGION_ID_FIELD && rows.some((row) => typeof row[header] === "number"),
        );
        const preferred = [
            ...numeric.filter((header) => header !== N_CELLS_FIELD && header !== N_TRANSCRIPTS_FIELD),
            ...numeric.filter((header) => header === N_CELLS_FIELD || header === N_TRANSCRIPTS_FIELD),
        ];
        if (!this.config.sample_x && preferred[0]) this.config.sample_x = preferred[0];
        if (!this.config.sample_y) this.config.sample_y = preferred[1] ?? preferred[0];
        return buildSampleCsv(rows);
    }

    ensureCountColumns() {
        const fields: string[] = [];
        const regionField = this.dataStore.regions?.region_field;
        if (regionField && !this.dataStore.columnIndex[regionField]?.data) fields.push(regionField);
        if (this.dataStore.columnIndex[TRANSCRIPT_COUNTS_FIELD] && !this.dataStore.columnIndex[TRANSCRIPT_COUNTS_FIELD]?.data) {
            fields.push(TRANSCRIPT_COUNTS_FIELD);
        }
        if (fields.length === 0) return Promise.resolve();
        return new Promise<void>((resolve) => {
            window.mdv.chartManager.loadColumnSet(fields, this.dataStore.name, () => resolve());
        });
    }

    cellCounts() {
        const nCells = new Map(this.regions.map((region) => [region.id, 0]));
        const transcriptColumn = this.dataStore.columnIndex[TRANSCRIPT_COUNTS_FIELD];
        const transcriptData = transcriptColumn?.data;
        const nTranscripts = transcriptData ? new Map(this.regions.map((region) => [region.id, 0])) : null;
        const membership = this.regionMembership();
        const n = this.dataStore.size;
        if (!membership) return { nCells, nTranscripts };
        if (membership.regionCode === null) {
            const id = this.regions[0]?.id;
            if (!id) return { nCells, nTranscripts };
            nCells.set(id, n);
            if (nTranscripts && transcriptData) {
                let sum = 0;
                for (let i = 0; i < n; i++) {
                    const value = Number(transcriptData[i]);
                    if (Number.isFinite(value)) sum += value;
                }
                nTranscripts.set(id, sum);
            }
            return { nCells, nTranscripts };
        }
        const codeToId = new Map<number, string>();
        membership.slotCodes.forEach((code, slot) => {
            if (code >= 0) codeToId.set(code, this.regions[slot].id);
        });
        for (let i = 0; i < n; i++) {
            const id = codeToId.get(membership.regionCode[i]);
            if (!id) continue;
            nCells.set(id, (nCells.get(id) ?? 0) + 1);
            if (nTranscripts && transcriptData) {
                const value = Number(transcriptData[i]);
                if (Number.isFinite(value)) nTranscripts.set(id, (nTranscripts.get(id) ?? 0) + value);
            }
        }
        return { nCells, nTranscripts };
    }

    async userColumnsByRegion() {
        const extras = new Map<string, Record<string, string | number>>();
        const regionIds = this.regions.map((region) => region.id);
        const ownName = this.dataStore.name;
        let chosen: { name: string; field: string } | null = null;
        let bestScore = 0;
        for (const source of window.mdv.chartManager.dataSources) {
            if (source.name === ownName) continue;
            if (source.dataStore?.columnIndex?.[N_CELLS_FIELD] && (source.name === "samples" || source.name === "sample_summary")) {
                continue;
            }
            const store = source.dataStore;
            if (!store) continue;
            const columns: { field: string; values: string[] }[] = [];
            for (const field of Object.keys(store.columnIndex)) {
                const column = store.columnIndex[field];
                if (!column) continue;
                if ((column.datatype === "text" || column.datatype === "text16") && Array.isArray(column.values)) {
                    columns.push({ field, values: column.values });
                }
            }
            const field = bestIdColumn(regionIds, columns);
            if (!field) continue;
            const values = columns.find((column) => column.field === field)?.values ?? [];
            const ids = new Set(regionIds);
            let score = 0;
            for (const value of values) if (ids.has(value)) score += 1;
            if (score > bestScore) {
                bestScore = score;
                chosen = { name: source.name, field };
            }
        }
        if (!chosen) return extras;
        const store = window.mdv.chartManager.dsIndex[chosen.name]?.dataStore;
        if (!store) return extras;
        const fields = Object.keys(store.columnIndex);
        const missing = fields.filter((field) => !store.columnIndex[field]?.data);
        if (missing.length) {
            await new Promise<void>((resolve) => {
                window.mdv.chartManager.loadColumnSet(missing, chosen.name, () => resolve());
            });
        }
        const idColumn = store.columnIndex[chosen.field];
        if (!idColumn?.data) return extras;
        const rowById = new Map<string, number>();
        for (let i = 0; i < store.size; i++) {
            const id = String(idColumn.getValue(i));
            if (!rowById.has(id)) rowById.set(id, i);
        }
        for (const region of this.regions) {
            const row = rowById.get(region.id);
            if (row === undefined) continue;
            const record: Record<string, string | number> = {};
            for (const field of fields) {
                if (field === chosen.field) continue;
                const column = store.columnIndex[field];
                if (!column?.data) continue;
                const value = column.getValue(row);
                if (value === undefined || value === "missing") continue;
                record[summaryFieldName(field)] = value;
            }
            extras.set(region.id, record);
        }
        return extras;
    }

    getSettings(): GuiSpecs {
        const summary = this.summaryDataStore();
        const axes = this.sampleAxes();
        const settings: GuiSpecs = [
            ...super.getSettings(),
        ];
        if (this.canEdit() || summary) {
            settings.push({
                type: "radiobuttons",
                label: "Layout",
                current_value: this.layoutMode(),
                choices: [
                    ["Grid", "grid"],
                    ["X/Y", "xy"],
                ],
                func: (v) => {
                    if (v === "xy") {
                        void this.openSampleSummary(true);
                        return;
                    }
                    this.config.layout = "grid";
                    this.applyLayout();
                },
            });
        }
        if (summary && axes && this.layoutMode() === "xy") settings.push(...this.summaryColumnSettings(summary, axes));
        if (this.canEdit()) {
            settings.push({
                type: "button",
                label: "Sample summary",
                current_value: null,
                func: () => {
                    void this.openSampleSummary(false);
                },
            });
        }
        settings.push(
            {
                type: "slider",
                label: "Size",
                current_value: this.size,
                min: 0.25,
                max: 4,
                step: 0.25,
                continuous: true,
                func: (v) => {
                    this.size = Number(v);
                    this.syncLayer();
                },
            },
            {
                type: "slider",
                label: "Opacity",
                current_value: this.opacity,
                min: 0,
                max: 255,
                step: 1,
                continuous: true,
                func: (v) => {
                    this.opacity = Number(v);
                    this.syncLayer();
                },
            },
            {
                type: "slider",
                label: "Saturation",
                current_value: this.saturation,
                min: 0,
                max: 1,
                step: 0.05,
                continuous: true,
                func: (v) => {
                    this.saturation = Number(v);
                    this.syncLayer();
                },
            },
            {
                type: "check",
                label: "Show points",
                current_value: this.showPoints(),
                func: (v) => {
                    this.config.show_points = v;
                    if (!v) this.cellPositions = null;
                    this.cellPositionsDirty = true;
                    this.syncLayer();
                },
            },
        );
        return settings;
    }

    summaryColumnSettings(
        store: { getColumnList: (filter: string) => { field?: string; name?: string }[] },
        axes: [string, string],
    ): GuiSpecs {
        const choices: { name: string; value: string }[] = [];
        for (const column of store.getColumnList("number")) {
            if (typeof column.field === "string" && typeof column.name === "string") {
                choices.push({ name: column.name, value: column.field });
            }
        }
        if (choices.length === 0) return [];
        const withCurrent = (field: string) => {
            if (!choices.some((choice) => choice.value === field)) return [{ name: field, value: field }, ...choices];
            return choices;
        };
        return [
            {
                type: "dropdown",
                label: "X",
                current_value: axes[0],
                values: [withCurrent(axes[0]), "name", "value"],
                func: (v) => this.setSampleAxis(0, v),
            },
            {
                type: "dropdown",
                label: "Y",
                current_value: axes[1],
                values: [withCurrent(axes[1]), "name", "value"],
                func: (v) => this.setSampleAxis(1, v),
            },
        ];
    }

    setSampleAxis(axis: 0 | 1, field: string) {
        if (axis === 0) this.config.sample_x = field;
        else this.config.sample_y = field;
        if (this.layoutMode() !== "xy") return;
        this.ensureSummaryColumns();
    }

    selectedChannelIndexes() {
        const channels = this.getChannels();
        return [...new Set(channels.selections.map((selection) => selection.c))];
    }

    async loadInitial() {
        for (const region of this.regions) {
            if (!region.zarr) {
                this.imageArray.uploadRgba(this.points.find((p) => p.label === region.id)?.index ?? 0, labelThumbnail(region.id, this.textureSize));
            }
        }
        this.imageArray.refreshMipmaps();
        await this.ensureChannels(this.selectedChannelIndexes());
        if (!this.disposed) this.loadingNote.textContent = "";
    }

    async ensureChannels(channelIndexes: number[], allowRetry = true) {
        const pending = this.regions.filter((region) => {
            if (!region.zarr || !region.url) return false;
            const have = this.planes.get(region.id);
            return channelIndexes.some((c) => !have?.has(c));
        });
        let done = 0;
        await mapPool(pending, LOAD_CONCURRENCY, async (region) => {
            if (this.disposed || !region.url) return;
            const missing = channelIndexes.filter((c) => !this.planes.get(region.id)?.has(c));
            const loaded = await loadRegionPyramid(region.url, this.textureSize, missing);
            if (this.disposed) return;
            if (!loaded) {
                region.zarr = false;
                const index = this.points.find((p) => p.label === region.id)?.index ?? 0;
                this.imageArray.uploadRgba(index, labelThumbnail(region.id, this.textureSize));
                return;
            }
            const bucket = this.planes.get(region.id) ?? new Map<number, Float32Array>();
            for (const [c, plane] of loaded.planes) bucket.set(c, plane);
            this.planes.set(region.id, bucket);
            if (!this.channelNames.length && loaded.channelNames.length) {
                this.channelNames = loaded.channelNames;
                this.channelColors = loaded.channelColors;
                this.adoptViewerChannels();
            }
            this.fillDomains(bucket);
            done += 1;
            this.loadingNote.textContent = `Loading pyramid levels ${done}/${pending.length}`;
            this.compositeRegion(region.id);
            this.imageArray.refreshMipmaps();
            this.deck.redraw("composite");
        });
        if (this.disposed) return;
        if (allowRetry) {
            const missing = channelIndexes.filter((c) =>
                this.regions.some((region) => region.zarr && region.url && !this.planes.get(region.id)?.has(c)),
            );
            if (missing.length) await this.ensureChannels(missing, false);
        }
        if (this.disposed) return;
        this.imageArray.refreshMipmaps();
        this.deck.redraw("composite");
        this.notify();
    }

    fillDomains(bucket: Map<number, Float32Array>) {
        const channels = this.getChannels();
        let changed = false;
        channels.selections.forEach((selection, index) => {
            const plane = bucket.get(selection.c);
            if (!plane) return;
            const domain = planeDomain(plane);
            const previous = channels.domains[index];
            if (!previous || (previous[0] === 0 && previous[1] === 1)) {
                channels.domains[index] = domain;
                if (channels.contrastLimits[index]?.[0] === 0 && channels.contrastLimits[index]?.[1] === 255) {
                    channels.contrastLimits[index] = domain;
                }
                changed = true;
            }
        });
        if (changed) this.notify();
    }

    compositeRegion(regionId: string) {
        const point = this.points.find((p) => p.label === regionId);
        const bucket = this.planes.get(regionId);
        if (!point || !bucket) return;
        const channels = this.getChannels();
        const tone = this.getTone();
        const draws: ChannelComposite[] = channels.selections.map((selection, index) => ({
            plane: bucket.get(selection.c) ?? new Float32Array(this.textureSize * this.textureSize),
            color: channels.colors[index] ?? [255, 255, 255],
            visible: channels.channelsVisible[index] !== false,
            contrastLimits: channels.contrastLimits[index] ?? [0, 1],
            brightness: tone.brightness[index] ?? DEFAULT_TONE,
            contrast: tone.contrast[index] ?? DEFAULT_TONE,
        }));
        this.imageArray.uploadRgba(point.index, compositePlanes(this.textureSize, draws));
    }

    compositeAll() {
        for (const region of this.regions) {
            if (this.planes.has(region.id)) this.compositeRegion(region.id);
        }
        this.imageArray.refreshMipmaps();
        // A reason forces a frame. Without one, deck waits until the next pointer event.
        this.deck?.redraw("composite");
    }

    getChannels(): ThumbnailChannels {
        return this.config.channels as ThumbnailChannels;
    }

    getTone(): ThumbnailTone {
        const count = this.getChannels().ids.length;
        const tone = (this.config.vivLayerProps ?? {}) as Partial<ThumbnailTone>;
        return {
            brightness: Array.from({ length: count }, (_, i) => tone.brightness?.[i] ?? DEFAULT_TONE),
            contrast: Array.from({ length: count }, (_, i) => tone.contrast?.[i] ?? DEFAULT_TONE),
        };
    }

    getChannelNames() {
        return this.channelNames;
    }

    histogramSlots() {
        const channels = this.getChannels();
        return channels.selections.map((selection, index) => {
            const plane = this.regions
                .map((region) => this.planes.get(region.id))
                .find((bucket) => bucket?.has(selection.c))
                ?.get(selection.c);
            return {
                raster: plane ? { width: this.textureSize, height: this.textureSize, data: plane } : null,
                domain: channels.domains[index] ?? planeDomain(plane ?? []),
            };
        });
    }

    setChannels(patch: Partial<ThumbnailChannels>) {
        this.markChannelsCustom();
        const channels = this.getChannels();
        Object.assign(channels, patch);
        const added = this.selectedChannelIndexes().filter((c) => {
            return this.regions.some((region) => region.zarr && !this.planes.get(region.id)?.has(c));
        });
        this.compositeAll();
        this.notify();
        if (added.length) void this.ensureChannels(added);
    }

    patchTone(patch: Record<string, unknown>) {
        const tone = this.getTone();
        if (Array.isArray(patch.brightness)) tone.brightness = patch.brightness as number[];
        if (Array.isArray(patch.contrast)) tone.contrast = patch.contrast as number[];
        this.config.vivLayerProps = tone;
        this.compositeAll();
        this.notify();
    }

    patchToneAtIndex(index: number, key: "brightness" | "contrast", value: number) {
        const tone = this.getTone();
        const next = key === "brightness" ? [...tone.brightness] : [...tone.contrast];
        next[index] = value;
        this.patchTone({ [key]: next });
    }

    addChannel() {
        this.markChannelsCustom();
        const channels = this.getChannels();
        const index = channels.ids.length;
        const used = new Set(channels.selections.map((selection) => selection.c));
        const nameCount = this.channelNames.length;
        let c = 0;
        const limit = nameCount > 0 ? nameCount : index + 1;
        for (; c < limit; c++) if (!used.has(c)) break;
        if (nameCount > 0 && c >= nameCount) return;
        const domain = channels.domains[0] ?? [0, 1];
        channels.ids = [...channels.ids, `thumb-ch-${index}`];
        channels.colors = [...channels.colors, paletteColor(index)];
        channels.channelsVisible = [...channels.channelsVisible, true];
        channels.selections = [...channels.selections, { z: 0, c, t: 0 }];
        channels.contrastLimits = [...channels.contrastLimits, [...domain] as [number, number]];
        channels.domains = [...channels.domains, [...domain] as [number, number]];
        const tone = this.getTone();
        this.config.vivLayerProps = {
            brightness: [...tone.brightness, DEFAULT_TONE],
            contrast: [...tone.contrast, DEFAULT_TONE],
        };
        this.notify();
        void this.ensureChannels([c]);
    }

    removeChannel(index: number) {
        this.markChannelsCustom();
        const channels = this.getChannels();
        if (channels.ids.length <= 1) return;
        const drop = <T,>(values: T[]) => values.filter((_, i) => i !== index);
        channels.ids = drop(channels.ids);
        channels.colors = drop(channels.colors);
        channels.channelsVisible = drop(channels.channelsVisible);
        channels.selections = drop(channels.selections);
        channels.contrastLimits = drop(channels.contrastLimits);
        channels.domains = drop(channels.domains);
        const tone = this.getTone();
        this.config.vivLayerProps = {
            brightness: drop(tone.brightness),
            contrast: drop(tone.contrast),
        };
        this.compositeAll();
        this.notify();
    }

    subscribe(listener: () => void) {
        this.channelListeners.add(listener);
        return () => {
            this.channelListeners.delete(listener);
        };
    }

    notify() {
        for (const listener of this.channelListeners) listener();
    }

    remove(notify?: boolean) {
        this.disposed = true;
        this.resizeObserver?.disconnect();
        this.deck?.finalize();
        super.remove(notify);
    }

}

BaseChart.types["spatial_thumbnail_scatter"] = {
    class: SpatialThumbnailScatter,
    name: "Image Scatter Plot (Spatial)",
    params: [],
    required: (ds: DataStore) => {
        const regions = ds.regions?.all_regions;
        if (!regions) return false;
        return Object.values(regions).some((region) => Boolean((region as { spatial?: unknown }).spatial));
    },
    init: (config, _dataSource, extraControls) => {
        config.type = "spatial_thumbnail_scatter";
        config.texture_size = Number(extraControls?.texture_size ?? config.texture_size ?? 128);
        config.grid_columns = Number(extraControls?.grid_columns ?? config.grid_columns ?? 8);
        config.gap = Number(extraControls?.gap ?? config.gap ?? 8);
        if (!config.channels) config.channels = defaultChannels();
        if (!config.vivLayerProps) config.vivLayerProps = defaultTone(config.channels.ids.length);
    },
    extra_controls: () => {
        const sizes = [32, 64, 128, 256, 512, 1024].map((s) => ({ name: `${s}`, value: `${s}` }));
        const columns = [4, 6, 8, 12].map((s) => ({ name: `${s}`, value: `${s}` }));
        const gaps = [0, 4, 8, 16].map((s) => ({ name: `${s}`, value: `${s}` }));
        return [
            { type: "dropdown", name: "texture_size", label: "Texture size", values: sizes, defaultVal: "128" },
            { type: "dropdown", name: "grid_columns", label: "Grid columns", values: columns, defaultVal: "8" },
            { type: "dropdown", name: "gap", label: "Gap", values: gaps, defaultVal: "8" },
        ];
    },
};

export { SpatialThumbnailScatter };
export default SpatialThumbnailScatter;
