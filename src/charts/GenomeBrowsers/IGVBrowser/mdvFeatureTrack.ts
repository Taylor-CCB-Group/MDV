type IgvTrackApi = {
    TrackBase: new (...args: any[]) => any;
    createTrack: (...args: any[]) => Promise<any>;
    registerTrackClass: (name: string, klass: any) => void;
};

const { default: igvBase } = await import("igv");
const igv = igvBase as typeof igvBase & IgvTrackApi;
import type { IGVBaseFeature } from "./igvUtils";
import { getStructuralVariantStyle } from "./igvUtils";
type MdvStyle = { name?: string; color?: string } | null;

class MdvFeatureTrack extends igv.TrackBase {
    config: any;
    browser: any;
    trackView: any;
    delegateTrack: any;
    height: number;
    isSvs: boolean;

    constructor(config: any, browser: any) {
        super(config, browser);
        this.config = config;
        this.browser = browser;
        this.height = config.height || 220;
        this.isSvs = Boolean(config.isSvs);
    }

    async postInit() {
        this.delegateTrack = await igv.createTrack(
            {
                type: "annotation",
                format: "bed",
                features: this.config.features || [],
                name: this.config.name || "MDV Features",
                displayMode: this.config.displayMode || "EXPANDED",
                searchable: false,
            },
            this.browser,
        );
    }

    private getFeature(feature: IGVBaseFeature) {
        return feature._f || feature;
    }

    private getDrawTop(options: any): number {
        if (Number.isFinite(options?.pixelTop)) return Number(options.pixelTop);
        if (Number.isFinite(options?.top)) return Number(options.top);
        return 0;
    }

    private getDrawHeight(options: any): number {
        if (Number.isFinite(options?.pixelHeight) && Number(options.pixelHeight) > 0) {
            return Number(options.pixelHeight);
        }
        if (Number.isFinite(options?.height) && Number(options.height) > 0) {
            return Number(options.height);
        }
        return this.height;
    }

    private getSvLaneIndex(svtype: unknown): number {
        const normalized = typeof svtype === "string" ? svtype.trim().toUpperCase() : "";
        switch (normalized) {
            case "DEL":
                return 0;
            case "INS":
                return 1;
            case "INV":
                return 2;
            case "DUP":
                return 3;
            case "TRA":
                return 4;
            case "BND":
                return 5;
            default:
                return 6;
        }
    }

    private isAllView(options: any): boolean {
        const chr = options?.referenceFrame?.chr ?? options?.referenceFrame?.chrName;
        return typeof chr === "string" && chr.toLowerCase() === "all";
    }

    private getRenderRange(feature: IGVBaseFeature, options: any): { start: number; end: number } {
        const proxy = this.getFeature(feature);
        let start = Number(proxy.start);
        let end = Number(proxy.end);

        if (!Number.isFinite(start) || !Number.isFinite(end)) {
            return { start: 0, end: 0 };
        }

        if (!this.isAllView(options)) {
            return {
                start: Math.min(start, end),
                end: Math.max(start, end),
            };
        }

        const genome = this.browser?.genome;
        const chr = proxy.chr;
        if (genome && typeof genome.getGenomeCoordinate === "function" && typeof chr === "string") {
            const gStart = Number(genome.getGenomeCoordinate(chr, start));
            const gEnd = Number(genome.getGenomeCoordinate(chr, end));
            if (Number.isFinite(gStart) && Number.isFinite(gEnd)) {
                start = gStart;
                end = gEnd;
            }
        }

        return {
            start: Math.min(start, end),
            end: Math.max(start, end),
        };
    }

    private drawSvGlyph(
        ctx: CanvasRenderingContext2D,
        feature: IGVBaseFeature,
        bpStart: number,
        bpPerPixel: number,
        top: number,
        height: number,
        pixelWidth: number,
        stackIdx: number = 0,
        renderStart?: number,
        renderEnd?: number,
    ) {
        const proxy = this.getFeature(feature);
        const normalizedSvType = typeof proxy.svtype === "string" ? proxy.svtype.trim().toUpperCase() : "";
        const style = getStructuralVariantStyle(proxy.svtype);
        const color = proxy.color || style.fillStyle;
        const stroke = proxy.color || style.strokeStyle;
        const featureStart = Number.isFinite(renderStart) ? Number(renderStart) : proxy.start;
        const featureEnd = Number.isFinite(renderEnd) ? Number(renderEnd) : proxy.end;
        const startPx = Math.round((featureStart - bpStart) / bpPerPixel);
        const endPx = Math.round((featureEnd - bpStart) / bpPerPixel);
        const viewRight = Math.max(1, Math.round(pixelWidth));
        if (endPx < 0 || startPx > viewRight) {
            return;
        }
        const x1 = Math.min(startPx, endPx);
        const x2 = Math.max(startPx, endPx);
        const centerX = Math.round((x1 + x2) / 2);
        const laneCount = 7;
        const laneIndex = this.getSvLaneIndex(proxy.svtype);
        const laneHeight = Math.max(10, Math.floor(height / laneCount));
        // Stacking: split lane into 3 sub-lanes if needed
        const stackMax = 3;
        const stackHeight = Math.max(6, Math.floor(laneHeight / stackMax));
        const laneTop = top + laneIndex * laneHeight;
        // Offset midY by stackIdx (-1, 0, 1 for 3 stacks)
        const midY = laneTop + Math.round(laneHeight / 2) + (stackIdx - 1) * Math.floor(stackHeight / 2);
        const halfH = Math.max(3, Math.round(laneHeight * 0.25));
        const spanLeft = Math.max(0, Math.min(viewRight, x1));
        const clampedRight = Math.max(0, Math.min(viewRight, x2));
        const spanRight = Math.max(spanLeft + 1, clampedRight);
        const spanWidth = Math.max(4, spanRight - spanLeft);

        ctx.save();
        ctx.lineWidth = 2;
        ctx.fillStyle = color;
        ctx.strokeStyle = stroke;

        if (style.glyph === "breakend") {
            ctx.setLineDash([5, 3]);
            ctx.beginPath();
            ctx.moveTo(spanLeft, midY);
            ctx.lineTo(spanRight, midY);
            ctx.stroke();
            ctx.setLineDash([]);
            this.drawCircle(ctx, spanLeft, midY, 4, stroke, color);
            this.drawCircle(ctx, spanRight, midY, 4, stroke, color);
        } else if (style.glyph === "inversion") {
            ctx.beginPath();
            ctx.moveTo(spanLeft, midY);
            ctx.lineTo(spanRight, midY);
            ctx.stroke();
            this.drawArrowHead(ctx, spanLeft, midY, 6, stroke, false);
            this.drawArrowHead(ctx, spanRight, midY, 6, stroke, true);
        } else {
            ctx.beginPath();
            ctx.moveTo(spanLeft, midY);
            ctx.lineTo(spanRight, midY);
            ctx.stroke();

            switch (style.glyph) {
                case "capped_line":
                    this.drawDeletionTicks(ctx, spanLeft, spanRight, midY, halfH, stroke);
                    break;
                case "diamond":
                    this.drawDiamond(ctx, centerX, midY, Math.max(5, Math.min(10, spanWidth / 3)), color, stroke);
                    break;
                case "triangle":
                    if (normalizedSvType === "INS") {
                        this.drawThinInsertionTriangle(ctx, centerX, midY, 3, 9, color, stroke);
                    } else {
                        this.drawTriangle(ctx, centerX, midY - halfH / 2, Math.max(5, Math.min(10, spanWidth / 3)), color, stroke);
                    }
                    break;
                case "double_bar":
                    this.drawDoubleBar(ctx, centerX, midY, Math.max(6, Math.min(12, spanWidth / 2)), halfH, color, stroke);
                    break;
                default:
                    this.drawDiamond(ctx, centerX, midY, Math.max(5, Math.min(10, spanWidth / 3)), color, stroke);
                    break;
            }
        }

        ctx.restore();
    }

    private drawDiamond(ctx: CanvasRenderingContext2D, x: number, y: number, radius: number, fillStyle: string, strokeStyle: string) {
        ctx.save();
        ctx.fillStyle = fillStyle;
        ctx.strokeStyle = strokeStyle;
        ctx.beginPath();
        ctx.moveTo(x, y - radius);
        ctx.lineTo(x + radius, y);
        ctx.lineTo(x, y + radius);
        ctx.lineTo(x - radius, y);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
        ctx.restore();
    }

    private drawTriangle(ctx: CanvasRenderingContext2D, x: number, y: number, radius: number, fillStyle: string, strokeStyle: string) {
        ctx.save();
        ctx.fillStyle = fillStyle;
        ctx.strokeStyle = strokeStyle;
        ctx.beginPath();
        ctx.moveTo(x, y - radius);
        ctx.lineTo(x + radius, y + radius);
        ctx.lineTo(x - radius, y + radius);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
        ctx.restore();
    }

    private drawThinInsertionTriangle(
        ctx: CanvasRenderingContext2D,
        x: number,
        y: number,
        halfWidth: number,
        height: number,
        fillStyle: string,
        strokeStyle: string,
    ) {
        ctx.save();
        ctx.fillStyle = fillStyle;
        ctx.strokeStyle = strokeStyle;
        ctx.beginPath();
        ctx.moveTo(x, y - height / 2);
        ctx.lineTo(x + halfWidth, y + height / 2);
        ctx.lineTo(x - halfWidth, y + height / 2);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
        ctx.restore();
    }

    private drawDoubleBar(ctx: CanvasRenderingContext2D, x: number, y: number, width: number, halfHeight: number, fillStyle: string, strokeStyle: string) {
        ctx.save();
        ctx.fillStyle = fillStyle;
        ctx.strokeStyle = strokeStyle;
        const barWidth = Math.max(2, Math.round(width / 6));
        ctx.fillRect(x - width / 2, y - halfHeight, barWidth, halfHeight * 2);
        ctx.fillRect(x + width / 2 - barWidth, y - halfHeight, barWidth, halfHeight * 2);
        ctx.strokeRect(x - width / 2, y - halfHeight, barWidth, halfHeight * 2);
        ctx.strokeRect(x + width / 2 - barWidth, y - halfHeight, barWidth, halfHeight * 2);
        ctx.restore();
    }

    private drawDeletionTicks(ctx: CanvasRenderingContext2D, x1: number, x2: number, y: number, halfHeight: number, strokeStyle: string) {
        ctx.save();
        ctx.strokeStyle = strokeStyle;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(x1, y - halfHeight);
        ctx.lineTo(x1, y + halfHeight);
        ctx.moveTo(x2, y - halfHeight);
        ctx.lineTo(x2, y + halfHeight);
        ctx.stroke();
        ctx.restore();
    }

    private drawArrowHead(ctx: CanvasRenderingContext2D, x: number, y: number, size: number, strokeStyle: string, rightFacing: boolean) {
        ctx.save();
        ctx.strokeStyle = strokeStyle;
        ctx.beginPath();
        if (rightFacing) {
            ctx.moveTo(x - size, y - size / 2);
            ctx.lineTo(x, y);
            ctx.lineTo(x - size, y + size / 2);
        } else {
            ctx.moveTo(x + size, y - size / 2);
            ctx.lineTo(x, y);
            ctx.lineTo(x + size, y + size / 2);
        }
        ctx.stroke();
        ctx.restore();
    }

    private drawCircle(ctx: CanvasRenderingContext2D, x: number, y: number, radius: number, strokeStyle: string, fillStyle: string) {
        ctx.save();
        ctx.fillStyle = fillStyle;
        ctx.strokeStyle = strokeStyle;
        ctx.beginPath();
        ctx.arc(x, y, radius, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
        ctx.restore();
    }

    async getFeatures(chr: string, start: number, end: number, bpPerPixel: number) {
        if (!this.delegateTrack) {
            return [];
        }
        const baseFeatures = (await this.delegateTrack.getFeatures(chr, start, end, bpPerPixel)) || [];
        const styleGetter = this.config?.__mdvGetFeatureStyle;
        const filterFn = this.config?.__mdvIsFeatureVisible;

        if (typeof styleGetter !== "function" && typeof filterFn !== "function") {
            return baseFeatures;
        }

     const out = baseFeatures
        .filter((feature: IGVBaseFeature) => {
            const rowIndex = Number.parseInt(`${feature?.id ?? ""}`, 10);
            if (!Number.isFinite(rowIndex)) return false;
            if (typeof filterFn === "function" && !filterFn(rowIndex)) return false;
            return true;
        })
        .map((feature: IGVBaseFeature) => {
            const proxy = feature._f?feature._f:feature;//sometimes igv wraps features in an object with a _f property
            const rowIndex = Number.parseInt(`${feature?.id ?? ""}`, 10);
            if (typeof styleGetter === "function") {
                const style: MdvStyle = styleGetter(rowIndex);
                if (style?.name !== undefined) {
                    proxy.name = style.name;
                } else {
                    proxy.name="";
                }
                if (style?.color !== undefined) {
                    proxy.color = style.color;
                } else {
                    delete proxy.color;
                }
            }
            return feature;
        });
            return out;
        }

    draw(options: any) {
        if (!this.isSvs) {
            if (this.delegateTrack?.draw) {
                return this.delegateTrack.draw(options);
            }
            return;
        }

        const ctx = options.context;
        const features = options.features || [];
        const bpStart = options.bpStart;
        const bpEnd = bpStart + options.pixelWidth * options.bpPerPixel + 1;
        const top = this.getDrawTop(options);
        const height = this.getDrawHeight(options);
        const styleGetter = this.config?.__mdvGetFeatureStyle;
        const filterFn = this.config?.__mdvIsFeatureVisible;

        ctx.save();
        ctx.clearRect(0, top, options.pixelWidth, height);

        // 1. Group features by SV type (lane)
        const laneCount = 7;
        type DrawFeature = {
            feature: IGVBaseFeature;
            renderStart: number;
            renderEnd: number;
        };
        const featuresByLane = Array.from({ length: laneCount }, () => [] as DrawFeature[]);
        for (const feature of features as IGVBaseFeature[]) {
            const proxy = this.getFeature(feature);
            const rowIndex = Number.parseInt(`${proxy?.id ?? ""}`, 10);
            if (!Number.isFinite(rowIndex)) continue;
            if (typeof filterFn === "function" && !filterFn(rowIndex)) continue;
            const { start: renderStart, end: renderEnd } = this.getRenderRange(proxy, options);
            if (renderEnd < bpStart || renderStart > bpEnd) continue;
            const laneIdx = this.getSvLaneIndex(proxy.svtype);
            featuresByLane[laneIdx].push({
                feature: proxy,
                renderStart,
                renderEnd,
            });
        }

        // 2. For each lane, sort and assign stack index for overlaps.
        // INS and breakpoints (TRA/BND) stay on a single centered level.
        const stackingInfo = new WeakMap<IGVBaseFeature, number>();
        for (let laneIdx = 0; laneIdx < laneCount; ++laneIdx) {
            const laneFeatures = featuresByLane[laneIdx];
            laneFeatures.sort((a, b) => (a.renderStart - b.renderStart) || (a.renderEnd - b.renderEnd));
            const singleLevelLane = laneIdx === 1 || laneIdx === 4 || laneIdx === 5; // INS, TRA, BND
            if (singleLevelLane) {
                for (const item of laneFeatures) {
                    stackingInfo.set(item.feature, 1); // centered row in drawSvGlyph
                }
                continue;
            }
            const active: { end: number, stack: number }[] = [];
            for (const item of laneFeatures) {
                for (let i = active.length - 1; i >= 0; --i) {
                    if (active[i].end < item.renderStart) active.splice(i, 1);
                }
                let stackIdx = 0;
                const used = new Set(active.map(a => a.stack));
                while (used.has(stackIdx) && stackIdx < 3) stackIdx++;
                if (stackIdx >= 3) stackIdx = 2;
                stackingInfo.set(item.feature, stackIdx);
                active.push({ end: item.renderEnd, stack: stackIdx });
            }
        }

        // 3. Draw features, offsetting by stack index (pass as last arg for now)
        for (let laneIdx = 0; laneIdx < laneCount; ++laneIdx) {
            const laneFeatures = featuresByLane[laneIdx];
            for (const item of laneFeatures) {
                const proxy = item.feature;
                const rowIndex = Number.parseInt(`${proxy?.id ?? ""}`, 10);
                if (typeof styleGetter === "function") {
                    const style: MdvStyle = styleGetter(rowIndex);
                    if (style?.name !== undefined) {
                        proxy.name = style.name;
                    }
                    if (style?.color !== undefined) {
                        proxy.color = style.color;
                    }
                }
                const stackIdx = stackingInfo.get(proxy) ?? 0;
                this.drawSvGlyph(
                    ctx,
                    proxy,
                    bpStart,
                    options.bpPerPixel,
                    top,
                    height,
                    options.pixelWidth,
                    stackIdx,
                    item.renderStart,
                    item.renderEnd,
                );
            }
        }

        ctx.restore();
        return top + height;
    }

    clickedFeatures(clickState: any) {
        if (!this.isSvs || !this.delegateTrack?.featureSource?.featureCache) {
            return this.delegateTrack?.clickedFeatures?.(clickState) || [];
        }

        const { genomicLocation, y, referenceFrame } = clickState || {};
        if (!referenceFrame || !Number.isFinite(genomicLocation)) {
            return [];
        }

        const bpPerPixel = referenceFrame.bpPerPixel || 1;
        const bpLeeway = Math.max(1, Math.round(4 * bpPerPixel));
        const laneCount = 7;
        const laneHeight = Math.max(1, Math.floor(this.getDrawHeight(clickState) / laneCount));
        const clickedLane = Math.max(0, Math.min(laneCount - 1, Math.floor((Number(y) || 0) / laneHeight)));

        const cache = this.delegateTrack.featureSource.featureCache;
        let candidates: IGVBaseFeature[] = [];
        const isAll = typeof referenceFrame?.chr === "string" && referenceFrame.chr.toLowerCase() === "all";

        if (isAll) {
            const genome = this.browser?.genome;
            if (genome && typeof genome.getChromosomeCoordinate === "function") {
                const chrCoord = genome.getChromosomeCoordinate(genomicLocation);
                const queryChr = chrCoord?.chr;
                const queryPos = Number(chrCoord?.position);
                if (typeof queryChr === "string" && Number.isFinite(queryPos)) {
                    candidates = cache.queryFeatures(
                        queryChr,
                        queryPos - bpLeeway,
                        queryPos + bpLeeway,
                    ) || [];
                }
            }

            // Fallback in case this feature cache stores all-view records directly.
            if (candidates.length === 0) {
                candidates = cache.queryFeatures(
                    referenceFrame.chr,
                    genomicLocation - bpLeeway,
                    genomicLocation + bpLeeway,
                ) || [];
            }
        } else {
            candidates = cache.queryFeatures(
                referenceFrame.chr,
                genomicLocation - bpLeeway,
                genomicLocation + bpLeeway,
            ) || [];
        }

        const hits = candidates.filter((feature: IGVBaseFeature) => {
            const proxy = this.getFeature(feature);
            return this.getSvLaneIndex(proxy.svtype) === clickedLane;
        });
        //remove any that have been filtered out by the user config 
        const filterFn = this.config?.__mdvIsFeatureVisible;
        const finalHits = typeof filterFn === "function"
            ? hits.filter(feature => {
                const proxy = this.getFeature(feature);
                const rowIndex = Number.parseInt(`${proxy?.id ?? ""}`, 10);
                return Number.isFinite(rowIndex) && filterFn(rowIndex);
            })
            : hits;        

        if (finalHits.length > 1) {
            finalHits.sort((a: IGVBaseFeature, b: IGVBaseFeature) => {
                const rangeA = this.getRenderRange(a, clickState);
                const rangeB = this.getRenderRange(b, clickState);
                const posDiffA = Math.min(
                    Math.abs(rangeA.start - genomicLocation),
                    Math.abs(rangeA.end - genomicLocation),
                );
                const posDiffB = Math.min(
                    Math.abs(rangeB.start - genomicLocation),
                    Math.abs(rangeB.end - genomicLocation),
                );
                return posDiffA - posDiffB;
            });
        }

        return finalHits;
    }
    
    popupData(clickState: any) {
        return this.clickedFeatures(clickState);
    }

    menuItemList() {
        const delegateItems = this.delegateTrack?.menuItemList?.() || [];
        if (delegateItems.length === 0) {
            return [];
        }

        const delegateTrack = this.delegateTrack;
        delegateTrack.trackView = this.trackView;
        delegateTrack.browser = this.browser;

        return delegateItems.map((item: any) => this.adaptMenuItem(item, delegateTrack));
    }

    adaptMenuItem(item: any, delegateTrack: any) {
        if (typeof item !== "object" || item === null) {
            return item;
        }

        const wrapped = { ...item };

        if (typeof item.click === "function") {
            wrapped.click = (e: any) => {
                item.click.call(delegateTrack, e);
            };
        }

        if (typeof item.dialog === "function") {
            wrapped.dialog = (e: any) => {
                item.dialog.call(delegateTrack, e);
            };
        }

        return wrapped;
    }

    async updateMDVFeatures() {
        const viewports = this.trackView?.viewports;
        if (!viewports || viewports.length === 0) {
            return;
        }
        await Promise.all(viewports.map(async (viewport: any) => {
            viewport.featureCache = undefined;
            await viewport.loadFeatures();
            viewport.repaint();
        }));
    }

    /*refresh() {
        const viewport = this.trackView?.viewports?.[0];
        viewport?.repaint?.();
    }*/

    get defaultHeight() {
        return this.height;
    }
}

igv.registerTrackClass("mdv_feature_track", MdvFeatureTrack);

class StandaloneCramTrack extends igv.TrackBase {
    config: any;
    browser: any;
    trackView: any;
    delegateTrack: any;

    constructor(config: any, browser: any) {
        super(config, browser);
        this.config = config;
        this.browser = browser;
        this.config.format = "cram"; // Explicitly set format to CRAM
        if (!config.indexURL && config.url) {
            // Auto-generate CRAI URL preserving query and hash
            try {
                const url = new URL(config.url, window.location.href);
                url.pathname = url.pathname + '.crai';
                this.config.indexURL = url.href;
            } catch (e) {
                // Fallback for relative URLs or parsing errors
                this.config.indexURL = `${config.url}.crai`;
            }
        } else {
            this.config.indexURL = config.indexURL;
        }
    }

    async postInit() {
        this.delegateTrack = await igv.createTrack(
            {
                type: "alignment",
                format: this.config.format,
                url: this.config.url,
                indexURL: this.config.indexURL,
                name: this.config.name || "Standalone CRAM Track",
                displayMode: this.config.displayMode || "EXPANDED",
                searchable: false,
            },
            this.browser,
        );
    }

    async getFeatures(chr: string, start: number, end: number) {
        return this.delegateTrack.getFeatures(chr, start, end);
    }

    draw(options: any) {
        this.delegateTrack.draw(options);
    }
}

igv.registerTrackClass("standalone_cram_track", StandaloneCramTrack);