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

type SupplementalSegment = {
    chr: string;
    start: number;
    strand: "+" | "-";
    lenOnRef: number;
    readStart: number;
    readEnd: number;
};





function parseCigarOps(cigar: unknown): Array<{ len: number; op: string }> {
    if (typeof cigar !== "string" || !cigar) return [];
    const out: Array<{ len: number; op: string }> = [];
    const matches = cigar.matchAll(/(\d+)([MIDNSHP=X])/g);
    for (const m of matches) {
        const len = Number.parseInt(m[1] || "0", 10);
        const op = m[2] || "";
        if (!Number.isFinite(len) || len <= 0 || !op) continue;
        out.push({ len, op });
    }
    return out;
}

const TERMINAL_MATCH_TOLERANCE_BP = 200;
const MIN_ROLE_CLIP_BP = 50;

function getTerminalClipTotals(cigar: unknown): { leftClip: number; rightClip: number } {
    const ops = parseCigarOps(cigar);
    if (ops.length === 0) return { leftClip: 0, rightClip: 0 };

    let leftClip = 0;
    let leftEdgeMatch = 0;
    let i = 0;
    for (; i < ops.length; i++) {
        const { op, len } = ops[i];
        if (op === "H" || op === "S") {
            leftClip += len;
        } else {
            break;
        }
    }
    // If edge op is match, compare it against immediately-adjacent clips behind it.
    // This handles patterns like 624M10475H (clip dominates).
    if (i < ops.length) {
        const edge = ops[i];
        if (
            edge &&
            (edge.op === "M" || edge.op === "=" || edge.op === "X")
        ) {
            leftEdgeMatch = edge.len;
            let k = i + 1;
            let adjacentClip = 0;
            for (; k < ops.length; k++) {
                const { op, len } = ops[k];
                if (op === "H" || op === "S") {
                    adjacentClip += len;
                } else {
                    break;
                }
            }
            if (adjacentClip > 0) {
                leftClip += adjacentClip;
            } else if (edge.len <= TERMINAL_MATCH_TOLERANCE_BP) {
                // No immediate clip after edge match: keep previous tolerance behavior.
                i += 1;
                for (; i < ops.length; i++) {
                    const { op, len } = ops[i];
                    if (op === "H" || op === "S") {
                        leftClip += len;
                    } else {
                        break;
                    }
                }
            }
        }
    }
    if (leftClip > 0 && leftEdgeMatch > 0) {
        leftClip = Math.max(leftClip, leftEdgeMatch);
    }

    let rightClip = 0;
    let rightEdgeMatch = 0;
    let j = ops.length - 1;
    for (; j >= 0; j--) {
        const { op, len } = ops[j];
        if (op === "H" || op === "S") {
            rightClip += len;
        } else {
            break;
        }
    }
    if (j >= 0) {
        const edge = ops[j];
        if (
            edge &&
            (edge.op === "M" || edge.op === "=" || edge.op === "X")
        ) {
            rightEdgeMatch = edge.len;
            let k = j - 1;
            let adjacentClip = 0;
            for (; k >= 0; k--) {
                const { op, len } = ops[k];
                if (op === "H" || op === "S") {
                    adjacentClip += len;
                } else {
                    break;
                }
            }
            if (adjacentClip > 0) {
                rightClip += adjacentClip;
            } else if (edge.len <= TERMINAL_MATCH_TOLERANCE_BP) {
                // No immediate clip before edge match: keep previous tolerance behavior.
                j -= 1;
                for (; j >= 0; j--) {
                    const { op, len } = ops[j];
                    if (op === "H" || op === "S") {
                        rightClip += len;
                    } else {
                        break;
                    }
                }
            }
        }
    }
    if (rightClip > 0 && rightEdgeMatch > 0) {
        rightClip = Math.max(rightClip, rightEdgeMatch);
    }
    return { leftClip, rightClip };
}

function getTerminalEdgeMatchTotals(cigar: unknown): { leftMatch: number; rightMatch: number } {
    const ops = parseCigarOps(cigar);
    if (ops.length === 0) return { leftMatch: 0, rightMatch: 0 };

    let leftMatch = 0;
    for (let i = 0; i < ops.length; i++) {
        const { op, len } = ops[i];
        if (op === "H" || op === "S") continue;
        if (op === "M" || op === "=" || op === "X") leftMatch = len;
        break;
    }

    let rightMatch = 0;
    for (let i = ops.length - 1; i >= 0; i--) {
        const { op, len } = ops[i];
        if (op === "H" || op === "S") continue;
        if (op === "M" || op === "=" || op === "X") rightMatch = len;
        break;
    }

    return { leftMatch, rightMatch };
}

function computeOriginalReadIntervalFromCigar(cigar: unknown, strandPlus: boolean): { readStart: number; readEnd: number } | null {
    const ops = parseCigarOps(cigar);
    if (ops.length === 0) return null;

    const { leftClip, rightClip } = getTerminalClipTotals(cigar);

    let alignedQueryLen = 0;
    for (const { len, op } of ops) {
        if (op === "M" || op === "I" || op === "=" || op === "X") {
            alignedQueryLen += len;
        }
    }
    if (alignedQueryLen <= 0) return null;

    // Coordinates on the original read sequence (0-based, end-exclusive).
    let qStart = leftClip;
    let qEnd = qStart + alignedQueryLen;
    const readLen = leftClip + alignedQueryLen + rightClip;
    if (readLen <= 0) return null;

    if (!strandPlus) {
        // Reverse-strand CIGAR clip sides are mirrored relative to original read orientation.
        const rStart = Math.max(0, readLen - qEnd);
        const rEnd = Math.max(rStart, readLen - qStart);
        qStart = rStart;
        qEnd = rEnd;
    }

    return { readStart: qStart, readEnd: qEnd };
}

type EndpointSide = "start" | "end";

function endpointCoord(start: number, end: number, strandPlus: boolean, side: EndpointSide): number {
    if (side === "start") {
        return strandPlus ? start : end;
    }
    return strandPlus ? end : start;
}



class MdvSplitAlignmentTrack extends igv.TrackBase {
    config: any;
    browser: any;
    trackView: any;
    delegateTrack: any;
    private _interOverlay: HTMLCanvasElement | null = null;
    private _activeReadName: string | null = null;
    private _hoverBoundViewports = new WeakSet<any>();
    private _hoverEvalPending = false;
    private _lastHoverEvalTs = 0;
    private _interactionVersion = 0;
    private _syncedInteractionVersion = 0;
    private _postFetchRepaintPending = false;
    private _dragOverlayRafPending = false;
    private _scrollOverlayRafPending = false;
    private _boundTrackDragHandler?: () => void;
    private _boundLocusChangeHandler?: () => void;

    constructor(config: any, browser: any) {
        config.visibilityWindow=100000;
        super(config, browser);
        this.config = config;
        this.browser = browser;
    }

    async postInit() {
        this.delegateTrack = await igv.createTrack(
            {
                ...this.config,
                type: "alignment",
                format: this.config.format,
                url: this.config.url,
                indexURL: this.config.indexURL,
                name: this.config.name || "Split Reads",
                searchable: false,
                // Keep rendering minimal: we care about blocks/gaps/splits, not base mismatch detail.
                showMismatches: false,
                showAllBases: false,
                showCoverage: false,
                hideSmallIndels: true,
                indelSizeThreshold: 5,
            },
            this.browser,
        );
     
        
       
        
        
        this._boundTrackDragHandler = this.onTrackDrag.bind(this);
        this.browser?.on?.("trackdrag", this._boundTrackDragHandler);
        this._boundLocusChangeHandler = this.onLocusChange.bind(this);
        this.browser?.on?.("locuschange", this._boundLocusChangeHandler);
    }



    private onTrackDrag() {
        if (this._dragOverlayRafPending) return;
        this._dragOverlayRafPending = true;
        const run = () => {
            this._dragOverlayRafPending = false;
            this.redrawInterFrameOverlayOnly();
        };
        if (typeof window !== "undefined" && typeof window.requestAnimationFrame === "function") {
            window.requestAnimationFrame(run);
        } else {
            setTimeout(run, 0);
        }
    }

    private scheduleOverlayRedrawFromScroll() {
        if (this._scrollOverlayRafPending) return;
        this._scrollOverlayRafPending = true;
        const run = () => {
            this._scrollOverlayRafPending = false;
            this.redrawInterFrameOverlayOnly();
        };
        if (typeof window !== "undefined" && typeof window.requestAnimationFrame === "function") {
            window.requestAnimationFrame(run);
        } else {
            setTimeout(run, 0);
        }
    }

    private onLocusChange() {
       
        // Ignore drag-time locus updates; keep links visible while panning.
        if (this.browser?.dragObject || this.browser?.isTrackPanning?.()) {
            return;
        }
        // Clear stale inter-frame links immediately while new locus data is loading.
        this.clearInterOverlay();
        this.delegateTrack.trackView = {
            ...this.trackView,
            axisCanvas: {
                style: { display: "none" },
            },
        };
    }

    private splitReadColor(readName: string): string {
        // Deterministic per-read hue for split reads.
        let hash = 0;
        for (let i = 0; i < readName.length; i++) {
            hash = ((hash << 5) - hash + readName.charCodeAt(i)) | 0;
        }
        const hue = Math.abs(hash) % 360;
        return `hsl(${hue}, 70%, 45%)`;
    }

    private splitReadColorAlpha(readName: string, alpha: number): string {
        let hash = 0;
        for (let i = 0; i < readName.length; i++) {
            hash = ((hash << 5) - hash + readName.charCodeAt(i)) | 0;
        }
        const hue = Math.abs(hash) % 360;
        return `hsla(${hue}, 70%, 45%, ${alpha})`;
    }




    private requestHoverRepaint(viewport?: any) {
         viewport?.repaint?.();
    }

    private pickSplitReadNameFromHits(hits: any[]): string | null {
        const first = Array.isArray(hits) && hits.length > 0 ? hits[0] : null;
        const candidate =
            first?.firstAlignment && this.hasSplitTag(first.firstAlignment) ? first.firstAlignment :
            first?.secondAlignment && this.hasSplitTag(first.secondAlignment) ? first.secondAlignment :
            this.hasSplitTag(first) ? first :
            null;
        return typeof candidate?.readName === "string" ? candidate.readName : null;
    }

    private buildProxyClickState(sourceViewport: any, event: MouseEvent): any {
        const delegateViewports = this.delegateTrack?.trackView?.viewports || [];
        let delegateViewport = delegateViewports.find((vp: any) => vp?.referenceFrame === sourceViewport?.referenceFrame);
        if (!delegateViewport && Array.isArray(this.trackView?.viewports)) {
            const idx = this.trackView.viewports.indexOf(sourceViewport);
            if (idx >= 0 && idx < delegateViewports.length) {
                delegateViewport = delegateViewports[idx];
            }
        }
        if (!delegateViewport || typeof delegateViewport.createClickState !== "function") return null;
        return delegateViewport.createClickState(event);
    }

    private markInteractionDirty() {
        this._interactionVersion += 1;
    }



    private isInteractionSynced() {
        return this._syncedInteractionVersion === this._interactionVersion;
    }

    private ensureHoverHandlers() {
        const viewports = this.trackView?.viewports || [];
        for (const vp of viewports) {
            if (!vp?.viewportElement || this._hoverBoundViewports.has(vp)) continue;
            this._hoverBoundViewports.add(vp);
            const evaluateHover = (event: MouseEvent) => {
                if (this.browser?.dragObject || this.browser?.isScrolling) return;
                const clickState = this.buildProxyClickState(vp, event)
                    || (typeof vp.createClickState === "function" ? vp.createClickState(event) : null);
                if (!clickState) return;
                // Normalize packed group bounds in-place for this exact hover state.
                // This fixes slight-pan cases where group.pixelBottom remains 0.
                //this.gatherPackedAlignments(clickState?.viewport?.cachedFeatures);
                const hits = this.delegateTrack?.clickedFeatures?.(clickState) || [];
                const readName = this.pickSplitReadNameFromHits(hits);
                if (readName !== this._activeReadName) {
                    this._activeReadName = readName;
                    this.requestHoverRepaint(vp);
                }
            };
            const move = (event: MouseEvent) => {
                if (this._activeReadName === null) {
                    this._lastHoverEvalTs = Date.now();
                    evaluateHover(event);
                    return;
                }
                const now = Date.now();
                if (now - this._lastHoverEvalTs < 90) return;
                if (this._hoverEvalPending) return;
                this._hoverEvalPending = true;
                const run = () => {
                    this._hoverEvalPending = false;
                    this._lastHoverEvalTs = Date.now();
                    evaluateHover(event);
                };
                if (typeof window !== "undefined" && typeof window.requestAnimationFrame === "function") {
                    window.requestAnimationFrame(run);
                } else {
                    run();
                }
            };
        
            const onMouseUp = (event: MouseEvent) => {
                if (this.isInteractionSynced()) evaluateHover(event);
            };
            const onMouseDown = (event: MouseEvent) => {
                this.markInteractionDirty();
            };
            const leave = () => {
                if (this._activeReadName !== null) {
                    this._activeReadName = null;
                    this.requestHoverRepaint(vp);
                }
            };
            const enter = (event: MouseEvent) => {
                evaluateHover(event);
            };
            const syncScroll = () => {
                this.scheduleOverlayRedrawFromScroll();
            };
            vp.viewportElement.addEventListener("mousedown", onMouseDown);
            vp.viewportElement.addEventListener("mouseup", onMouseUp);
            vp.viewportElement.addEventListener("mouseenter", enter);
            vp.viewportElement.addEventListener("mousemove", move);
            vp.viewportElement.addEventListener("mouseleave", leave);
            vp.viewportElement.addEventListener("scroll", syncScroll, { passive: true });
            vp.viewportElement.addEventListener("wheel", syncScroll, { passive: true });
        }
    }

    private schedulePostFetchRepaint() {
        if (this._postFetchRepaintPending) return;
        this._postFetchRepaintPending = true;
        const run = () => {
            this._postFetchRepaintPending = false;
            this.trackView?.repaintViews?.();
        };
        if (typeof window !== "undefined" && typeof window.requestAnimationFrame === "function") {
            window.requestAnimationFrame(run);
        } else {
            setTimeout(run, 0);
        }
    }

    async getFeatures(chr: string, start: number, end: number, bpPerPixel: number) {
        const features = await this.delegateTrack.getFeatures(chr, start, end, bpPerPixel);
        // IGV can repack on short pans during getFeatures; force our overlay pass to resync.
        this.markInteractionDirty();
        this.schedulePostFetchRepaint();
        return features;
    }

    clickedFeatures(clickState: any) {
        return this.delegateTrack?.clickedFeatures?.(clickState) || [];
    }

    async popupData(clickState: any) {
        const delegate = this.delegateTrack;
        if (!delegate) return [];
        return await delegate.popupData?.(clickState) || [];
    }
 

    private getVisibleReferenceFrames() {
        const viewports = this.delegateTrack?.trackView?.viewports || [];
        const frames: Array<{ chr: string; start: number; end: number; bpPerPixel: number; viewport: any }> = [];
        for (const vp of viewports) {
            if (!vp?.isVisible?.()) continue;
            const rf = vp.referenceFrame;
            const chr = rf?.chr;
            const start = Number(rf?.start);
            const width = Number(vp?.getWidth?.());
            const bpp = Number(rf?.bpPerPixel);
            if (typeof chr !== "string" || !Number.isFinite(start) || !Number.isFinite(width) || !Number.isFinite(bpp)) {
                continue;
            }
            frames.push({
                chr,
                start,
                end: start + width * bpp,
                bpPerPixel: bpp,
                viewport: vp,
            });
        }
        return frames;
    }

    private ensureInterOverlayCanvas(): HTMLCanvasElement | null {
        if (typeof document === "undefined") return null;
        const host: HTMLElement | null = this.browser?.columnContainer || this.browser?.root || null;
        if (!host) return null;
        if (getComputedStyle(host).position === "static") {
            host.style.position = "relative";
        }
        if (!this._interOverlay || !host.contains(this._interOverlay)) {
            const existing = host.querySelector("canvas[data-mdv-split-overlay='inter']") as HTMLCanvasElement | null;
            if (existing) {
                this._interOverlay = existing;
            } else {
                const c = document.createElement("canvas");
                c.dataset.mdvSplitOverlay = "inter";
                c.style.position = "absolute";
                c.style.left = "0";
                c.style.top = "0";
                c.style.width = "100%";
                c.style.height = "100%";
                c.style.pointerEvents = "none";
                c.style.zIndex = "9999";
                host.appendChild(c);
                this._interOverlay = c;
            }
        }
        const rect = host.getBoundingClientRect();
        const w = Math.max(1, Math.round(rect.width));
        const h = Math.max(1, Math.round(rect.height));
        if (this._interOverlay.width !== w || this._interOverlay.height !== h) {
            this._interOverlay.width = w;
            this._interOverlay.height = h;
        }
        return this._interOverlay;
    }

    private clearInterOverlay() {
        const host: HTMLElement | null = this.browser?.columnContainer || this.browser?.root || null;
        if (!host) return;
        const overlays = host.querySelectorAll("canvas[data-mdv-split-overlay='inter']");
        overlays.forEach((c) => {
            const canvas = c as HTMLCanvasElement;
            const ctx = canvas.getContext("2d");
            if (ctx) ctx.clearRect(0, 0, canvas.width, canvas.height);
            canvas.remove();
        });
        this._interOverlay = null;
    }

    clearSplitOverlay() {
        this.clearInterOverlay();
    }

    clearSplitReadFilter() {
        this._activeReadName = null;
        if (this.delegateTrack?.alignmentTrack) {
            this.delegateTrack.alignmentTrack.selectedReadName = undefined;
        }
        this.trackView?.repaintViews?.();
    }






    private parseAlignmentReadInterval(alignment: any): {
        readStart: number;
        readEnd: number;
        strandPlus: boolean;
        clipAtReadStart: boolean;
        clipAtReadEnd: boolean;
    } | null {
        const strandPlus = alignment?.strand !== false;
        const interval = computeOriginalReadIntervalFromCigar(alignment?.cigar, strandPlus);
        if (!interval) return null;
        const { leftClip, rightClip } = getTerminalClipTotals(alignment?.cigar);
        const { leftMatch, rightMatch } = getTerminalEdgeMatchTotals(alignment?.cigar);
        const leftClipEffective = leftClip >= MIN_ROLE_CLIP_BP && leftClip > leftMatch;
        const rightClipEffective = rightClip >= MIN_ROLE_CLIP_BP && rightClip > rightMatch;
        const clipAtReadStart = strandPlus ? leftClipEffective : rightClipEffective;
        const clipAtReadEnd = strandPlus ? rightClipEffective : leftClipEffective;
        return {
            readStart: interval.readStart,
            readEnd: interval.readEnd,
            strandPlus,
            clipAtReadStart,
            clipAtReadEnd,
        };
    }

    private hasSplitTag(alignment: any): boolean {
        const tags = typeof alignment?.tags === "function" ? alignment.tags() : undefined;
        return Boolean(tags?.SA ?? tags?.["SA"]);
    }

    private orderNodesByReadPath<T extends {
        readStart: number;
        readEnd: number;
        clipAtReadStart?: boolean;
        clipAtReadEnd?: boolean;
    }>(nodes: T[]): T[] {
        // Build a path using explicit read-coordinate continuity:
        // choose the next segment whose readStart is closest to previous readEnd.
        const remaining = [...nodes].sort((a, b) => (a.readStart - b.readStart) || (a.readEnd - b.readEnd));
        if (remaining.length <= 1) return remaining;

        // If present, seed with segment clipped at read end (first segment),
        // else earliest by readStart.
        let seedIdx = remaining.findIndex(n => Boolean(n.clipAtReadEnd) && !Boolean(n.clipAtReadStart));
        if (seedIdx < 0) {
            seedIdx = remaining.findIndex(n => Boolean(n.clipAtReadEnd));
        }
        if (seedIdx < 0) seedIdx = 0;
        const ordered: T[] = [remaining.splice(seedIdx, 1)[0] as T];
        while (remaining.length > 0) {
            const prev = ordered[ordered.length - 1];
            let bestIdx = 0;
            let bestDist = Number.POSITIVE_INFINITY;
            for (let i = 0; i < remaining.length; i++) {
                const cand = remaining[i];
                const dist = Math.abs(cand.readStart - prev.readEnd);
                if (dist < bestDist) {
                    bestDist = dist;
                    bestIdx = i;
                } else if (dist === bestDist) {
                    const cur = remaining[bestIdx];
                    if (cand.readStart < cur.readStart || (cand.readStart === cur.readStart && cand.readEnd < cur.readEnd)) {
                        bestIdx = i;
                    }
                }
            }
            ordered.push(remaining.splice(bestIdx, 1)[0] as T);
        }
        return ordered;
    }

    private gatherPackedAlignments(alignmentContainer: any) {
        const packedGroups = alignmentContainer?.packedGroups;
        const rows: Array<{ y: number; alignments: any[] }> = [];
        if (!packedGroups || typeof packedGroups.values !== "function") {
            return rows;
        }

        const alignmentRowHeight = this.delegateTrack?.displayMode === "SQUISHED"
            ? this.delegateTrack?.squishedRowHeight || 8
            : this.delegateTrack?.alignmentRowHeight || 14;
        const yInset = Number(this.delegateTrack?.alignmentsYOffset) || 0;
        let y = yInset;
        for (const group of packedGroups.values()) {
            const groupRows = group?.rows || [];
            for (const row of groupRows) {
                const rowTop = Number(row?.pixelTop);
                const rowCenter = Number.isFinite(rowTop)
                    ? rowTop + alignmentRowHeight / 2
                    : y + alignmentRowHeight / 2;
                rows.push({
                    y: rowCenter,
                    alignments: row?.alignments || [],
                });
                y += alignmentRowHeight;
            }
            y += this.delegateTrack?.groupBy ? 5 : 0;
        }
        return rows;
    }

    private buildRowsByViewport(
        visibleFrames: Array<{ chr: string; start: number; end: number; bpPerPixel: number; viewport: any }>,
    ) {
        const rowsByViewport = new Map<any, Array<{ y: number; alignments: any[] }>>();
        for (const frame of visibleFrames) {
            const cached = frame.viewport?.cachedFeatures;
            rowsByViewport.set(frame.viewport, this.gatherPackedAlignments(cached));
        }
        return rowsByViewport;
    }

    private redrawInterFrameOverlayOnly() {
        const visibleFrames = this.getVisibleReferenceFrames();
        if (visibleFrames.length < 2) {
            this.clearInterOverlay();
            return;
        }
        const rowsByViewport = this.buildRowsByViewport(visibleFrames);
        this.drawInterFrameOverlay(visibleFrames, rowsByViewport);
    }

    private drawInterFrameOverlay(
        visibleFrames: Array<{ chr: string; start: number; end: number; bpPerPixel: number; viewport: any }>,
        rowsByViewport: Map<any, Array<{ y: number; alignments: any[] }>>,
    ) {
        if (visibleFrames.length < 2) {
            this.clearInterOverlay();
            return;
        }
        const overlay = this.ensureInterOverlayCanvas();
        const octx = overlay?.getContext("2d");
        if (!overlay || !octx) return;
        octx.clearRect(0, 0, overlay.width, overlay.height);
        const host: HTMLElement | null = this.browser?.columnContainer || this.browser?.root || null;
        const hostRect = host?.getBoundingClientRect?.();
        if (!hostRect) return;

        type Node = {
            readName: string;
            frame: { chr: string; start: number; end: number; bpPerPixel: number; viewport: any };
            y: number;
            aStart: number;
            aEnd: number;
            readStart: number;
            readEnd: number;
            strandPlus: boolean;
            clipAtReadStart: boolean;
            clipAtReadEnd: boolean;
            hasSplit: boolean;
        };

        const nodesByRead = new Map<string, Node[]>();
        for (const frame of visibleFrames) {
            const rows = rowsByViewport.get(frame.viewport) || [];
            for (const row of rows) {
                for (const item of row.alignments) {
                    const alignments = item?.paired ? [item.firstAlignment, item.secondAlignment].filter(Boolean) : [item];
                    for (const alignment of alignments) {
                        const readName = typeof alignment?.readName === "string" ? alignment.readName : "";
                        if (!readName) continue;
                        const aStart = Number(alignment.start);
                        const aEnd = aStart + Number(alignment.lengthOnRef || 0);
                        if (!Number.isFinite(aStart) || !Number.isFinite(aEnd)) continue;
                        const interval = this.parseAlignmentReadInterval(alignment);
                        if (!interval) continue;
                        const list = nodesByRead.get(readName) || [];
                        list.push({
                            readName,
                            frame,
                            y: row.y,
                            aStart,
                            aEnd,
                            readStart: interval.readStart,
                            readEnd: interval.readEnd,
                            strandPlus: interval.strandPlus,
                            clipAtReadStart: interval.clipAtReadStart,
                            clipAtReadEnd: interval.clipAtReadEnd,
                            hasSplit: this.hasSplitTag(alignment),
                        });
                        nodesByRead.set(readName, list);
                    }
                }
            }
        }

        octx.save();
        octx.lineWidth = this.config.interSplitLineWidth || 1.8;
        for (const [readName, nodes] of nodesByRead.entries()) {
            if (this._activeReadName && readName !== this._activeReadName) continue;
            if (nodes.length < 2) continue;
            if (!nodes.some(n => n.hasSplit)) continue;
            const orderedNodes = this.orderNodesByReadPath(nodes);
            octx.strokeStyle = this.splitReadColor(readName);
            for (let i = 0; i < orderedNodes.length - 1; i++) {
                const source = orderedNodes[i];
                const target = orderedNodes[i + 1];
                if (source.frame.viewport === target.frame.viewport) continue;

                const sourceRect = source.frame.viewport?.viewportElement?.getBoundingClientRect?.();
                const targetRect = target.frame.viewport?.viewportElement?.getBoundingClientRect?.();
                if (!sourceRect || !targetRect) continue;

                const sourceBp = endpointCoord(source.aStart, source.aEnd, source.strandPlus, "end");
                const targetBp = endpointCoord(target.aStart, target.aEnd, target.strandPlus, "start");

                const sourceContentTop = Number(source.frame.viewport?.getContentTop?.() || 0);
                const targetContentTop = Number(target.frame.viewport?.getContentTop?.() || 0);
                const x1 = (sourceRect.left - hostRect.left) + ((sourceBp - source.frame.start) / source.frame.bpPerPixel);
                const y1 = (sourceRect.top - hostRect.top) + (source.y - sourceContentTop);
                const x2 = (targetRect.left - hostRect.left) + ((targetBp - target.frame.start) / target.frame.bpPerPixel);
                const y2 = (targetRect.top - hostRect.top) + (target.y - targetContentTop);

                if (x1 < (sourceRect.left - hostRect.left) || x1 > (sourceRect.right - hostRect.left)) continue;
                if (x2 < (targetRect.left - hostRect.left) || x2 > (targetRect.right - hostRect.left)) continue;

                const cx = (x1 + x2) / 2;
                const cy = Math.min(y1, y2) - Math.max(16, Math.min(72, Math.abs(x2 - x1) * 0.12));
                const sourceLeft = sourceRect.left - hostRect.left;
                const sourceTop = sourceRect.top - hostRect.top;
                const sourceWidth = Math.max(1, sourceRect.width);
                const sourceHeight = Math.max(1, sourceRect.height);
                const targetLeft = targetRect.left - hostRect.left;
                const targetTop = targetRect.top - hostRect.top;
                const targetWidth = Math.max(1, targetRect.width);
                const targetHeight = Math.max(1, targetRect.height);

                octx.save();
                octx.beginPath();
                octx.rect(sourceLeft, sourceTop, sourceWidth, sourceHeight);
                octx.rect(targetLeft, targetTop, targetWidth, targetHeight);
                octx.clip();
                octx.beginPath();
                octx.moveTo(x1, y1);
                octx.quadraticCurveTo(cx, cy, x2, y2);
                octx.stroke();
                octx.restore();
            }
        }
        octx.restore();
    }

    private drawSplitOverlays(options: any) {
        this.ensureHoverHandlers();
        const ctx: CanvasRenderingContext2D | undefined = options?.context;
        const referenceFrame = options?.referenceFrame;
        const bpStart = Number(options?.bpStart);
        const bpPerPixel = Number(options?.bpPerPixel);
        const pixelWidth = Number(options?.pixelWidth);
        const alignmentContainer = options?.features;
        if (!ctx || !referenceFrame || !Number.isFinite(bpStart) || !Number.isFinite(bpPerPixel) || !Number.isFinite(pixelWidth)) {
            return;
        }

        const visibleFrames = this.getVisibleReferenceFrames();
        if (visibleFrames.length < 2) {
            this.clearInterOverlay();
        }
        const rowsByViewport = this.buildRowsByViewport(visibleFrames);

        let rows = rowsByViewport.get(options?.viewport) || [];
        if (rows.length === 0) {
            rows = this.gatherPackedAlignments(alignmentContainer);
        }
        type LocalNode = {
            readName: string;
            y: number;
            aStart: number;
            aEnd: number;
            readStart: number;
            readEnd: number;
            strandPlus: boolean;
            clipAtReadStart: boolean;
            clipAtReadEnd: boolean;
            hasSplit: boolean;
        };
        const rowEntries: Array<{ y: number; alignment: any }> = [];
        const nodesByRead = new Map<string, LocalNode[]>();
        if (rows.length > 0) {
            for (const row of rows) {
                for (const a0 of row.alignments) {
                    const alignments = a0?.paired ? [a0.firstAlignment, a0.secondAlignment].filter(Boolean) : [a0];
                    for (const alignment of alignments) {
                        rowEntries.push({ y: row.y, alignment });
                        const readName = typeof alignment?.readName === "string" ? alignment.readName : "";
                        if (!readName) continue;
                        const aStart = Number(alignment.start);
                        const aEnd = aStart + Number(alignment.lengthOnRef || 0);
                        if (!Number.isFinite(aStart) || !Number.isFinite(aEnd)) continue;
                        if (aEnd < bpStart || aStart > bpStart + pixelWidth * bpPerPixel) continue;
                        const interval = this.parseAlignmentReadInterval(alignment);
                        if (!interval) continue;
                        const list = nodesByRead.get(readName) || [];
                        list.push({
                            readName,
                            y: row.y,
                            aStart,
                            aEnd,
                            readStart: interval.readStart,
                            readEnd: interval.readEnd,
                            strandPlus: interval.strandPlus,
                            clipAtReadStart: interval.clipAtReadStart,
                            clipAtReadEnd: interval.clipAtReadEnd,
                            hasSplit: this.hasSplitTag(alignment),
                        });
                        nodesByRead.set(readName, list);
                    }
                }
            }
        }

        ctx.save();
        ctx.lineWidth = this.config.splitLineWidth || 1.5;
        for (const [readName, nodes] of nodesByRead.entries()) {
            if (this._activeReadName && readName !== this._activeReadName) continue;
            if (nodes.length < 2) continue;
            if (!nodes.some(n => n.hasSplit)) continue;
            const orderedNodes = this.orderNodesByReadPath(nodes);
            ctx.strokeStyle = this.splitReadColor(readName);
            for (let i = 0; i < orderedNodes.length - 1; i++) {
                const source = orderedNodes[i];
                const target = orderedNodes[i + 1];
                const sourceBp = endpointCoord(source.aStart, source.aEnd, source.strandPlus, "end");
                const targetBp = endpointCoord(target.aStart, target.aEnd, target.strandPlus, "start");
                const sourceX = (sourceBp - bpStart) / bpPerPixel;
                const targetX = (targetBp - bpStart) / bpPerPixel;
                if (sourceX < 0 || sourceX > pixelWidth) continue;
                if (targetX < 0 || targetX > pixelWidth) continue;
                const midX = (sourceX + targetX) / 2;
                const midY = Math.min(source.y, target.y) - Math.max(8, Math.min(32, Math.abs(targetX - sourceX) * 0.14));
                ctx.beginPath();
                ctx.moveTo(sourceX, source.y);
                ctx.quadraticCurveTo(midX, midY, targetX, target.y);
                ctx.stroke();
            }
        }

        for (const { y, alignment } of rowEntries) {
            const readName = typeof alignment?.readName === "string" ? alignment.readName : "";
            if (!readName || !this.hasSplitTag(alignment)) continue;
            const aStart = Number(alignment.start);
            const aEnd = aStart + Number(alignment.lengthOnRef || 0);
            if (!Number.isFinite(aStart) || !Number.isFinite(aEnd)) continue;
            const left = Math.max(0, (aStart - bpStart) / bpPerPixel);
            const right = Math.min(pixelWidth, (aEnd - bpStart) / bpPerPixel);
            const width = right - left;
            if (width <= 1) continue;
            const rowHeight = this.delegateTrack?.displayMode === "SQUISHED"
                ? this.delegateTrack?.squishedRowHeight || 8
                : this.delegateTrack?.alignmentRowHeight || 14;
            const outlineTop = y - rowHeight / 2 + 0.5;
            const outlineHeight = Math.max(4, rowHeight - 2);
            ctx.save();
            const isHovered = this._activeReadName === readName;
            const hasHover = this._activeReadName !== null;
            if (hasHover && !isHovered) {
                ctx.fillStyle = "rgba(150, 150, 150, 0.16)";
                ctx.strokeStyle = "rgba(120, 120, 120, 0.9)";
            } else {
                ctx.fillStyle = this.splitReadColorAlpha(readName, 0.20);
                ctx.strokeStyle = this.splitReadColor(readName);
            }
            ctx.fillRect(left, outlineTop, width, outlineHeight);
            ctx.lineWidth = 1;
            ctx.strokeRect(left, outlineTop, width, outlineHeight);
            ctx.restore();
        }
        ctx.restore();
        this.drawInterFrameOverlay(visibleFrames, rowsByViewport);
    }

    computePixelHeight(features: any) {
        if (typeof this.delegateTrack?.computePixelHeight === "function") {
            return this.delegateTrack.computePixelHeight(features);
        }
        return this.config.height || 220;
    }

    draw(options: any) {
        this.delegateTrack.draw(options);
        this.drawSplitOverlays(options);
    }
}

igv.registerTrackClass("mdv_split_alignment_track", MdvSplitAlignmentTrack);
