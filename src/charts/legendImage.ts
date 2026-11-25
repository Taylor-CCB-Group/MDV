const SVG_NS = "http://www.w3.org/2000/svg";

/**
 * Where `element` can be seen on screen: its box, cut to each element around it, up to and
 * including `outermost`, that hides what overflows it (e.g. a long legend list that scrolls).
 */
export function getVisibleBox(element: Element, outermost: Element): DOMRect {
    const box = element.getBoundingClientRect();
    let { left, top, right, bottom } = box;
    const win = element.ownerDocument.defaultView;
    for (let el = element.parentElement; el && win; el = el.parentElement) {
        const style = win.getComputedStyle(el);
        const clipsX = clips(style.overflowX || style.overflow);
        const clipsY = clips(style.overflowY || style.overflow);
        if (clipsX || clipsY) {
            const clip = el.getBoundingClientRect();
            if (clipsX) {
                left = Math.max(left, clip.left);
                right = Math.min(right, clip.right);
            }
            if (clipsY) {
                top = Math.max(top, clip.top);
                bottom = Math.min(bottom, clip.bottom);
            }
        }
        if (el === outermost) break;
    }
    return DOMRect.fromRect({ x: left, y: top, width: Math.max(0, right - left), height: Math.max(0, bottom - top) });
}

/** Whether an `overflow` value hides what doesn't fit. */
function clips(overflow: string): boolean {
    return ["hidden", "clip", "auto", "scroll", "overlay"].includes(overflow);
}

/** The outermost element around `legend` that doesn't also hold the chart's SVG. */
function outermostAround(legend: Element, svg: SVGSVGElement): Element {
    let outermost = legend;
    while (outermost.parentElement && !outermost.parentElement.contains(svg)) {
        outermost = outermost.parentElement;
    }
    return outermost;
}

/**
 * Puts the drawings of a chart's legends into its SVG for an exported image, where each is shown
 * on screen, and returns a function that puts everything back afterwards.
 *
 * Each drawing is placed by where its own SVG sits relative to the chart's, so this works for
 * legends floating over the plot and for those in the strip above it, and is cut to the part
 * that can be seen (a long list scrolled in its legend shows as it does on screen). When legends
 * sit above the plot, the image is made taller by that much and they're drawn above the plot.
 *
 * `legends` are the elements charts keep for their legends: a `LegendWrapper`'s host, or an older
 * legend. A minimised legend has no drawing, so it's left out.
 */
export function addLegendsToSvg(svg: SVGSVGElement, legends: (Element | null | undefined)[]): () => void {
    const restorers: (() => void)[] = [];
    const origin = svg.getBoundingClientRect();
    let above = 0;
    for (const legend of legends) {
        if (!legend) continue;
        const drawing = legend.querySelector("g");
        if (!drawing) continue;
        const own = drawing.parentElement?.closest("svg") ?? legend;
        const box = own.getBoundingClientRect();
        const visible = getVisibleBox(own, outermostAround(legend, svg));
        if (visible.width === 0 || visible.height === 0) continue;
        const left = Math.round(visible.left - origin.left);
        const top = Math.round(visible.top - origin.top);
        above = Math.max(above, -top);
        // a nested SVG the size of the visible part, which cuts off the rest
        const frame = svg.ownerDocument.createElementNS(SVG_NS, "svg");
        frame.setAttribute("x", String(left));
        frame.setAttribute("y", String(top));
        frame.setAttribute("width", String(Math.round(visible.width)));
        frame.setAttribute("height", String(Math.round(visible.height)));
        const parent = drawing.parentNode;
        const next = drawing.nextSibling;
        const transform = drawing.getAttribute("transform");
        // the drawing's own offset inside its SVG still applies
        const offset = `translate(${Math.round(box.left - visible.left)},${Math.round(box.top - visible.top)})`;
        drawing.setAttribute("transform", transform ? `${offset} ${transform}` : offset);
        frame.append(drawing);
        svg.append(frame);
        restorers.push(() => {
            setOrRemove(drawing, "transform", transform);
            parent?.insertBefore(drawing, next);
            frame.remove();
        });
    }
    if (above > 0) restorers.push(extendUpwards(svg, above, origin));
    return () => {
        for (const restore of restorers) restore();
    };
}

/** Makes room above the chart's drawing by growing the SVG; returns a function that undoes it. */
function extendUpwards(svg: SVGSVGElement, extra: number, box: DOMRect): () => void {
    const height = svg.getAttribute("height");
    const viewBox = svg.getAttribute("viewBox");
    const width = Number.parseFloat(svg.getAttribute("width") ?? "") || box.width;
    const fullHeight = (Number.parseFloat(height ?? "") || box.height) + extra;
    svg.setAttribute("height", String(fullHeight));
    svg.setAttribute("viewBox", `0 ${-extra} ${width} ${fullHeight}`);
    return () => {
        setOrRemove(svg, "height", height);
        setOrRemove(svg, "viewBox", viewBox);
    };
}

function setOrRemove(element: Element, name: string, value: string | null) {
    if (value === null) element.removeAttribute(name);
    else element.setAttribute(name, value);
}
