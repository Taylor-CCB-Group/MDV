/**
 * Puts the drawings of a chart's legends into its SVG for an exported image, where each is shown
 * on screen, and returns a function that puts everything back afterwards.
 *
 * Each drawing is placed by where its own SVG sits relative to the chart's, so this works for
 * legends floating over the plot and for those in the strip above it. When legends sit above the
 * plot, the image is made taller by that much and they're drawn above the plot, as on screen.
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
        const box = (drawing.parentElement?.closest("svg") ?? legend).getBoundingClientRect();
        const left = Math.round(box.left - origin.left);
        const top = Math.round(box.top - origin.top);
        above = Math.max(above, -top);
        const parent = drawing.parentNode;
        const next = drawing.nextSibling;
        const transform = drawing.getAttribute("transform");
        // the drawing's own offset inside its SVG still applies
        const offset = `translate(${left},${top})`;
        drawing.setAttribute("transform", transform ? `${offset} ${transform}` : offset);
        svg.append(drawing);
        restorers.push(() => {
            setOrRemove(drawing, "transform", transform);
            parent?.insertBefore(drawing, next);
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
