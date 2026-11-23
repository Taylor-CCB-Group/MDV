/**
 * Puts the drawings of a chart's legends into its SVG for an exported image, where each sits over
 * the plot, and returns a function that puts them back afterwards.
 *
 * `legends` are the elements charts keep for their legends: a `LegendWrapper`'s host, whose frame
 * (`.legend-container`) is what sits over the plot, or an older legend that is positioned itself.
 * A minimised legend has no drawing, so it's left out.
 */
export function addLegendsToSvg(svg: SVGSVGElement, legends: (Element | null | undefined)[]): () => void {
    const restorers: (() => void)[] = [];
    for (const legend of legends) {
        if (!legend) continue;
        const drawing = legend.querySelector("g");
        if (!drawing) continue;
        const frame = legend.querySelector(".legend-container") ?? legend;
        const left = frame instanceof HTMLElement ? frame.offsetLeft : 0;
        const top = frame instanceof HTMLElement ? frame.offsetTop : 0;
        const parent = drawing.parentNode;
        const next = drawing.nextSibling;
        const transform = drawing.getAttribute("transform");
        drawing.setAttribute("transform", `translate(${left},${top})`);
        svg.append(drawing);
        restorers.push(() => {
            if (transform === null) drawing.removeAttribute("transform");
            else drawing.setAttribute("transform", transform);
            parent?.insertBefore(drawing, next);
        });
    }
    return () => {
        for (const restore of restorers) restore();
    };
}
