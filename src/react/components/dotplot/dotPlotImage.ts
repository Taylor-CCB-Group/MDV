const SVG_NS = "http://www.w3.org/2000/svg";

/**
 * Computed style properties that decide how the chart looks. Inlining them makes the image
 * standalone: CSS variables, Tailwind classes and the page theme don't exist outside the page.
 */
const SVG_STYLE_PROPERTIES = [
    "fill",
    "fill-opacity",
    "stroke",
    "stroke-width",
    "stroke-opacity",
    "stroke-dasharray",
    "opacity",
    "visibility",
    "overflow",
    "display",
    "font-family",
    "font-size",
    "font-weight",
    "font-style",
    "text-anchor",
    "dominant-baseline",
    "alignment-baseline",
    "stop-color",
    "stop-opacity",
];

/** Selector for legend titles drawn as HTML rather than SVG (e.g. the fraction legend's). */
const HTML_TITLE_SELECTOR = ".legend-title";

/** Replaces each element's class and style with its computed style, recursing in step. */
function inlineStyles(copy: Element, source: Element) {
    const style = getComputedStyle(source);
    const declarations = SVG_STYLE_PROPERTIES.flatMap((property) => {
        const value = style.getPropertyValue(property);
        return value ? [`${property}:${value}`] : [];
    });
    copy.setAttribute("style", declarations.join(";"));
    copy.removeAttribute("class");
    for (let i = 0; i < source.children.length; i++) {
        const child = copy.children[i];
        if (child) inlineStyles(child, source.children[i]);
    }
}

/**
 * One standalone SVG of everything drawn in `container`: each SVG in it (plot and legends)
 * nested at its place on screen, and HTML legend titles as text. Styles are inlined, so it
 * looks the same outside the page.
 */
export function buildChartSvg(container: HTMLElement): SVGSVGElement {
    const box = container.getBoundingClientRect();
    const doc = container.ownerDocument;
    const out = doc.createElementNS(SVG_NS, "svg");
    out.setAttribute("xmlns", SVG_NS);
    out.setAttribute("width", String(box.width));
    out.setAttribute("height", String(box.height));
    out.setAttribute("viewBox", `0 0 ${box.width} ${box.height}`);

    for (const svg of container.querySelectorAll("svg")) {
        // nested SVGs come along with their outer one
        if (svg.parentElement?.closest("svg")) continue;
        const rect = svg.getBoundingClientRect();
        if (rect.width === 0 || rect.height === 0) continue;
        const copy = svg.cloneNode(true);
        if (!(copy instanceof Element)) continue;
        inlineStyles(copy, svg);
        copy.setAttribute("x", String(rect.left - box.left));
        copy.setAttribute("y", String(rect.top - box.top));
        copy.setAttribute("width", String(rect.width));
        copy.setAttribute("height", String(rect.height));
        out.append(copy);
    }

    for (const title of container.querySelectorAll(HTML_TITLE_SELECTOR)) {
        const text = title.textContent?.trim();
        if (!text || title.closest("svg")) continue;
        const rect = title.getBoundingClientRect();
        if (rect.width === 0 || rect.height === 0) continue;
        const style = getComputedStyle(title);
        const label = doc.createElementNS(SVG_NS, "text");
        label.setAttribute("x", String(rect.left - box.left));
        label.setAttribute("y", String(rect.top - box.top + rect.height / 2));
        label.setAttribute(
            "style",
            [
                `fill:${style.color}`,
                `font-family:${style.fontFamily}`,
                `font-size:${style.fontSize}`,
                `font-weight:${style.fontWeight}`,
                "dominant-baseline:middle",
            ].join(";"),
        );
        label.textContent = text;
        out.append(label);
    }
    return out;
}
