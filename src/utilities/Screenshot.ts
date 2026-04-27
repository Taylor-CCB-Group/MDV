import { toPng } from "html-to-image";

const CAPTURE_EXCLUDED_CLASSES = [
    "ciview-ctm-main",
    "ciview-tooltip",
    "tippy-box",
    "tippy-arrow",
    "tippy-backdrop",
    "tippy-content",
];

const rectsIntersect = (a: DOMRect, b: DOMRect) =>
    a.left < b.right && a.right > b.left && a.top < b.bottom && a.bottom > b.top;

const shouldIncludeCaptureNode = (node: HTMLElement | SVGElement, root: HTMLElement, rootRect: DOMRect) => {
    if (node === root) return true;

    if (node instanceof HTMLElement) {
        if (node.hidden) return false;
        if (node.getAttribute("aria-hidden") === "true") return false;
        if (CAPTURE_EXCLUDED_CLASSES.some((cls) => node.classList.contains(cls))) {
            return false;
        }
    }

    const computedStyle = window.getComputedStyle(node);
    if (computedStyle.display === "none" || computedStyle.visibility === "hidden" || computedStyle.opacity === "0") {
        return false;
    }

    // Only do geometry checks for subtree roots so we prune expensive branches early
    // without paying getBoundingClientRect() on every individual SVG shape.
    const shouldCheckBounds =
        node instanceof HTMLElement || (node instanceof SVGSVGElement && node.ownerSVGElement === null);
    if (!shouldCheckBounds) {
        return true;
    }

    const rect = node.getBoundingClientRect();
    const hasNoBox =
        rect.width === 0 &&
        rect.height === 0 &&
        computedStyle.position !== "fixed" &&
        computedStyle.position !== "sticky";
    if (hasNoBox) return false;

    return rectsIntersect(rect, rootRect);
};

/**
 * Creates a png dataUrl image from a given DOM element, with some filtering of
 * out-of-bounds elements and tooltips etc.
 */
export default async (root: HTMLDivElement) => {
    const bounds = root.getBoundingClientRect();
    const aspect = bounds.width / bounds.height;
    const dataUrl = await toPng(root, {
        canvasWidth: 250,
        canvasHeight: 250 / aspect,
        // pixelRatio: 1, //makes very blurry
        filter: (node) => {
            if (!(node instanceof HTMLElement)) {
                // ostensible 'never' according to the types, but that's not accurate
                return true;
            }
            return shouldIncludeCaptureNode(node, root, bounds);
        },
    });
    return dataUrl;
};
