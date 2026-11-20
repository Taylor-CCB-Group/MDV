import type { Viewport } from "@deck.gl/core";

/** The parts of deck's `TextLayer` read here (structurally, as more than one copy of its type can be installed). */
type TextLayerLike = {
    props: {
        data?: unknown;
        getPosition?: unknown;
        getSize?: unknown;
        sizeScale?: number;
        backgroundPadding?: readonly number[];
        getPixelOffset?: unknown;
    };
    state?: unknown;
};

type Box = [left: number, top: number, right: number, bottom: number];

/** `backgroundPadding` as [left, top, right, bottom]. */
function paddingOf(padding: readonly number[] | undefined): Box {
    if (!padding || padding.length < 2) return [0, 0, 0, 0];
    if (padding.length < 4) return [padding[0], padding[1], padding[0], padding[1]];
    return [padding[0], padding[1], padding[2], padding[3]];
}

/** A prop that may be an accessor, for one object. */
function accessorValue(prop: unknown, object: unknown, index: number): unknown {
    return typeof prop === "function" ? prop(object, { index }) : prop;
}

/** `object[key]`, if `object` is an object. */
function field(object: unknown, key: string): unknown {
    return object && typeof object === "object" ? Reflect.get(object, key) : undefined;
}

const isNumbers = (value: unknown, length: number): value is number[] =>
    Array.isArray(value) && value.length >= length && value.every((n) => typeof n === "number");

/**
 * The index of the text label whose background box contains `point` (pixels in `viewport`), or -1;
 * null if the layer can't tell yet (its font atlas isn't laid out).
 *
 * Worked out on the CPU, as deck draws the background (`TextBackgroundLayer`): the layer's own
 * text box for each label (`getBoundingRect`, in font-atlas units, a private accessor), scaled by
 * the label's pixel size over the atlas font size, plus the padding. So a press can be matched to a
 * label straight away, on WebGPU too, where deck's own pick of the press comes back later.
 * Assumes pixel `sizeUnits` and no rotation, as MDV's labels are.
 */
export function textLabelAt(layer: TextLayerLike, [x, y]: readonly number[], viewport: Viewport): number | null {
    const getBoundingRect: unknown = Reflect.get(layer, "getBoundingRect");
    const fontSize = field(field(field(layer.state, "fontAtlasManager"), "props"), "fontSize");
    const data: unknown = layer.props.data;
    if (typeof getBoundingRect !== "function" || typeof fontSize !== "number" || !Array.isArray(data)) return null;
    const { getPosition, getSize, sizeScale = 1, backgroundPadding, getPixelOffset } = layer.props;
    const [padL, padT, padR, padB] = paddingOf(backgroundPadding);
    for (let index = data.length - 1; index >= 0; index--) {
        const object: unknown = data[index];
        const position = accessorValue(getPosition, object, index);
        const rect: unknown = getBoundingRect.call(layer, object, { index, data, target: [] });
        if (!isNumbers(position, 2) || !isNumbers(rect, 4)) continue;
        const [cx, cy] = viewport.project(position);
        const offset = accessorValue(getPixelOffset, object, index);
        const [ox, oy] = isNumbers(offset, 2) ? offset : [0, 0];
        const size = accessorValue(getSize, object, index);
        const scale = ((typeof size === "number" ? size : 32) * sizeScale) / fontSize;
        const [rx, ry, rw, rh] = rect;
        const left = cx + ox + rx * scale - padL;
        const top = cy + oy + ry * scale - padT;
        if (x >= left && x <= left + rw * scale + padL + padR && y >= top && y <= top + rh * scale + padT + padB) {
            return index;
        }
    }
    return -1;
}
