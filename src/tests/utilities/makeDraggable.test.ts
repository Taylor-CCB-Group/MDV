import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import { makeDraggable } from "@/utilities/Elements.js";

/** jsdom has no layout, so read positions and sizes from inline styles. */
function px(value: string) {
    return Number.parseFloat(value) || 0;
}

function rect(width: number, height: number): DOMRect {
    return { x: 0, y: 0, left: 0, top: 0, width, height, right: width, bottom: height, toJSON: () => ({}) };
}

/** A box of the given size at (left, top), inside a parent of the given size. */
function setup(
    parentSize: [number, number],
    box: { left: number; top: number; width: number; height: number },
    contain: HTMLElement | "topleft" | undefined,
) {
    const parent = document.createElement("div");
    const el = document.createElement("div");
    parent.append(el);
    document.body.append(parent);
    el.style.left = `${box.left}px`;
    el.style.top = `${box.top}px`;
    el.style.width = `${box.width}px`;
    el.style.height = `${box.height}px`;
    vi.spyOn(parent, "getBoundingClientRect").mockReturnValue(rect(...parentSize));
    vi.spyOn(el, "getBoundingClientRect").mockReturnValue(rect(box.width, box.height));
    makeDraggable(el, { contain, doc: document });
    return el;
}

/** Drags by (dx, dy) in one mouse move. */
function drag(el: HTMLElement, dx: number, dy: number) {
    el.dispatchEvent(new MouseEvent("mousedown", { bubbles: true, clientX: 0, clientY: 0 }));
    document.dispatchEvent(new MouseEvent("mousemove", { clientX: dx, clientY: dy }));
    document.dispatchEvent(new MouseEvent("mouseup"));
    return [el.style.left, el.style.top];
}

beforeEach(() => {
    vi.spyOn(HTMLElement.prototype, "offsetLeft", "get").mockImplementation(function (this: HTMLElement) {
        return px(this.style.left);
    });
    vi.spyOn(HTMLElement.prototype, "offsetTop", "get").mockImplementation(function (this: HTMLElement) {
        return px(this.style.top);
    });
});

afterEach(() => {
    document.body.innerHTML = "";
    vi.restoreAllMocks();
});

describe("makeDraggable with a container", () => {
    const parentSize: [number, number] = [300, 200];
    const box = { left: 150, top: 100, width: 100, height: 80 };

    test("moves freely inside the container", () => {
        const el = setup(parentSize, box, document.body);
        expect(drag(el, -50, -40)).toEqual(["100px", "60px"]);
    });

    test("stops at the bottom-right edge instead of ignoring the move", () => {
        const el = setup(parentSize, box, document.body);
        expect(drag(el, 500, 500)).toEqual(["200px", "120px"]);
    });

    test("stops at the top-left edge", () => {
        const el = setup(parentSize, box, document.body);
        expect(drag(el, -500, -500)).toEqual(["0px", "0px"]);
    });

    test("an element that starts partly outside can be dragged back in", () => {
        // e.g. a legend saved near the edge of a chart that has since shrunk
        const el = setup(parentSize, { ...box, left: 250, top: 150 }, document.body);
        expect(drag(el, -10, -10)).toEqual(["200px", "120px"]);
    });
});

describe("makeDraggable with 'topleft' containment", () => {
    test("only stops at the top and left edges", () => {
        const el = setup([300, 200], { left: 150, top: 100, width: 100, height: 80 }, "topleft");
        expect(drag(el, 500, 500)).toEqual(["650px", "600px"]);
        expect(drag(el, -1000, -1000)).toEqual(["0px", "0px"]);
    });
});

describe("makeDraggable without a container", () => {
    test("is not limited at all", () => {
        const el = setup([300, 200], { left: 150, top: 100, width: 100, height: 80 }, undefined);
        expect(drag(el, -200, -200)).toEqual(["-50px", "-100px"]);
    });
});
