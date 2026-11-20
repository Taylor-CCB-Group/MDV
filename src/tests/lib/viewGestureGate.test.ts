import { describe, expect, test } from "vitest";
import { hasViewModifier, isViewDrag } from "@/lib/viewGestureGate";

/** A drag as mjolnir reports it, started by a pointer move with these keys and buttons held. */
function drag(
    src: { metaKey?: boolean; ctrlKey?: boolean; shiftKey?: boolean; buttons?: number },
    rightButton = false,
) {
    return { srcEvent: new MouseEvent("pointermove", src), rightButton, middleButton: false };
}

describe("on a Mac", () => {
    test("⌘ moves the view", () => {
        expect(isViewDrag(drag({ metaKey: true, buttons: 1 }), true)).toBe(true);
    });

    test("Ctrl does not, as Ctrl-click is a right-click there", () => {
        expect(hasViewModifier({ ctrlKey: true }, true)).toBe(false);
        expect(isViewDrag(drag({ ctrlKey: true, buttons: 1 }), true)).toBe(false);
    });
});

describe("off a Mac", () => {
    test("Ctrl moves the view", () => {
        expect(isViewDrag(drag({ ctrlKey: true, buttons: 1 }), false)).toBe(true);
    });
});

describe("whatever the platform", () => {
    test.each([true, false])("a plain drag, or one with only Shift, draws (mac: %s)", (mac) => {
        expect(isViewDrag(drag({ buttons: 1 }), mac)).toBe(false);
        expect(isViewDrag(drag({ shiftKey: true, buttons: 1 }), mac)).toBe(false);
    });

    test.each([true, false])("the right or middle button moves the view (mac: %s)", (mac) => {
        expect(isViewDrag(drag({ buttons: 2 }), mac)).toBe(true);
        expect(isViewDrag(drag({ buttons: 4 }), mac)).toBe(true);
        expect(isViewDrag(drag({}, true), mac)).toBe(true);
    });
});
