import { beforeEach, describe, expect, test } from "vitest";
import { needsOwnTab, showOwnTabNotice, type FrameState } from "@/utils/ownTabNotice";

function frameState(framed: boolean, crossOriginIsolated: boolean): FrameState {
    const state: FrameState = { top: null, crossOriginIsolated };
    state.top = framed ? {} : state;
    return state;
}

describe("needsOwnTab", () => {
    test("asks for a tab when framed without cross-origin isolation", () => {
        expect(needsOwnTab(frameState(true, false))).toBe(true);
    });

    test("leaves a framed page alone when it is isolated", () => {
        expect(needsOwnTab(frameState(true, true))).toBe(false);
    });

    // A tab of its own cannot fix a top-level page without isolation, so do not offer one.
    test("leaves a top-level page alone, isolated or not", () => {
        expect(needsOwnTab(frameState(false, true))).toBe(false);
        expect(needsOwnTab(frameState(false, false))).toBe(false);
    });
});

describe("showOwnTabNotice", () => {
    beforeEach(() => {
        document.body.innerHTML = "";
    });

    test("links to the same address in a new tab", () => {
        showOwnTabNotice(document, "https://app.example.org/project/7");

        const link = document.querySelector("a");
        expect(link?.getAttribute("href")).toBe("https://app.example.org/project/7");
        expect(link?.getAttribute("target")).toBe("_blank");
        expect(link?.getAttribute("rel")).toBe("noopener");
        expect(document.body.textContent).toContain("new tab");
    });

    test("shows one notice however often it is called", () => {
        showOwnTabNotice(document, "https://app.example.org/");
        showOwnTabNotice(document, "https://app.example.org/");

        expect(document.querySelectorAll("a")).toHaveLength(1);
    });
});
