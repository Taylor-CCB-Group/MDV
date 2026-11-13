import { act, renderHook } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import { useTweened } from "@/spatial_stats/useTweened";

/** A hand-cranked animation-frame clock: `frame(ms)` advances time and runs the queued callbacks. */
function installFrameClock() {
    let now = 0;
    let queue: FrameRequestCallback[] = [];
    vi.spyOn(performance, "now").mockImplementation(() => now);
    vi.stubGlobal("requestAnimationFrame", (cb: FrameRequestCallback) => queue.push(cb));
    vi.stubGlobal("cancelAnimationFrame", () => {
        queue = [];
    });
    return {
        frame(ms: number) {
            now += ms;
            const run = queue;
            queue = [];
            act(() => {
                for (const cb of run) cb(now);
            });
        },
        pending: () => queue.length,
    };
}

describe("useTweened", () => {
    beforeEach(() => vi.stubGlobal("matchMedia", () => ({ matches: false })));
    afterEach(() => {
        vi.unstubAllGlobals();
        vi.restoreAllMocks();
    });

    test("eases from the old value to the new one over the duration, then stops", () => {
        const clock = installFrameClock();
        const { result, rerender } = renderHook(({ v }) => useTweened(v, 300), { initialProps: { v: 2 } });
        expect(result.current).toBe(2);
        rerender({ v: 15 });
        expect(result.current).toBe(2); // first frame not yet run
        clock.frame(100);
        const a = result.current;
        clock.frame(100);
        const b = result.current;
        expect(a).toBeGreaterThan(2);
        expect(b).toBeGreaterThan(a);
        expect(b).toBeLessThan(15);
        clock.frame(200);
        expect(result.current).toBe(15);
        expect(clock.pending()).toBe(0);
    });

    test("a new target mid-flight starts from where the value had got to", () => {
        const clock = installFrameClock();
        const { result, rerender } = renderHook(({ v }) => useTweened(v, 300), { initialProps: { v: 2 } });
        rerender({ v: 15 });
        clock.frame(100);
        const mid = result.current;
        rerender({ v: 3 });
        clock.frame(1);
        expect(Math.abs(result.current - mid)).toBeLessThan(1); // no jump back to 2 or on to 15
        clock.frame(400);
        expect(result.current).toBe(3);
    });

    test("jumps straight to the target under prefers-reduced-motion", () => {
        installFrameClock();
        vi.stubGlobal("matchMedia", () => ({ matches: true }));
        const { result, rerender } = renderHook(({ v }) => useTweened(v, 300), { initialProps: { v: 2 } });
        rerender({ v: 15 });
        expect(result.current).toBe(15);
    });
});
