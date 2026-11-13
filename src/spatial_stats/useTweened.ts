import { useEffect, useRef, useState } from "react";

/** Cubic ease-out from `from` to `to` at progress `t ∈ [0, 1]`. */
export function tweenValue(from: number, to: number, t: number): number {
    const u = Math.min(1, Math.max(0, t));
    return from + (to - from) * (1 - (1 - u) ** 3);
}

const reducedMotion = () =>
    typeof window !== "undefined" && (window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ?? false);

/**
 * `target`, eased over `ms` whenever it changes (an axis range, say), starting from wherever the
 * previous transition had got to. Jumps straight there under `prefers-reduced-motion`.
 */
export function useTweened(target: number, ms = 300): number {
    const [value, setValue] = useState(target);
    const current = useRef(target);
    useEffect(() => {
        const from = current.current;
        if (from === target || reducedMotion() || ms <= 0) {
            current.current = target;
            setValue(target);
            return;
        }
        const start = performance.now();
        let frame = requestAnimationFrame(function step(now) {
            const t = (now - start) / ms;
            current.current = tweenValue(from, target, t);
            setValue(current.current);
            if (t < 1) frame = requestAnimationFrame(step);
        });
        return () => cancelAnimationFrame(frame);
    }, [target, ms]);
    return value;
}
