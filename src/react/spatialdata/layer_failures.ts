/**
 * Turning a SpatialData chart's render errors into something a user can see.
 *
 * Errors arrive two ways:
 * - deck's `onError`, for a layer whose lifecycle threw. deck passes the layer, so the
 *   render-stack entry is known.
 * - luma's device `onError`, for WebGPU validation errors. These are asynchronous and
 *   name GPU resources, not layers. One invalid resource also invalidates the frame's
 *   whole command buffer, so a single bad layer blanks every layer in the Deck, every
 *   frame, until it stops drawing.
 *
 * For the second kind, the first message is usually the cause ("Texture size … exceeded")
 * and names no layer; the ones that follow name luma resources that do carry the layer
 * id ("Model(<layer id>).draw"). So an unattributed message is held briefly as the
 * likely cause of the next attributed one. Consequences ("… is invalid due to a previous
 * error") are not reported on their own.
 */

/** A render-stack entry that failed and was turned off. */
export type LayerFailure = {
    /** The root cause, as the first message the failure produced. */
    reason: string;
    /** The message that tied the failure to this layer, when it differs from the reason. */
    detail?: string;
    source: "webgpu" | "deck";
    /**
     * Where it was thrown, for errors thrown in JS (deck's). WebGPU validation errors arrive
     * asynchronously, so their stack only shows how luma reported them, and is left out.
     */
    stack?: string;
    /** How many errors were attributed to it. */
    count: number;
    /** Dismissed from the chart's error count (its menu pip); the layer dialog still shows the failure. */
    dismissed: boolean;
};

/** An error that could not be tied to a layer. `stack` as for `LayerFailure`. */
export type CanvasError = { message: string; stack?: string; count: number; dismissed: boolean };

/** How long an unattributed message waits to be claimed as the cause of an attributed one. */
export const CAUSE_WINDOW_MS = 1000;
const CONSEQUENCE = /is invalid due to a previous error/;
const MAX_LOGGED = 50;

function isBoundary(ch: string | undefined, chars: RegExp): boolean {
    return ch === undefined || chars.test(ch);
}

/**
 * The render-stack entry a deck layer id or luma resource label belongs to: the longest
 * entry id that appears as a whole label prefix. Layer ids are `<entry id>[-<sublayer>…]`.
 */
export function entryForLabel(text: string, entryIds: readonly string[]): string | undefined {
    let best: string | undefined;
    for (const id of entryIds) {
        if (!id || (best && id.length <= best.length)) continue;
        for (let at = text.indexOf(id); at >= 0; at = text.indexOf(id, at + 1)) {
            if (isBoundary(text[at - 1], /["(\s]/) && isBoundary(text[at + id.length], /[-"/)\s]/)) {
                best = id;
                break;
            }
        }
    }
    return best;
}

/** Collapses resource numbering so "Texture-11" and "Texture-12" log as one message. */
function logKey(message: string) {
    return message.replace(/\d+/g, "#");
}

export type LayerFailureSink = {
    entryIds: () => readonly string[];
    fail: (entryId: string, failure: Omit<LayerFailure, "count" | "dismissed">) => void;
    canvas: (message: string, stack?: string) => void;
};

/** One per chart: attributes errors, pairs them with their cause, and logs each kind once. */
export class LayerFailureTracker {
    private pending: { message: string; at: number } | null = null;
    private timer: ReturnType<typeof setTimeout> | null = null;
    private readonly logged = new Set<string>();

    constructor(
        private readonly sink: LayerFailureSink,
        private readonly now: () => number = () => performance.now(),
        private readonly log: (error: Error) => void = (error) => console.error(error),
    ) {}

    /** From luma's device `onError`. */
    device(error: Error): void {
        this.logOnce(error);
        const { message } = error;
        const entryId = entryForLabel(message, this.sink.entryIds());
        if (entryId) {
            const cause = this.takeCause();
            this.sink.fail(entryId, {
                reason: cause ?? message,
                ...(cause ? { detail: message } : {}),
                source: "webgpu",
            });
            return;
        }
        if (CONSEQUENCE.test(message) || this.pending) return;
        this.pending = { message, at: this.now() };
        this.timer = setTimeout(() => this.flush(), CAUSE_WINDOW_MS);
    }

    /** From deck's `onError`, which passes the failing layer. */
    layer(error: Error, layerId: string | undefined): void {
        this.logOnce(error);
        const entryId = layerId ? entryForLabel(layerId, this.sink.entryIds()) : undefined;
        const stack = error.stack ? { stack: error.stack } : {};
        if (entryId) this.sink.fail(entryId, { reason: error.message, source: "deck", ...stack });
        else this.sink.canvas(error.message, error.stack);
    }

    /** Report a held message nothing claimed as a canvas-level error. */
    flush(): void {
        const pending = this.pending;
        this.dispose();
        if (pending) this.sink.canvas(pending.message);
    }

    dispose(): void {
        if (this.timer) clearTimeout(this.timer);
        this.timer = null;
        this.pending = null;
    }

    private takeCause(): string | undefined {
        const pending = this.pending;
        this.dispose();
        if (!pending || this.now() - pending.at > CAUSE_WINDOW_MS) return undefined;
        return pending.message;
    }

    private logOnce(error: Error) {
        const key = logKey(error.message);
        if (this.logged.has(key) || this.logged.size >= MAX_LOGGED) return;
        this.logged.add(key);
        this.log(error);
    }
}
