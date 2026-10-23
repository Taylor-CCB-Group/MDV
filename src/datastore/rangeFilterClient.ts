import { evaluateFilterSpec, type FilterRequest, type FilterResponse, type FilterSpec } from "./rangeFilterEval";

type Job = {
    spec: FilterSpec;
    size: number;
    resolve: (mask: Uint8Array | null) => void;
    reject: (error: unknown) => void;
};

/**
 * Evaluates filter specs on a dedicated worker for one dimension, latest-wins: at most one
 * request is in flight, and while it runs only the newest further request is kept. During rapid
 * brushing the worker therefore never works through a backlog of stale specs.
 *
 * A request that is replaced before it is dispatched resolves to `null`, which `commitFilter`
 * treats as no result. The in-flight one still lands and, being newer than what is applied, is
 * shown while the newest runs.
 *
 * Each request gets a fresh output buffer, so a result is never overwritten while it is being
 * committed. If the worker fails to load or crashes, outstanding and later requests are evaluated
 * on the main thread instead; so is a request that can't be posted.
 */
export class RangeFilterClient {
    private worker: Worker | null = null;
    private inFlight: (Job & { id: number; output: SharedArrayBuffer }) | null = null;
    private queued: Job | null = null;
    private nextId = 0;
    private failed = false;

    constructor(private readonly createWorker: () => Worker) {}

    /** Resolves to an include mask for rows `[0, size)`, or `null` if superseded first. */
    evaluate(spec: FilterSpec, size: number): Promise<Uint8Array | null> {
        return new Promise((resolve, reject) => {
            const job: Job = { spec, size, resolve, reject };
            if (this.failed) {
                evaluateOnMainThread(job);
            } else if (this.inFlight) {
                this.queued?.resolve(null);
                this.queued = job;
            } else {
                this.dispatch(job);
            }
        });
    }

    /** Drops the queued request, if any; an in-flight one cannot be recalled. */
    cancelQueued() {
        this.queued?.resolve(null);
        this.queued = null;
    }

    /** Stops the worker; outstanding requests resolve to `null`. */
    terminate() {
        this.worker?.terminate();
        this.worker = null;
        this.inFlight?.resolve(null);
        this.inFlight = null;
        this.cancelQueued();
    }

    private dispatch(job: Job) {
        const id = ++this.nextId;
        const output = new SharedArrayBuffer(job.size);
        this.inFlight = { ...job, id, output };
        const request: FilterRequest = { id, spec: job.spec, size: job.size, output };
        try {
            this.getWorker().postMessage(request);
        } catch (e) {
            // e.g. a spec that can't be cloned: nothing was sent, so evaluate it here instead
            console.warn("could not post range filter to its worker; filtering on the main thread", e);
            this.inFlight = null;
            evaluateOnMainThread(job);
            this.dispatchQueued();
        }
    }

    private dispatchQueued() {
        const next = this.queued;
        this.queued = null;
        if (next) this.dispatch(next);
    }

    private getWorker() {
        if (!this.worker) {
            const worker = this.createWorker();
            worker.onmessage = (event: MessageEvent<FilterResponse>) => this.onResponse(event.data);
            worker.onerror = (event) => this.onWorkerError(event);
            this.worker = worker;
        }
        return this.worker;
    }

    private onResponse(response: FilterResponse) {
        const job = this.inFlight;
        if (!job || response.id !== job.id) return;
        this.inFlight = null;
        if (response.error !== undefined) {
            job.reject(new Error(response.error));
        } else {
            job.resolve(new Uint8Array(job.output));
        }
        this.dispatchQueued();
    }

    private onWorkerError(event: ErrorEvent) {
        console.warn("range filter worker failed; filtering on the main thread instead", event.message);
        event.preventDefault();
        this.failed = true;
        this.worker?.terminate();
        this.worker = null;
        const inFlight = this.inFlight;
        const queued = this.queued;
        this.inFlight = null;
        this.queued = null;
        if (inFlight) {
            // with a newer request queued, the in-flight one is already stale
            if (queued) inFlight.resolve(null);
            else evaluateOnMainThread(inFlight);
        }
        if (queued) evaluateOnMainThread(queued);
    }
}

function evaluateOnMainThread(job: Job) {
    try {
        const mask = new Uint8Array(job.size);
        evaluateFilterSpec(job.spec, job.size, mask);
        job.resolve(mask);
    } catch (e) {
        job.reject(e);
    }
}
