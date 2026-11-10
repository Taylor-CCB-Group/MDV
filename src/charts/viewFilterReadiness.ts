/**
 * While a view loads, an active link is only worth resolving once the dataSource it links to has the filters
 * saved with the view. Those filters aren't view state: each is held in the config of the chart that applies it
 * (a Selection Dialog, a histogram's brush...), and is applied by that chart as part of its own initialisation, at
 * no particular time relative to the charts that consume the link. A chart that resolves the link before then is
 * built on - and loads columns for - rows that are about to be filtered out.
 *
 * So charts report when their stored filters have been applied ({@link StoredFilters}), a view load keeps track of
 * which dataSources still have charts that haven't ({@link ViewFilterReadiness}), and a link waits for the
 * dataSource it links to before it is resolved for a chart that is being added.
 *
 * Nothing here orders the charts: a chart that applies a filter can itself depend on a link (a Selection Dialog
 * filtering on a linked column), and it simply waits like any other, so the order follows from who is waiting on
 * whom. That only fails if the waiting goes round in a circle ({@link linkCycleMembers}), and for anything that
 * never reports there is a time limit. In either case charts carry on as they did before any of this existed:
 * resolving links straight away, and changing if a filter then lands.
 */

/** How long a chart is given to apply its stored filters before it is treated as having done so. */
export const STORED_FILTERS_TIMEOUT_MS = 5000;
/** How long a view load waits for a dataSource's charts before links to it are resolved regardless. */
export const VIEW_FILTERS_TIMEOUT_MS = 10000;

/**
 * Waits for every promise to settle, but no longer than `timeoutMs`.
 * @returns whether they all settled in time
 */
export function allSettledWithin(promises: Iterable<PromiseLike<unknown>>, timeoutMs: number) {
    return new Promise<boolean>((resolve) => {
        const timeout = setTimeout(() => resolve(false), timeoutMs);
        Promise.allSettled(promises).then(() => {
            clearTimeout(timeout);
            resolve(true);
        });
    });
}

/**
 * The dataSources that are on a cycle of links.
 *
 * @param waitsOn for each dataSource, the dataSources it has rows-as-columns links to: its charts may wait for
 * those to be ready, and it is not ready itself until its charts are. Waiting on a dataSource that (directly or not)
 * waits on you would never end, so such dataSources must not be waited for.
 */
export function linkCycleMembers(waitsOn: Record<string, string[]>) {
    const members = new Set<string>();
    for (const start of Object.keys(waitsOn)) {
        // `start` is on a cycle if following links from it leads back to it
        const seen = new Set<string>();
        const queue = [...(waitsOn[start] ?? [])];
        while (queue.length > 0) {
            const name = queue.pop();
            if (name === undefined || seen.has(name)) continue;
            if (name === start) {
                members.add(start);
                break;
            }
            seen.add(name);
            queue.push(...(waitsOn[name] ?? []));
        }
    }
    return members;
}

/**
 * The filters saved in one chart's config, and whether they have been applied to the dataStore yet.
 *
 * A chart that applies its filters in its constructor has nothing to do here. One that applies them later (in
 * React effects, say) declares each with {@link expect} while it is being constructed, and reports each with
 * {@link applied} once it has called `Dimension.filter()` for it.
 */
export class StoredFilters {
    /** filters that have not been applied yet; `reported` once the chart has asked for it to be */
    private outstanding = new Map<string, PromiseWithResolvers<void> & { reported: boolean }>();

    /** Declares a filter that is saved in the chart's config and will be applied later. Call during construction. */
    expect(key: string) {
        if (this.outstanding.has(key)) return;
        this.outstanding.set(key, { ...Promise.withResolvers<void>(), reported: false });
    }

    /**
     * Reports that the filter declared with {@link expect} has been applied.
     * @param landed what `Dimension.filter()` returned: a filter evaluated on a worker is not applied until then
     */
    applied(key: string, landed?: PromiseLike<unknown>) {
        const filter = this.outstanding.get(key);
        if (!filter || filter.reported) return;
        filter.reported = true;
        const done = () => {
            this.outstanding.delete(key);
            filter.resolve();
        };
        // whether or not it landed: a filter that was superseded or failed is not going to be applied by waiting
        Promise.resolve(landed).then(done, done);
    }

    /** Resolves once every expected filter has been applied - or the time limit has passed. */
    whenApplied(timeoutMs = STORED_FILTERS_TIMEOUT_MS) {
        const waitingFor = [...this.outstanding.values()].map((filter) => filter.promise);
        if (waitingFor.length === 0) return Promise.resolve();
        return new Promise<void>((resolve) => {
            const timeout = setTimeout(() => {
                console.warn("Chart did not apply its stored filters in time:", [...this.outstanding.keys()]);
                resolve();
            }, timeoutMs);
            Promise.all(waitingFor).then(() => {
                clearTimeout(timeout);
                resolve();
            });
        });
    }
}

/**
 * For one load of a view: which of its charts have yet to apply their stored filters, by dataSource.
 */
export class ViewFilterReadiness {
    /** the view's charts that have not settled, by dataSource name */
    private unsettled = new Map<string, Set<unknown>>();
    private waiting = new Map<string, PromiseWithResolvers<void>>();
    private timeout: ReturnType<typeof setTimeout> | undefined;

    /**
     * @param charts the configs of the charts in the view, by dataSource name - identified by object identity
     * @param ungated dataSources that are never waited for: see {@link linkCycleMembers}
     */
    constructor(
        charts: Record<string, unknown[]>,
        private ungated: ReadonlySet<string> = new Set(),
        timeoutMs = VIEW_FILTERS_TIMEOUT_MS,
    ) {
        for (const [dataSource, configs] of Object.entries(charts)) {
            if (configs.length > 0) this.unsettled.set(dataSource, new Set(configs));
        }
        if (this.unsettled.size > 0) {
            this.timeout = setTimeout(() => this.release(), timeoutMs);
        }
    }

    /**
     * Records that a chart of the view has applied its stored filters - or isn't going to, having failed to load.
     * Ignored for a chart that isn't one the view was loaded with (one the user added since).
     */
    settle(dataSource: string, chart: unknown) {
        const charts = this.unsettled.get(dataSource);
        if (!charts?.delete(chart) || charts.size > 0) return;
        this.unsettled.delete(dataSource);
        this.waiting.get(dataSource)?.resolve();
        this.waiting.delete(dataSource);
        if (this.unsettled.size === 0) clearTimeout(this.timeout);
    }

    /** Whether the filters saved with the view have been applied to `dataSource` - see {@link settle}. */
    isReady(dataSource: string) {
        return this.ungated.has(dataSource) || !this.unsettled.has(dataSource);
    }

    /** Resolves once {@link isReady}. */
    whenReady(dataSource: string) {
        if (this.isReady(dataSource)) return Promise.resolve();
        let waiting = this.waiting.get(dataSource);
        if (!waiting) {
            waiting = Promise.withResolvers<void>();
            this.waiting.set(dataSource, waiting);
        }
        return waiting.promise;
    }

    /** Stops waiting for anything, because this load of the view has been replaced by another. */
    abandon() {
        clearTimeout(this.timeout);
        this.unsettled.clear();
        for (const waiting of this.waiting.values()) waiting.resolve();
        this.waiting.clear();
    }

    /** Stops waiting for anything: whatever has not settled by now is not going to hold up the rest of the view. */
    private release() {
        console.warn("Charts did not apply their stored filters in time, for dataSources:", [...this.unsettled.keys()]);
        this.unsettled.clear();
        for (const waiting of this.waiting.values()) waiting.resolve();
        this.waiting.clear();
    }
}
