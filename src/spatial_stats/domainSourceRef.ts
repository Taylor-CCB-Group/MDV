import type { SimpleFeatureCollection } from "@deck.gl-community/editable-layers";
import type { Gate } from "@/react/gates/types";
import { type DomainSource, featureRings, gateDomainSource, selectionDomainSource } from "./domainSource";

/**
 * A serialisable pointer to a domain source — what a chart (or later a dataflow node) stores in its
 * config. Resolving it reads MobX state (a chart's `selectionFeatureCollection`, the gate map), so a
 * `reaction`/`autorun` that resolves it re-runs when the selection or gate changes.
 */
export type DomainSourceRef =
    | { kind: "selection"; chartId: string }
    | { kind: "gate"; gateId: string }
    /** The whole region (its coordinate-system extent); `region` unset = the datastore's only region. */
    | { kind: "region"; region?: string };

/** A chart as far as domain sources care: one that may hold a live selection. */
export type SelectionHost = {
    /** Compared by identity only. */
    dataStore: object;
    config: {
        id: string;
        title?: string;
        /** Position fields; non-string field specs mean there is no usable x/y pair. */
        param: readonly unknown[];
        region?: string;
        selectionFeatureCollection?: SimpleFeatureCollection;
    };
};

export type DomainSourceEnv = {
    dataStore: object;
    charts: Iterable<SelectionHost>;
    gates: ReadonlyMap<string, Gate>;
    /** SpatialData region keys; empty when the datastore has none (the whole extent is then one source). */
    regions: readonly string[];
    /** The whole-region source; must return the same object for the same region. */
    regionSource(region?: string): DomainSource;
};

export function refKey(ref: DomainSourceRef) {
    if (ref.kind === "selection") return `selection:${ref.chartId}`;
    if (ref.kind === "gate") return `gate:${ref.gateId}`;
    return `region:${ref.region ?? ""}`;
}

const regionRefs = (regions: readonly string[]): DomainSourceRef[] =>
    regions.length > 0 ? regions.map((region) => ({ kind: "region", region })) : [{ kind: "region" }];

/** What a chart uses when it hasn't been pointed anywhere: a drawn selection, else the whole region. */
export function defaultDomainSourceRef(env: DomainSourceEnv): DomainSourceRef {
    return listDomainSources(env).find((s) => s.ready && s.ref.kind === "selection")?.ref ?? regionRefs(env.regions)[0];
}

/** Every source available for `dataStore`: charts on it that have a selection layer, and its gates. */
export function listDomainSources({ dataStore, charts, gates, regions }: DomainSourceEnv) {
    const out: { ref: DomainSourceRef; label: string; ready: boolean }[] = [];
    for (const c of charts) {
        if (c.dataStore !== dataStore || !("selectionFeatureCollection" in c.config)) continue;
        const fc = c.config.selectionFeatureCollection;
        out.push({
            ref: { kind: "selection", chartId: c.config.id },
            label: `Selection in "${c.config.title || c.config.id}"`,
            ready: !!fc && featureRings(fc).length > 0,
        });
    }
    for (const g of gates.values())
        out.push({ ref: { kind: "gate", gateId: g.id }, label: `Gate "${g.name}"`, ready: true });
    for (const ref of regionRefs(regions))
        out.push({
            ref,
            label: ref.kind === "region" && ref.region ? `Whole region "${ref.region}"` : "Whole extent",
            ready: true,
        });
    return out;
}

/** One `DomainSource` per underlying geometry object, so consumers can key work on its identity: MDV
 *  replaces the selection collection (it is `observable.ref`) and the gate object on every edit. */
const resolved = new WeakMap<object, DomainSource>();
function cached(key: object, make: () => DomainSource | undefined) {
    let s = resolved.get(key);
    if (!s) {
        s = make();
        if (s) resolved.set(key, s);
    }
    return s;
}

/**
 * The source `ref` points at now. A selection with nothing drawn, or a gate that has been deleted,
 * falls back to the whole region (the selection's chart's region, else the datastore's), so there is
 * always a window. Same input geometry → same object.
 */
export function resolveDomainSourceRef(ref: DomainSourceRef, env: DomainSourceEnv): DomainSource {
    if (ref.kind === "region") return env.regionSource(ref.region);
    if (ref.kind === "gate") {
        const gate = env.gates.get(ref.gateId);
        return gate ? (cached(gate, () => gateDomainSource(gate)) ?? env.regionSource()) : env.regionSource();
    }
    for (const c of env.charts) {
        if (c.config.id !== ref.chartId) continue;
        const fc = c.config.selectionFeatureCollection;
        return (fc && cached(fc, () => selectionDomainSource(c.config))) || env.regionSource(c.config.region);
    }
    return env.regionSource();
}
