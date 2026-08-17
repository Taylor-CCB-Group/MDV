import { autorun, runInAction, untracked } from "mobx";
import { useEffect } from "react";

import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import { flattenFields } from "@/lib/columnTypeHelpers";
import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";

import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "../components/SpatialDataMDVReact";

/**
 * The column choices a spatial layer holds in MDV's own vocabulary, before they
 * are reduced to the concrete column names the viewer understands.
 *
 * MDV's column picker can return a `RowsAsColsQuery` — an "active link", whose
 * columns are whatever the linked datasource currently has selected, and which
 * therefore CHANGE while the chart is open. The viewer's `LayerConfig` has no
 * vocabulary for that: `fillColorByColumn.columnName` and `tooltipFields` are
 * plain strings, by design, because a saved Render Stack has to mean the same
 * thing to a reader that has never heard of MDV's links.
 *
 * So the query lives here, on MDV's side of the layer props, and
 * {@link useProjectMdvFieldSpecs} keeps the viewer-facing fields written from
 * it. The spec is the source of truth; the concrete fields are derived, which is
 * why the picker writes a spec even for an ordinary column — a config where the
 * two could disagree would be one where the last writer wins at random.
 *
 * Absent on a layer that predates this (or one whose column was never chosen
 * through the picker), in which case the concrete fields stand on their own and
 * nothing overwrites them.
 */
export type MdvFieldSpecs = {
    fillColorByColumn?: FieldSpec;
    tooltipFields?: FieldSpecs;
};

/**
 * Where {@link MdvFieldSpecs} sits on a spatial entry's props.
 *
 * Serialisation needs nothing: `RowsAsColsQuery.toJSON` puts it in the saved
 * config, and `deserialiseValueRecursive` walks the whole config on load and
 * rebuilds any `{type: "RowsAsColsQuery"}` it finds, wherever it finds it.
 */
export const MDV_FIELD_SPECS_PROP = "mdvFieldSpecs";

/** A layer config as MDV holds it: the viewer's fields, plus the specs they came from. */
export type WithMdvFieldSpecs<T> = T & { mdvFieldSpecs?: MdvFieldSpecs };

/** Accepts either an entry's `props` or the layer config built from them. */
export function mdvFieldSpecsOf(propsOrConfig: unknown): MdvFieldSpecs | undefined {
    if (!propsOrConfig || typeof propsOrConfig !== "object") return undefined;
    const specs = (propsOrConfig as Record<string, unknown>)[MDV_FIELD_SPECS_PROP];
    if (!specs || typeof specs !== "object" || Array.isArray(specs)) return undefined;
    return specs as MdvFieldSpecs;
}

/**
 * Strip MDV's specs from a config on its way to the viewer.
 *
 * A `RowsAsColsQuery` is a live MobX object, and the viewer treats a layer
 * config as inert serialisable state — it may hold it, diff it, or write it into
 * a saved stack. None of those are things to do with a query, and the viewer has
 * no use for one, so the boundary is where it comes off.
 */
export function withoutMdvFieldSpecs<T extends object>(config: T): T {
    if (!(MDV_FIELD_SPECS_PROP in config)) return config;
    const { mdvFieldSpecs: _removed, ...rest } = config as T & WithMdvFieldSpecs<object>;
    return rest as T;
}

function spatialEntries(stack: RenderStack | undefined) {
    return (stack?.entries ?? []).filter(
        (entry): entry is Extract<RenderStackEntry, { kind: "spatial" }> => entry.kind === "spatial",
    );
}

function sameFields(a: string[] | undefined, b: string[]): boolean {
    return !!a && a.length === b.length && a.every((field, index) => field === b[index]);
}

/**
 * One pass: write every spatial layer's concrete column fields from its
 * {@link MdvFieldSpecs}. Returns whether anything actually moved.
 *
 * Called inside the `autorun` below, so resolving a spec is what subscribes to it
 * — `query.fields` is a computed over the linked datasource's current selection.
 * The comparison reads are `untracked` because this pass writes to props it would
 * otherwise be observing: that converges (the second pass writes nothing), but it
 * converges by running twice on every change for no reason.
 */
export function projectMdvFieldSpecs(stack: RenderStack | undefined): boolean {
    let changed = false;

    for (const entry of spatialEntries(stack)) {
        const specs = mdvFieldSpecsOf(entry.props);
        if (!specs) continue;

        const fillColorFields = specs.fillColorByColumn ? flattenFields(specs.fillColorByColumn) : undefined;
        const tooltipFields = specs.tooltipFields ? flattenFields(specs.tooltipFields) : undefined;

        untracked(() => {
            runInAction(() => {
                const props = entry.props as Record<string, unknown>;

                // An empty resolution means the link has not finished initialising,
                // not that the user chose nothing. Leaving the previous column in
                // place keeps the canvas as it was rather than flashing through
                // uncoloured on every reload.
                const nextColumnName = fillColorFields?.[0];
                const fillColorByColumn = props.fillColorByColumn as { columnName?: string } | undefined;
                if (nextColumnName && fillColorByColumn?.columnName !== nextColumnName) {
                    props.fillColorByColumn = { ...fillColorByColumn, columnName: nextColumnName };
                    changed = true;
                }

                if (tooltipFields && !sameFields(props.tooltipFields as string[] | undefined, tooltipFields)) {
                    props.tooltipFields = tooltipFields;
                    changed = true;
                }
            });
        });
    }

    return changed;
}

/**
 * Keep each spatial layer's concrete column fields written from its
 * {@link MdvFieldSpecs}, for as long as the chart is open.
 *
 * This is what makes an active link active: choosing a different gene in the
 * linked datasource lands here as a new column name, and the existing fill-colour
 * path reloads and recolours from it exactly as it would for a column the user had
 * picked by hand.
 */
export function useProjectMdvFieldSpecs(config: SpatialDataMdvReactConfig, chart: SpatialDataMdvReact) {
    useEffect(
        () =>
            autorun(() => {
                if (projectMdvFieldSpecs(config.renderStack)) {
                    chart.bumpRenderStackPropsGeneration();
                }
            }),
        [config, chart],
    );
}
