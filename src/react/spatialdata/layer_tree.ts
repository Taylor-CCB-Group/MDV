import type { RenderStackEntry } from "@spatialdata/layers";

/** One row of the layer dialog's tree, in display order (which is draw order, bottom first). */
export type LayerTreeRowInfo = {
    id: string;
    depth: number;
    parentId: string | null;
    isGroup: boolean;
    /** A group showing its children. */
    expanded: boolean;
};

/**
 * The rows for `entries`, a group followed by its children when it is expanded. Children
 * of `hideChildrenOf` are left out, as they are while that group is being dragged: they
 * move with it.
 */
export function flattenLayerTree(
    entries: readonly RenderStackEntry[],
    collapsed: ReadonlySet<string>,
    hideChildrenOf: string | null = null,
): LayerTreeRowInfo[] {
    const rows: LayerTreeRowInfo[] = [];
    const visit = (list: readonly RenderStackEntry[], depth: number, parentId: string | null) => {
        for (const entry of list) {
            const isGroup = entry.kind === "group";
            const expanded = isGroup && !collapsed.has(entry.id) && entry.id !== hideChildrenOf;
            rows.push({ id: entry.id, depth, parentId, isGroup, expanded });
            if (entry.kind === "group" && expanded) visit(entry.children, depth + 1, entry.id);
        }
    };
    visit(entries, 0, null);
    return rows;
}

export type LayerTreeDrop = { depth: number; parentId: string | null; afterId: string | null };

/**
 * Where a row dropped over `overId` goes, with `depthOffset` levels of sideways drag
 * (dnd-kit's sortable-tree approach). It takes `overId`'s place in the list, and its
 * depth is the dragged one, clamped to what the rows around it allow: no deeper than
 * one inside an expanded group just above it, no shallower than the row below it.
 * Returns null when the row is not in the list.
 */
export function projectLayerTreeDrop(
    rows: readonly LayerTreeRowInfo[],
    activeId: string,
    overId: string,
    depthOffset: number,
): LayerTreeDrop | null {
    const activeIndex = rows.findIndex((row) => row.id === activeId);
    const overIndex = rows.findIndex((row) => row.id === overId);
    if (activeIndex === -1 || overIndex === -1) return null;
    const active = rows[activeIndex];
    const moved = rows.slice();
    moved.splice(activeIndex, 1);
    moved.splice(overIndex, 0, active);

    const previous = moved[overIndex - 1];
    const next = moved[overIndex + 1];
    const maxDepth = previous ? (previous.isGroup && previous.expanded ? previous.depth + 1 : previous.depth) : 0;
    const minDepth = next ? next.depth : 0;
    const depth = Math.max(minDepth, Math.min(maxDepth, active.depth + depthOffset));

    let parentId: string | null = null;
    let afterId: string | null = null;
    for (let i = overIndex - 1; i >= 0; i--) {
        const row = moved[i];
        if (row.depth === depth && afterId === null) afterId = row.id;
        if (row.depth < depth) {
            parentId = row.id;
            break;
        }
    }
    return { depth, parentId, afterId };
}
