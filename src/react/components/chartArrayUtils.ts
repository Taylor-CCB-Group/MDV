export type ChartArrayRect = {
    x: number;
    y: number;
    width: number;
    height: number;
};

export function roundRect(rect: ChartArrayRect): ChartArrayRect {
    const x = Math.round(rect.x);
    const y = Math.round(rect.y);
    const right = Math.round(rect.x + rect.width);
    const bottom = Math.round(rect.y + rect.height);
    return {
        x,
        y,
        width: Math.max(0, right - x),
        height: Math.max(0, bottom - y),
    };
}

export function measureCellRectsRelativeToRoot(
    root: HTMLElement,
    cells: readonly (HTMLElement | null)[],
): ChartArrayRect[] {
    const rootRect = root.getBoundingClientRect();
    return cells.map((cell) => {
        if (!cell) {
            return { x: 0, y: 0, width: 0, height: 0 };
        }
        const cellRect = cell.getBoundingClientRect();
        return roundRect({
            x: cellRect.left - rootRect.left,
            y: cellRect.top - rootRect.top,
            width: cellRect.width,
            height: cellRect.height,
        });
    });
}

export function measureRootSize(root: HTMLElement): ChartArrayRect {
    return roundRect({
        x: 0,
        y: 0,
        width: root.offsetWidth,
        height: root.offsetHeight,
    });
}

/** Cell indices with any intersection in the scroll root. */
export function indicesFromIntersectionEntries(
    entries: IntersectionObserverEntry[],
    cellCount: number,
): number[] {
    const visible = new Set<number>();
    for (const entry of entries) {
        const index = Number((entry.target as HTMLElement).dataset.chartArrayIndex);
        if (!Number.isInteger(index) || index < 0 || index >= cellCount) continue;
        if (entry.isIntersecting) {
            visible.add(index);
        }
    }
    return [...visible].sort((a, b) => a - b);
}
