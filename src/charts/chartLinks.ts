/** The parts of a chart link that matter here; links have more fields, which differ by type. */
export type ChartLink = {
    id?: string;
    type?: string;
    source_chart?: string;
    target_charts?: string[];
    linked_charts?: string[];
};

const targetsOf = (link: ChartLink) => link.target_charts ?? link.linked_charts;

/** A copy of the links that a chart is the source or a target of. */
export function snapshotChartLinks(links: ChartLink[] | undefined, chartId: string): ChartLink[] {
    return (links ?? [])
        .filter((link) => link.source_chart === chartId || targetsOf(link)?.includes(chartId))
        .map((link) => JSON.parse(JSON.stringify(link)));
}

const isSameLink = (a: ChartLink, b: ChartLink) =>
    a.type === b.type && a.source_chart === b.source_chart && (a.id === undefined || b.id === undefined || a.id === b.id);

/**
 * Puts back links that were dropped while a chart was replaced by one with the same id: a link that is
 * still there gets the chart back as a target, and one that was removed is added again and set up.
 */
export function restoreChartLinks(
    links: ChartLink[],
    snapshot: ChartLink[],
    chartId: string,
    setUpLink: (link: ChartLink) => void,
) {
    for (const saved of snapshot) {
        const existing = links.find((link) => isSameLink(link, saved));
        if (existing) {
            const targets = targetsOf(existing);
            if (targets && targetsOf(saved)?.includes(chartId) && !targets.includes(chartId)) {
                targets.push(chartId);
            }
            continue;
        }
        links.push(saved);
        try {
            setUpLink(saved);
        } catch (error) {
            console.error(`Could not set up a link for chart '${chartId}' again`, error, saved);
        }
    }
}
