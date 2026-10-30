import { restoreChartLinks, snapshotChartLinks, type ChartLink } from "@/charts/chartLinks";
import { describe, expect, test } from "vitest";

describe("chart links across replacing a chart", () => {
    test("snapshots only the links the chart takes part in", () => {
        const links: ChartLink[] = [
            { id: "a", type: "color_by_column", source_chart: "c1", target_charts: ["c2"] },
            { id: "b", type: "color_by_column", source_chart: "c3", target_charts: ["c4", "c1"] },
            { id: "c", type: "viewState", linked_charts: ["c1", "c5"] },
            { id: "d", type: "color_by_column", source_chart: "c3", target_charts: ["c4"] },
        ];
        expect(snapshotChartLinks(links, "c1").map((l) => l.id)).toEqual(["a", "b", "c"]);
        expect(snapshotChartLinks(undefined, "c1")).toEqual([]);
    });

    test("adds the chart back to a link that is still there", () => {
        const links: ChartLink[] = [{ id: "b", type: "color_by_column", source_chart: "c3", target_charts: ["c4", "c1"] }];
        const snapshot = snapshotChartLinks(links, "c1");
        links[0].target_charts?.splice(1, 1);
        const setUp: ChartLink[] = [];
        restoreChartLinks(links, snapshot, "c1", (l) => setUp.push(l));
        expect(links).toHaveLength(1);
        expect(links[0].target_charts).toEqual(["c4", "c1"]);
        expect(setUp).toEqual([]);
    });

    test("restores and sets up a link that was removed", () => {
        const links: ChartLink[] = [{ id: "a", type: "color_by_column", source_chart: "c1", target_charts: ["c2"] }];
        const snapshot = snapshotChartLinks(links, "c1");
        links.length = 0;
        const setUp: ChartLink[] = [];
        restoreChartLinks(links, snapshot, "c1", (l) => setUp.push(l));
        expect(links).toEqual(snapshot);
        expect(setUp).toEqual(snapshot);
    });

    test("does not duplicate a link when its chart was not removed from it", () => {
        const links: ChartLink[] = [{ id: "c", type: "viewState", linked_charts: ["c1", "c5"] }];
        const snapshot = snapshotChartLinks(links, "c1");
        restoreChartLinks(links, snapshot, "c1", () => {});
        expect(links).toHaveLength(1);
        expect(links[0].linked_charts).toEqual(["c1", "c5"]);
    });
});
