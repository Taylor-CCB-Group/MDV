import { describe, expect, it } from "vitest";
import { bestIdColumn, buildSampleCsv, summaryDatasourceName, summaryFieldName } from "./sampleSummary";

describe("summaryDatasourceName", () => {
    it("uses samples when that table is not already present", () => {
        expect(summaryDatasourceName(["cells", "genes"], false)).toBe("samples");
    });

    it("keeps an existing user table named samples", () => {
        expect(summaryDatasourceName(["cells", "samples"], false)).toBe("sample_summary");
    });

    it("reuses samples when it is already the summary", () => {
        expect(summaryDatasourceName(["cells", "samples"], true)).toBe("samples");
    });
});

describe("bestIdColumn", () => {
    it("picks the column that matches the most region ids", () => {
        expect(
            bestIdColumn(["a", "b"], [
                { field: "gene", values: ["TP53"] },
                { field: "sample", values: ["a", "b", "c"] },
            ]),
        ).toBe("sample");
    });
});

describe("buildSampleCsv", () => {
    it("quotes commas and renames reserved user fields", () => {
        expect(summaryFieldName("n_cells")).toBe("sample_n_cells");
        expect(
            buildSampleCsv([
                { region_id: "a,b", n_cells: 2 },
                { region_id: "c", n_cells: 0, note: "ok" },
            ]),
        ).toBe('region_id,n_cells,note\n"a,b",2,\nc,0,ok');
    });
});
