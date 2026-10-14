import { describe, expect, test } from "vitest";
import {
    type JobsDataSource,
    initialValues,
    paramOptions,
    setParam,
    toJobsDataSource,
    toSubmitParams,
} from "./jobForm";
import type { JobParamSpec, JobTool } from "./jobsApi";

function param(p: Partial<JobParamSpec> & Pick<JobParamSpec, "name" | "type">): JobParamSpec {
    return { label: p.name, options_from: null, default: null, applies_to: null, ...p };
}

const UMAP: JobTool = {
    id: "umap",
    name: "UMAP",
    description: "",
    input_shape: "matrix",
    output: { shape: "column", datasource_param: "datasource", columns_param: "output_name" },
    params: [
        param({ name: "datasource", type: "datasource" }),
        param({ name: "layer", type: "subgroup", options_from: "datasource", default: "gs" }),
        param({ name: "output_name", type: "text", default: "UMAP" }),
        param({ name: "n_neighbors", type: "int", default: 15 }),
        param({ name: "min_dist", type: "float", default: 0.5 }),
    ],
};

const CONCAT: JobTool = {
    ...UMAP,
    id: "concat_columns",
    params: [
        param({ name: "datasource", type: "datasource" }),
        param({ name: "column_a", type: "column", options_from: "datasource" }),
        param({ name: "output_name", type: "text" }),
    ],
};

const noRows = async () => new Uint32Array();

const SOURCES: JobsDataSource[] = [
    {
        name: "cells",
        columns: [{ field: "sample", name: "Sample" }],
        subgroups: ["gs", "rna"],
        size: 3,
        filterSize: 3,
        getFilteredIndices: noRows,
    },
    {
        name: "genes",
        columns: [{ field: "name", name: "Name" }],
        subgroups: [],
        size: 2,
        filterSize: 2,
        getFilteredIndices: noRows,
    },
];

const UMAP_ON_SELECTION: JobTool = {
    ...UMAP,
    params: [
        UMAP.params[0],
        param({ name: "selection", type: "selection", label: "Cells", options_from: "datasource" }),
        ...UMAP.params.slice(1),
    ],
};

describe("paramOptions", () => {
    test("a datasource param offers every datasource", () => {
        expect(paramOptions(UMAP.params[0], {}, SOURCES)).toEqual([
            { value: "cells", label: "cells" },
            { value: "genes", label: "genes" },
        ]);
    });

    test("a subgroup param offers the matrices of the datasource it is bound to", () => {
        expect(paramOptions(UMAP.params[1], { datasource: "cells" }, SOURCES)).toEqual([
            { value: "gs", label: "gs" },
            { value: "rna", label: "rna" },
        ]);
    });

    test("a column param offers the columns of the datasource it is bound to, by display name", () => {
        expect(paramOptions(CONCAT.params[1], { datasource: "genes" }, SOURCES)).toEqual([
            { value: "name", label: "Name" },
        ]);
    });

    test("free-entry params have no options", () => {
        expect(paramOptions(UMAP.params[2], {}, SOURCES)).toBeNull();
        expect(paramOptions(UMAP.params[3], {}, SOURCES)).toBeNull();
    });
});

describe("initialValues", () => {
    test("picks the first datasource, the default matrix and the default scalars", () => {
        expect(initialValues(UMAP, SOURCES)).toEqual({
            datasource: "cells",
            layer: "gs",
            output_name: "UMAP",
            n_neighbors: "15",
            min_dist: "0.5",
        });
    });

    test("falls back to the first matrix when the default is not on the datasource", () => {
        const sources = [{ ...SOURCES[0], subgroups: ["rna"] }];
        expect(initialValues(UMAP, sources).layer).toBe("rna");
    });

    test("leaves columns unchosen", () => {
        expect(initialValues(CONCAT, SOURCES).column_a).toBe("");
    });
});

describe("setParam", () => {
    test("changing the datasource resets the params bound to it", () => {
        const values = { ...initialValues(CONCAT, SOURCES), column_a: "sample" };

        const next = setParam(CONCAT, values, "datasource", "genes", SOURCES);

        expect(next.datasource).toBe("genes");
        expect(next.column_a).toBe("");
    });

    test("changing any other param leaves the rest alone", () => {
        const values = initialValues(UMAP, SOURCES);

        expect(setParam(UMAP, values, "output_name", "U", SOURCES)).toEqual({ ...values, output_name: "U" });
    });
});

describe("toSubmitParams", () => {
    test("sends numbers for numeric params and strings for the rest", () => {
        expect(toSubmitParams(UMAP, initialValues(UMAP, SOURCES))).toEqual({
            datasource: "cells",
            layer: "gs",
            output_name: "UMAP",
            n_neighbors: 15,
            min_dist: 0.5,
        });
    });

    test("leaves out an empty numeric param so the worker uses its default", () => {
        const params = toSubmitParams(UMAP, { ...initialValues(UMAP, SOURCES), n_neighbors: "" });
        expect(params).not.toHaveProperty("n_neighbors");
    });

    test("rejects a picker param left empty", () => {
        expect(() => toSubmitParams(CONCAT, initialValues(CONCAT, SOURCES))).toThrow("Choose a value for column_a");
    });

    test("rejects a non-whole number for an int param", () => {
        expect(() => toSubmitParams(UMAP, { ...initialValues(UMAP, SOURCES), n_neighbors: "2.5" })).toThrow(
            "n_neighbors must be a whole number",
        );
    });

    test("rejects text for a float param", () => {
        expect(() => toSubmitParams(UMAP, { ...initialValues(UMAP, SOURCES), min_dist: "abc" })).toThrow(
            "min_dist must be a number",
        );
    });
});

describe("selection", () => {
    test("starts unnamed", () => {
        expect(initialValues(UMAP_ON_SELECTION, SOURCES).selection).toBe("");
    });

    test("the default output name follows the selection's name", () => {
        let values = initialValues(UMAP_ON_SELECTION, SOURCES);

        values = setParam(UMAP_ON_SELECTION, values, "selection", "T cells", SOURCES);
        expect(values.output_name).toBe("UMAP_T_cells");

        values = setParam(UMAP_ON_SELECTION, values, "selection", "CD4+ T cells", SOURCES);
        expect(values.output_name).toBe("UMAP_CD4_T_cells");

        values = setParam(UMAP_ON_SELECTION, values, "selection", "", SOURCES);
        expect(values.output_name).toBe("UMAP");
    });

    test("an output name the user typed stays when the selection is renamed", () => {
        let values = initialValues(UMAP_ON_SELECTION, SOURCES);
        values = setParam(UMAP_ON_SELECTION, values, "output_name", "MyUMAP", SOURCES);

        values = setParam(UMAP_ON_SELECTION, values, "selection", "T cells", SOURCES);

        expect(values.output_name).toBe("MyUMAP");
    });

    test("toSubmitParams sends the trimmed name and leaves an empty one out", () => {
        const values = initialValues(UMAP_ON_SELECTION, SOURCES);

        expect(toSubmitParams(UMAP_ON_SELECTION, { ...values, selection: " T cells " }).selection).toBe("T cells");
        expect(toSubmitParams(UMAP_ON_SELECTION, values)).not.toHaveProperty("selection");
    });
});

describe("toJobsDataSource", () => {
    test("keeps real columns, collects the subgroup keys of every rows_as_columns link and reads the filter", async () => {
        const ds = {
            name: "cells",
            dataStore: {
                columns: [
                    { field: "sample", name: "Sample" },
                    { field: "gs|CD4|12", name: "CD4", subgroup: "gs" },
                ],
                links: {
                    genes: { rows_as_columns: { subgroups: { gs: {}, rna: {} } } },
                    other: {},
                },
                size: 3,
                filterSize: 2,
                getFilteredIndices: async () => new Uint32Array([0, 2]),
            },
        };

        const source = toJobsDataSource(ds);

        expect(source).toMatchObject({
            name: "cells",
            columns: [{ field: "sample", name: "Sample" }],
            subgroups: ["gs", "rna"],
            size: 3,
            filterSize: 2,
        });
        expect(Array.from(await source.getFilteredIndices())).toEqual([0, 2]);
    });
});
