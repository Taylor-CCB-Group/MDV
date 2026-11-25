import { describe, expect, test } from "vitest";
import { getDataSourceDisplayName, getDataSourceNameProblem } from "@/datastore/dataSourceDisplayName";

/** shaped like a DataStore: the identifier, and the datasources.json entry it was built from */
function store(name: string, display_name?: unknown) {
    return { name, config: display_name === undefined ? { name } : { name, display_name } };
}

describe("getDataSourceDisplayName", () => {
    test.each([
        ["no display name", store("upload.csv"), "upload.csv"],
        ["a display name", store("upload.csv", "Samples"), "Samples"],
        ["surrounding whitespace", store("upload.csv", "  Samples "), "Samples"],
        ["an empty display name", store("upload.csv", ""), "upload.csv"],
        ["a whitespace-only display name", store("upload.csv", "   "), "upload.csv"],
        ["a null display name", store("upload.csv", null), "upload.csv"],
        ["a display name that is not text", store("upload.csv", 42), "upload.csv"],
    ])("with %s", (_, dataStore, shown) => {
        expect(getDataSourceDisplayName(dataStore)).toBe(shown);
    });

    test("reads a datasources.json entry as well as a DataStore", () => {
        expect(getDataSourceDisplayName({ name: "upload.csv", display_name: "Samples" })).toBe("Samples");
        expect(getDataSourceDisplayName({ name: "upload.csv" })).toBe("upload.csv");
    });
});

describe("getDataSourceNameProblem", () => {
    const upload = store("upload.csv");
    const cells = store("cells");
    const all = [cells, upload];

    test.each(["", "   ", "\t\n"])("requires a name: %j", (name) => {
        expect(getDataSourceNameProblem(name, upload, all)).toBe("Name is required");
    });

    test.each(["cells", "CELLS", "  Cells "])("rejects another datasource's name in any case: %j", (name) => {
        expect(getDataSourceNameProblem(name, upload, all)).toBe(`A datasource called '${name.trim()}' already exists`);
    });

    test("rejects another datasource's display name in any case", () => {
        const renamed = [store("cells", "Samples"), upload];
        expect(getDataSourceNameProblem("samples", upload, renamed)).toBe(
            "A datasource called 'samples' already exists",
        );
    });

    test("keeps a renamed datasource's identifier reserved", () => {
        const renamed = [store("cells", "Cell table"), upload];
        expect(getDataSourceNameProblem("cells", upload, renamed)).not.toBeNull();
    });

    test.each(["Samples", "upload.csv", "UPLOAD.CSV", "Upload.csv"])("allows %j, including its own name", (name) => {
        expect(getDataSourceNameProblem(name, upload, all)).toBeNull();
    });

    test("allows changing the case of its own display name", () => {
        const renamed = [cells, store("upload.csv", "Samples")];
        expect(getDataSourceNameProblem("SAMPLES", renamed[1], renamed)).toBeNull();
    });

    test("allows a name another datasource is no longer shown as", () => {
        const renamed = [store("cells", "Cell table"), upload];
        expect(getDataSourceNameProblem("Samples", upload, renamed)).toBeNull();
    });
});
