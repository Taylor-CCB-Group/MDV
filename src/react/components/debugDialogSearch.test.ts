import { describe, expect, test } from "vitest";
import { filterJSON, filterJsonSection, matchesQuery, splitJsonSections } from "./debugDialogSearch";

describe("debugDialogSearch", () => {
    test("matchesQuery requires every word, case-insensitively, across texts", () => {
        expect(matchesQuery("", "anything")).toBe(true);
        expect(matchesQuery("convert LEGACY", "Convert legacy charts (experimental)")).toBe(true);
        expect(matchesQuery("legacy fault", "Convert legacy charts")).toBe(false);
        expect(matchesQuery("legacy fault", "Convert legacy charts", "fault")).toBe(true);
    });

    test("filterJSON keeps matching subtrees and ignores extra whitespace", () => {
        const obj = { a: { name: "cells" }, b: { name: "genes" } };
        expect(filterJSON(obj, "  genes ")).toEqual({ b: { name: "genes" } });
        expect(filterJSON(obj, "none")).toEqual({});
        expect(filterJSON(obj, "")).toBe(obj);
    });

    test("splitJsonSections gives object values their own section and collects primitives", () => {
        const sections = splitJsonSections({ buildInfo: { v: 1 }, datasources: [], mode: "dev" });
        expect(sections.map((s) => s.title)).toEqual(["General", "buildInfo", "datasources"]);
        expect(sections[0].json).toEqual({ mode: "dev" });
    });

    test("filterJsonSection keeps the whole section on a title match, else matching descendants", () => {
        const section = { title: "datasources", json: [{ name: "cells" }, { name: "genes" }] };
        expect(filterJsonSection(section, "datasource")).toBe(section);
        expect(filterJsonSection(section, "genes")?.json).toEqual({ 1: { name: "genes" } });
        expect(filterJsonSection(section, "nope")).toBeNull();
    });

    test("groupKeys turns a key into a group of nested sections, filtered by child", () => {
        const json = { buildInfo: { v: 1 }, project: { datasources: [{ name: "cells" }], views: { a: 1 } } };
        const sections = splitJsonSections(json, ["project"]);
        const project = sections.find((s) => s.title === "project");
        expect(project?.children?.map((c) => c.title)).toEqual(["datasources", "views"]);

        expect(filterJsonSection(project as never, "project")).toBe(project);
        expect(filterJsonSection(project as never, "views")?.children?.map((c) => c.title)).toEqual(["views"]);
        expect(filterJsonSection(project as never, "cells")?.children?.map((c) => c.title)).toEqual(["datasources"]);
        expect(filterJsonSection(project as never, "nope")).toBeNull();
        // the group's title counts as part of the path of what's inside it
        const both = filterJsonSection(project as never, "project cells");
        expect(both?.children?.map((c) => c.title)).toEqual(["datasources"]);
        expect(filterJsonSection(project as never, "project views")?.children?.map((c) => c.title)).toEqual(["views"]);
    });
});
