import axios from "axios";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import { getProject, getProjectInfoBase } from "@/modules/ProjectContext";

describe("ProjectContext", () => {
    let originalMdv: unknown;

    beforeEach(() => {
        originalMdv = (window as any).mdv;
        window.history.replaceState({}, "", "/mdv/project/174");
        (window as any).mdv = {
            ChartManager: class {} as any,
            chartManager: {
                config: {},
            } as any,
            chartTypes: {},
        };
    });

    afterEach(() => {
        if (originalMdv === undefined) {
            delete (window as any).mdv;
            return;
        }

        (window as any).mdv = originalMdv;
    });

    test("getProjectInfoBase derives project info without React hooks", () => {
        const info = getProjectInfoBase();

        expect(info.projectName).toBe("174");
        expect(info.root).toBe(`${window.location.origin}/mdv/project/174`);
        expect(info.mainApiRoute).toBe("/mdv/");
        expect(info.projectApiRoute).toBe("/mdv/project/174/");
    });

    test("getProject fetches just that project", async () => {
        const project = { id: 7, name: "Lung atlas", lastModified: "", owner: [], writable: true, permissions: null };
        const get = vi.spyOn(axios, "get").mockResolvedValue({ data: project });
        try {
            expect(await getProject(7, "/mdv/")).toEqual(project);
            expect(get).toHaveBeenCalledWith("/mdv/projects/7");
        } finally {
            get.mockRestore();
        }
    });

    test("getProject returns null when the request fails", async () => {
        const error = vi.spyOn(console, "error").mockImplementation(() => {});
        const get = vi.spyOn(axios, "get").mockRejectedValue(new Error("Not found"));
        try {
            expect(await getProject(99, "/mdv/")).toBeNull();
        } finally {
            get.mockRestore();
            error.mockRestore();
        }
    });
});
