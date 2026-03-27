import { beforeEach, describe, expect, test } from "vitest";
import { getProjectInfoBase } from "@/modules/ProjectContext";

describe("ProjectContext", () => {
    beforeEach(() => {
        window.history.replaceState({}, "", "/mdv/project/174");
        (window as any).mdv = {
            ChartManager: class {} as any,
            chartManager: {
                config: {},
            } as any,
            chartTypes: {},
        };
    });

    test("getProjectInfoBase derives project info without React hooks", () => {
        const info = getProjectInfoBase();

        expect(info.projectName).toBe("174");
        expect(info.root).toBe("http://localhost:3000/mdv/project/174");
        expect(info.mainApiRoute).toBe("/mdv");
        expect(info.projectApiRoute).toBe("/mdv/project/174/");
    });
});
