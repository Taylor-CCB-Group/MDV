import { beforeEach, describe, expect, test } from "vitest";
import {
    buildApiUrl,
    buildDashboardUrl,
    buildProjectUrl,
    getAppMountPath,
    getApiRootFromDir,
    getDashboardApiRoot,
    getProjectDirFromLocation,
    isDefaultPreviewApiRoot,
    isProjectPath,
    shouldRenderDashboard,
} from "@/utils/mdvRouting";

function setUrl(path: string) {
    window.history.replaceState({}, "", path);
}

describe("mdvRouting", () => {
    beforeEach(() => {
        setUrl("/");
    });

    test("detects mounted app root from login path", () => {
        setUrl("/mdv/login_dev?callback=1");

        expect(getAppMountPath()).toBe("/mdv/");
        expect(buildDashboardUrl("/")).toBe("/mdv/");
    });

    test("builds mounted project urls for same-origin deployments", () => {
        setUrl("/mdv/");

        expect(buildProjectUrl("174", "/")).toBe("/mdv/project/174");
        expect(buildApiUrl("projects", "/")).toBe("/mdv/projects");
    });

    test("handles mounted project pathname routes directly", () => {
        setUrl("/mdvmount/project/174");

        expect(getAppMountPath()).toBe("/mdvmount/");
        expect(isProjectPath()).toBe(true);
        expect(shouldRenderDashboard()).toBe(false);
        expect(getProjectDirFromLocation()).toBe(
            `${window.location.origin}/mdvmount/project/174`,
        );
        expect(buildApiUrl("projects", "/")).toBe("/mdvmount/projects");
        expect(buildProjectUrl("174", "/")).toBe("/mdvmount/project/174");
    });

    test("keeps mounted shell path when building explicit dir urls", () => {
        setUrl("/mdv/");

        expect(buildDashboardUrl("/prefix/")).toBe("/mdv/");
        expect(buildProjectUrl("174", "/prefix/")).toBe("/mdv/project/174");
    });

    test("preserves explicit dir routing when user supplied dir", () => {
        setUrl("/mdv/?dir=%2Fprefix%2F");

        expect(buildDashboardUrl("/prefix/")).toBe("/mdv/?dir=%2Fprefix%2F");
        expect(buildProjectUrl("174", "/prefix/")).toBe("/mdv/?dir=%2Fprefix%2Fproject%2F174");
    });

    test("detects project routes from pathname and query-dir routes", () => {
        setUrl("/project/174");
        expect(isProjectPath()).toBe(true);
        expect(shouldRenderDashboard()).toBe(false);

        setUrl("/?dir=http://localhost:5055/project/174");
        expect(getDashboardApiRoot()).toBe("http://localhost:5055");
        expect(getApiRootFromDir("http://localhost:5055/project/174")).toBe(
            "http://localhost:5055",
        );
        expect(shouldRenderDashboard()).toBe(false);
    });

    test("renders dashboard for root route", () => {
        setUrl("/");
        expect(isProjectPath()).toBe(false);
        expect(shouldRenderDashboard()).toBe(true);
    });

    test("treats preview api roots as equivalent with or without trailing slash", () => {
        expect(isDefaultPreviewApiRoot("http://localhost:5055")).toBe(true);
        expect(isDefaultPreviewApiRoot("http://localhost:5055/")).toBe(true);
        expect(isDefaultPreviewApiRoot("http://localhost:5056/")).toBe(false);
    });
});
