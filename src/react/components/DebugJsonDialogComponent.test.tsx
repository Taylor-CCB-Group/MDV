import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import type { ReactNode } from "react";
import { beforeEach, describe, expect, test, vi } from "vitest";
import { browserSettings } from "@/lib/browserSettings";
import ValidationFindingsStore from "@/lib/ValidationFindingsStore";
import DebugJsonDialogComponent from "./DebugJsonDialogComponent";

type MockChartManager = {
    validationFindings?: ValidationFindingsStore;
};

const mocks = vi.hoisted<{
    chartManager: MockChartManager;
    jsonViewSources: unknown[];
}>(() => ({
    chartManager: {},
    jsonViewSources: [],
}));

vi.mock("../hooks", () => ({
    useChartManager: () => mocks.chartManager,
}));

vi.mock("./CustomTooltip", () => ({
    default: ({ children }: { children: ReactNode }) => <>{children}</>,
}));

vi.mock("react18-json-view", () => ({
    default: ({ src }: { src: unknown }) => {
        mocks.jsonViewSources.push(src);
        return <pre data-testid="json-view">{JSON.stringify(src)}</pre>;
    },
}));

function buildStoreWithFindings() {
    const store = new ValidationFindingsStore();
    store.addChartFinding("scatter_plot v2", {
        issues: [
            {
                path: "param.0",
                message: "Expected numeric column",
                code: "invalid_type",
            },
        ],
        rawConfig: {
            type: "scatter_plot",
            param: ["cell_type"],
        },
        details: {
            chartId: "chart-1",
        },
    });
    store.addDatasourceFinding("cells", {
        issues: [
            {
                path: "columns.2.datatype",
                message: "Required",
            },
        ],
        rawConfig: {
            name: "cells",
        },
    });
    return store;
}

describe("DebugJsonDialogComponent", () => {
    beforeEach(() => {
        mocks.jsonViewSources.length = 0;
        mocks.chartManager.validationFindings = new ValidationFindingsStore();
    });

    test("renders grouped validation issue details outside the generic JSON view", () => {
        mocks.chartManager.validationFindings = buildStoreWithFindings();

        render(
            <DebugJsonDialogComponent
                json={{
                    chartTypes: [],
                    validationFindings: { old: "payload should not be supplied here" },
                }}
                showValidationSection
            />,
        );

        expect(screen.getByText("Validation issues")).toBeTruthy();
        expect(screen.getByText("Total:")).toBeTruthy();
        expect(screen.getByText("Charts (1)")).toBeTruthy();
        expect(screen.getByText("Datasources (1)")).toBeTruthy();
        expect(screen.getByText("scatter_plot v2")).toBeTruthy();
        expect(screen.getByText("cells")).toBeTruthy();
        expect(screen.getByText("param.0 (invalid_type)")).toBeTruthy();
        expect(screen.getByText("Expected numeric column")).toBeTruthy();
        expect(screen.getByText("columns.2.datatype")).toBeTruthy();
        expect(screen.getByText("Required")).toBeTruthy();
        expect(screen.getAllByText("Raw config")).toHaveLength(2);
        expect(screen.getByText("Details")).toBeTruthy();

        // each top-level key is its own section; the supplied validationFindings payload is omitted
        expect(mocks.jsonViewSources).toContainEqual([]);
        expect(screen.queryByText("validationFindings")).toBeNull();
        expect(JSON.stringify(mocks.jsonViewSources)).not.toContain("payload should not be supplied here");
    });

    test("renders the empty validation state", () => {
        render(<DebugJsonDialogComponent json={{ chartTypes: [] }} showValidationSection />);

        expect(screen.getByText("No validation issues recorded.")).toBeTruthy();
    });

    test("clear button clears all validation findings", () => {
        const store = buildStoreWithFindings();
        mocks.chartManager.validationFindings = store;

        render(<DebugJsonDialogComponent json={{ chartTypes: [] }} showValidationSection />);

        fireEvent.click(screen.getByRole("button", { name: "Clear", hidden: true }));

        expect(store.hasAny).toBe(false);
    });

    describe("initial layout", () => {
        const expandedOf = (title: string) =>
            screen.getByRole("button", { name: title, hidden: true }).getAttribute("aria-expanded");

        test("opens Preferences when there are no validation issues", () => {
            render(<DebugJsonDialogComponent json={{ chartTypes: [] }} showValidationSection />);

            expect(expandedOf("Preferences")).toBe("true");
            expect(expandedOf("Debug and diagnostics")).toBe("false");
        });

        test("opens the validation issues instead when there are some, and keeps that after clearing", () => {
            mocks.chartManager.validationFindings = buildStoreWithFindings();
            render(<DebugJsonDialogComponent json={{ chartTypes: [] }} showValidationSection />);

            expect(expandedOf("Preferences")).toBe("false");
            expect(expandedOf("Debug and diagnostics")).toBe("true");
            expect(expandedOf("Validation issues")).toBe("true");

            fireEvent.click(screen.getByRole("button", { name: "Clear", hidden: true }));
            expect(expandedOf("Debug and diagnostics")).toBe("true");
        });
    });

    describe("search", () => {
        const json = {
            buildInfo: { version: "1.2.3", commit: "abc" },
            datasources: [{ name: "cells" }],
        };
        const search = (value: string) =>
            fireEvent.change(screen.getByLabelText("Search Settings, Sections and Data"), { target: { value } });

        test("splits json into a section per top-level key", () => {
            render(<DebugJsonDialogComponent json={json} />);

            expect(screen.getByText("buildInfo")).toBeTruthy();
            expect(screen.getByText("datasources")).toBeTruthy();
        });

        test("finds a section by title and a setting by label", async () => {
            render(<DebugJsonDialogComponent json={json} />);

            search("convert legacy");
            // the search is debounced
            await waitFor(() => expect(screen.queryByText("buildInfo")).toBeNull());
            expect(screen.getByText("Preferences")).toBeTruthy();
            expect(screen.queryByText("Fault Injection")).toBeNull();

            search("buildinfo");
            await waitFor(() => expect(screen.queryByText("Preferences")).toBeNull());
            expect(screen.getByText("buildInfo")).toBeTruthy();
            expect(screen.queryByText("datasources")).toBeNull();
            expect(screen.queryByText("Preferences")).toBeNull();
        });

        test("reports when nothing matches", async () => {
            render(<DebugJsonDialogComponent json={json} />);

            search("zzz-no-such-thing");
            expect(await screen.findByText("Nothing matches your search.")).toBeTruthy();
        });

        test("groups nested sections and gives each json its own filter", async () => {
            render(
                <DebugJsonDialogComponent
                    json={{ project: { datasources: [{ name: "cells" }, { name: "genes" }], views: { a: 1 } } }}
                    groupKeys={["project"]}
                />,
            );

            expect(screen.getByText("project")).toBeTruthy();
            expect(screen.getByText("datasources")).toBeTruthy();
            expect(screen.getByText("views")).toBeTruthy();
            expect(screen.getAllByLabelText("Filter")).toHaveLength(2);

            fireEvent.change(screen.getAllByLabelText("Filter")[0], { target: { value: "genes" } });
            // assert on what is rendered now, not on the history of renders
            await waitFor(() => {
                const rendered = screen.getAllByTestId("json-view").map((node) => node.textContent ?? "");
                expect(rendered.filter((text) => text.includes("genes") && !text.includes("cells"))).toHaveLength(1);
            });
        });

        test("shows Preferences first and groups debug content under Debug and diagnostics", () => {
            render(<DebugJsonDialogComponent json={json} showValidationSection />);

            const headings = screen.getAllByRole("heading", { level: 3, hidden: true }).map((h) => h.textContent);
            expect(headings[0]).toBe("Preferences");
            expect(headings[1]).toBe("Debug and diagnostics");
            expect(headings.indexOf("buildInfo")).toBeGreaterThan(1);
        });

        test("every browser setting appears in Preferences", () => {
            render(<DebugJsonDialogComponent json={json} />);

            for (const setting of browserSettings) {
                expect(screen.getAllByText(setting.label).length).toBeGreaterThan(0);
            }
        });
    });
});
