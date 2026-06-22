import { fireEvent, render, screen } from "@testing-library/react";
import type { ReactNode } from "react";
import { beforeEach, describe, expect, test, vi } from "vitest";
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

vi.mock("../viv_loader_cache", () => ({
    vivLoaderCacheTelemetryObservable: {
        snapshot: { cacheEntries: 3 },
    },
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

        const topLevelJsonSource = mocks.jsonViewSources.at(-1);
        expect(topLevelJsonSource).toMatchObject({
            chartTypes: [],
            vivLoaderCacheTelemetry: { cacheEntries: 3 },
        });
        expect(topLevelJsonSource).not.toHaveProperty("validationFindings");
    });

    test("renders the empty validation state", () => {
        render(
            <DebugJsonDialogComponent
                json={{ chartTypes: [] }}
                showValidationSection
            />,
        );

        expect(screen.getByText("No validation issues recorded.")).toBeTruthy();
    });

    test("clear button clears all validation findings", () => {
        const store = buildStoreWithFindings();
        mocks.chartManager.validationFindings = store;

        render(
            <DebugJsonDialogComponent
                json={{ chartTypes: [] }}
                showValidationSection
            />,
        );

        fireEvent.click(screen.getByRole("button", { name: "Clear" }));

        expect(store.hasAny).toBe(false);
    });
});
