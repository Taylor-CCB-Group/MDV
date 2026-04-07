import {
    AnnotationDialogComponent,
    getAnnotationColumnChoiceState,
    getEditableAnnotationColumnNames,
} from "@/charts/dialogs/AnnotationDialogReact";
import TagModel from "@/table/TagModel";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { beforeEach, describe, expect, test, vi } from "vitest";

vi.mock("@/table/TagModel", () => ({
    default: {
        create: vi.fn(),
    },
}));

function createMockTagModel(
    dataStore: { isFiltered?: () => boolean; size?: number },
    selectionScope: "filtered" | "highlighted",
    selectionCount: number,
) {
    const getAnnotationViewState = () => ({
        tagList: new Set(["a", "b"]),
        tagsInSelection: new Set(["a"]),
        selectionCount,
    });
    const getAnnotationWarningFlags = () => {
        const showNoSelectionWarning = selectionCount === 0;
        let showFilteredScopeCoversWholeTableWarning = false;
        if (selectionScope === "filtered" && selectionCount > 0) {
            if (typeof dataStore.isFiltered === "function") {
                showFilteredScopeCoversWholeTableWarning = !dataStore.isFiltered();
            } else {
                const totalRows = dataStore.size ?? selectionCount;
                showFilteredScopeCoversWholeTableWarning =
                    selectionCount === totalRows;
            }
        }
        return {
            showFilteredScopeCoversWholeTableWarning,
            showNoSelectionWarning,
        };
    };
    return {
        addListener: vi.fn(),
        dispose: vi.fn(),
        entireSelectionHasTag: (tag: string) => tag === "a",
        getAnnotationViewState,
        getAnnotationWarningFlags,
        getMatchingRowIndices: vi.fn(() => [0, 2]),
        getSelectionLength: vi.fn(() => selectionCount),
        getTags: vi.fn(() => new Set(["a", "b"])),
        getTagsInSelection: vi.fn(() => new Set(["a"])),
        removeListener: vi.fn(),
        setTag: vi.fn(),
    };
}

function createMockDataStore() {
    const columns = [
        {
            datatype: "multitext",
            editable: false,
            field: "__tags",
            name: "__tags",
        },
        {
            datatype: "multitext",
            editable: true,
            field: "editable_tags",
            name: "editable_tags",
        },
        {
            datatype: "multitext",
            editable: false,
            field: "readonly_tags",
            name: "readonly_tags",
        },
        {
            datatype: "text",
            editable: true,
            field: "labels",
            name: "labels",
        },
    ];

    return {
        columnIndex: Object.fromEntries(
            columns.map((column) => [column.field, column]),
        ),
        columns,
        dataHighlighted: vi.fn(),
        getHighlightedData: vi.fn(() => []),
        isFiltered: vi.fn(() => false),
        name: "test-store",
        size: 100,
    };
}

describe("AnnotationDialogReact", () => {
    beforeEach(() => {
        vi.clearAllMocks();
        vi.mocked(TagModel.create).mockImplementation(
            async (dataStore, _columnName, selectionScope = "filtered") =>
                createMockTagModel(
                    dataStore as { isFiltered?: () => boolean; size?: number },
                    selectionScope,
                    selectionScope === "highlighted" ? 1 : 2,
                ) as never,
        );
    });

    test("filters the chooser to editable multitext columns", () => {
        const dataStore = createMockDataStore();

        expect(getEditableAnnotationColumnNames(dataStore.columns)).toEqual([
            "editable_tags",
        ]);
        expect(
            getAnnotationColumnChoiceState(dataStore as never, "readonly_tags"),
        ).toBe("readonly");
        expect(
            getAnnotationColumnChoiceState(dataStore as never, "labels"),
        ).toBe("clash");
        expect(
            getAnnotationColumnChoiceState(dataStore as never, "new_annotations"),
        ).toBe("ok");
    });

    test("renders tabbed scopes and shows read-only validation feedback", async () => {
        const dataStore = createMockDataStore();

        render(<AnnotationDialogComponent dataStore={dataStore as never} />);

        expect(screen.getByRole("tab", { name: "Filtered" })).toBeTruthy();
        expect(screen.getByRole("tab", { name: "Highlighted" })).toBeTruthy();

        fireEvent.click(screen.getByText("Annotation Setup"));

        const input = screen.getByLabelText("Column name");
        fireEvent.change(input, { target: { value: "readonly_tags" } });

        await waitFor(() => {
            const columnField = screen.getByLabelText("Column name");
            expect(columnField.getAttribute("aria-invalid")).toBe("true");
            expect(screen.getByText(/read-only/i)).toBeTruthy();
        });

        expect(vi.mocked(TagModel.create)).toHaveBeenCalledWith(
            dataStore,
            "editable_tags",
            "filtered",
        );
        expect(vi.mocked(TagModel.create)).toHaveBeenCalledWith(
            dataStore,
            "editable_tags",
            "highlighted",
        );
    });

    test("defaults to highlighted scope when highlighted rows exist and warns only when filtered would hit all rows", async () => {
        const dataStore = createMockDataStore();
        dataStore.getHighlightedData = vi.fn(() => [0] as never);

        render(<AnnotationDialogComponent dataStore={dataStore as never} />);

        await waitFor(() => {
            expect(screen.getByText(/Annotating 1 .*highlighted row/)).toBeTruthy();
        });

        fireEvent.click(screen.getByRole("tab", { name: "Filtered" }));

        expect(
            screen.getByTestId("annotation-whole-table-filter-warning"),
        ).toBeTruthy();

        dataStore.isFiltered = vi.fn(() => true);
        fireEvent.click(screen.getByRole("tab", { name: "Highlighted" }));
        fireEvent.click(screen.getByRole("tab", { name: "Filtered" }));

        expect(
            screen.queryByTestId("annotation-whole-table-filter-warning"),
        ).toBeNull();
    });

    test("shows a warning when the active scope has no matching rows", async () => {
        const dataStore = createMockDataStore();
        dataStore.getHighlightedData = vi.fn(() => [0] as never);
        vi.mocked(TagModel.create).mockImplementation(
            async (dataStore, _columnName, selectionScope = "filtered") =>
                createMockTagModel(
                    dataStore as { isFiltered?: () => boolean; size?: number },
                    selectionScope,
                    selectionScope === "highlighted" ? 0 : 2,
                ) as never,
        );

        render(<AnnotationDialogComponent dataStore={dataStore as never} />);

        await waitFor(() => {
            expect(screen.getByTestId("annotation-no-rows-warning")).toBeTruthy();
        });

        expect(screen.getByTestId("annotation-no-rows-warning").textContent).toMatch(
            /highlighted/i,
        );
    });

    test("shows tag mutation errors and keeps the draft value when add fails", async () => {
        const dataStore = createMockDataStore();
        const tagModel = createMockTagModel(dataStore, "filtered", 2);
        tagModel.setTag.mockImplementation(() => {
            throw new Error("Tag update failed");
        });
        vi.mocked(TagModel.create).mockResolvedValue(tagModel as never);

        render(<AnnotationDialogComponent dataStore={dataStore as never} />);

        const input = await screen.findByLabelText("Add annotation item");
        fireEvent.change(input, { target: { value: "new-tag" } });
        fireEvent.click(screen.getByRole("button", { name: "Add" }));

        await waitFor(() => {
            expect(screen.getByText("Tag update failed")).toBeTruthy();
        });

        expect(tagModel.setTag).toHaveBeenCalledWith("new-tag", true);
        expect(screen.getByDisplayValue("new-tag")).toBeTruthy();
    });

    test("shows tag mutation errors when toggling an available tag row fails", async () => {
        const dataStore = createMockDataStore();
        const tagModel = createMockTagModel(dataStore, "filtered", 2);
        tagModel.setTag.mockImplementation(() => {
            throw new Error("Row toggle failed");
        });
        vi.mocked(TagModel.create).mockResolvedValue(tagModel as never);

        render(<AnnotationDialogComponent dataStore={dataStore as never} />);

        fireEvent.click(await screen.findByText("b"));

        await waitFor(() => {
            expect(screen.getByText("Row toggle failed")).toBeTruthy();
        });

        expect(tagModel.setTag).toHaveBeenCalledWith("b", true);
    });
});
