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

function createMockTagModel(selectionCount: number) {
    return {
        addListener: vi.fn(),
        dispose: vi.fn(),
        entireSelectionHasTag: (tag: string) => tag === "a",
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
        name: "test-store",
    };
}

describe("AnnotationDialogReact", () => {
    beforeEach(() => {
        vi.clearAllMocks();
        vi.mocked(TagModel.create).mockImplementation(
            async (_dataStore, _columnName, selectionScope) =>
                createMockTagModel(selectionScope === "highlighted" ? 1 : 2) as never,
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
        expect(screen.queryByText("Apply actions to")).toBeNull();
        expect(
            screen.queryByText(
                "Switch between filtered and highlighted scopes without leaving the dialog.",
            ),
        ).toBeNull();

        fireEvent.click(screen.getByText("Annotation Setup"));

        const input = screen.getByLabelText("Column name");
        fireEvent.change(input, { target: { value: "readonly_tags" } });

        await waitFor(() => {
            expect(
                screen.getByText(
                    "That multitext column is read-only and can't be edited here.",
                ),
            ).toBeTruthy();
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

    test("shows tag mutation errors and keeps the draft value when add fails", async () => {
        const dataStore = createMockDataStore();
        const tagModel = createMockTagModel(2);
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
        const tagModel = createMockTagModel(2);
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
