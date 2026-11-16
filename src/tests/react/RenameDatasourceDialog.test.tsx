import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import DataStore from "@/datastore/DataStore";

vi.mock("@/dataloaders/DataLoaderUtil", () => ({ getProjectRoot: () => "/project/p1" }));

import {
    RenameDatasourceDialogComponent,
    RenameDatasourceFlow,
    renameDatasource,
} from "@/charts/dialogs/RenameDatasourceDialog";

function createDataStore(name: string, displayName?: string): DataStore {
    return Object.assign(Object.create(DataStore.prototype), {
        name,
        config: { name, display_name: displayName },
    });
}

function jsonResponse(status: number, body: unknown) {
    return new Response(JSON.stringify(body), { status, headers: { "Content-Type": "application/json" } });
}

const fetchMock = vi.fn<typeof fetch>();
const upload = createDataStore("upload.csv");
const cells = createDataStore("cells", "Cell table");
const onClose = vi.fn();
const onRenamed = vi.fn();

function renderDialog(dataStore = upload) {
    render(
        <RenameDatasourceDialogComponent
            dataStore={dataStore}
            allDataStores={[cells, dataStore]}
            onClose={onClose}
            onRenamed={onRenamed}
        />,
    );
    return {
        input: screen.getByLabelText<HTMLInputElement>("Datasource Name"),
        renameButton: screen.getByRole<HTMLButtonElement>("button", { name: "Rename" }),
    };
}

beforeEach(() => {
    vi.stubGlobal("fetch", fetchMock);
    fetchMock.mockResolvedValue(jsonResponse(200, { success: true, changed: true, display_name: "Samples" }));
});

afterEach(() => {
    cleanup();
    vi.unstubAllGlobals();
    vi.clearAllMocks();
});

describe("RenameDatasourceDialogComponent", () => {
    test("starts from the name the datasource is shown under, with nothing to rename yet", () => {
        const { input, renameButton } = renderDialog(createDataStore("upload.csv", "Patients"));

        expect(input.value).toBe("Patients");
        expect(renameButton.disabled).toBe(true);
        expect(screen.getByText(/code that use "upload.csv" keep working/)).toBeTruthy();
    });

    test("requires a name", () => {
        const { input, renameButton } = renderDialog();

        fireEvent.change(input, { target: { value: "   " } });

        expect(screen.getByText("Name is required")).toBeTruthy();
        expect(renameButton.disabled).toBe(true);
    });

    test.each(["cells", "Cell Table"])("refuses a name another datasource has or is shown as: %j", (name) => {
        const { input, renameButton } = renderDialog();

        fireEvent.change(input, { target: { value: name } });

        expect(screen.getByText(`A datasource called '${name}' already exists`)).toBeTruthy();
        expect(renameButton.disabled).toBe(true);
    });

    test("treats the current name with extra whitespace as unchanged", () => {
        const { input, renameButton } = renderDialog();

        fireEvent.change(input, { target: { value: "  upload.csv " } });

        expect(renameButton.disabled).toBe(true);
    });

    test("renames by identifier with the trimmed name, then reports it", async () => {
        const { input, renameButton } = renderDialog();

        fireEvent.change(input, { target: { value: "  Samples " } });
        fireEvent.click(renameButton);

        await waitFor(() => expect(onRenamed).toHaveBeenCalledOnce());
        expect(fetchMock).toHaveBeenCalledOnce();
        const [url, init] = fetchMock.mock.calls[0];
        expect(url).toBe("/project/p1/rename_datasource");
        expect(init?.method).toBe("POST");
        expect(JSON.parse(String(init?.body))).toEqual({ datasource: "upload.csv", name: "Samples" });
    });

    test("shows the server's reason when the rename is refused, and lets the user try again", async () => {
        fetchMock.mockResolvedValue(
            jsonResponse(400, { success: false, error: "A datasource called 'X' already exists" }),
        );
        const { input, renameButton } = renderDialog();

        fireEvent.change(input, { target: { value: "X" } });
        fireEvent.click(renameButton);

        expect(await screen.findByText("A datasource called 'X' already exists")).toBeTruthy();
        expect(onRenamed).not.toHaveBeenCalled();
        expect(renameButton.disabled).toBe(false);
    });
});

describe("RenameDatasourceFlow", () => {
    const saveView = vi.fn<() => Promise<boolean>>();

    function renderFlow(hasUnsavedChanges: boolean) {
        render(
            <RenameDatasourceFlow
                dataStore={upload}
                allDataStores={[cells, upload]}
                hasUnsavedChanges={hasUnsavedChanges}
                saveView={saveView}
                onClose={onClose}
                onRenamed={onRenamed}
            />,
        );
    }

    const askedToSave = () => screen.queryByText("Save the current view?");
    const askedForName = () => screen.queryByLabelText("Datasource Name");

    test("goes straight to the new name when nothing is unsaved", () => {
        renderFlow(false);

        expect(askedToSave()).toBeNull();
        expect(askedForName()).toBeTruthy();
    });

    test("saves the view first when asked to, then asks for the new name", async () => {
        saveView.mockResolvedValue(true);
        renderFlow(true);
        expect(askedForName()).toBeNull();

        fireEvent.click(screen.getByRole("button", { name: "Save" }));

        await waitFor(() => expect(askedForName()).toBeTruthy());
        expect(saveView).toHaveBeenCalledOnce();
    });

    test("can go on without saving", () => {
        renderFlow(true);

        fireEvent.click(screen.getByRole("button", { name: "Don't save" }));

        expect(askedForName()).toBeTruthy();
        expect(saveView).not.toHaveBeenCalled();
    });

    test.each([
        ["reports it was not saved", () => saveView.mockResolvedValue(false)],
        ["throws", () => saveView.mockRejectedValue(new Error("offline"))],
    ])("stays put and says so when the save %s", async (_, failSave) => {
        failSave();
        renderFlow(true);

        fireEvent.click(screen.getByRole("button", { name: "Save" }));

        expect(await screen.findByText("The view could not be saved.")).toBeTruthy();
        expect(askedForName()).toBeNull();
        expect(screen.getByRole<HTMLButtonElement>("button", { name: "Don't save" }).disabled).toBe(false);
    });

    test("cancelling closes without saving or renaming", () => {
        renderFlow(true);

        fireEvent.click(screen.getByRole("button", { name: "Cancel" }));

        expect(onClose).toHaveBeenCalledOnce();
        expect(saveView).not.toHaveBeenCalled();
        expect(fetchMock).not.toHaveBeenCalled();
    });
});

describe("renameDatasource", () => {
    test("reports whether the name changed", async () => {
        fetchMock.mockResolvedValueOnce(jsonResponse(200, { success: true, changed: false, display_name: "Samples" }));

        await expect(renameDatasource("upload.csv", "Samples")).resolves.toBe(false);
    });

    test("explains a failure with no JSON body", async () => {
        fetchMock.mockResolvedValueOnce(new Response("<html>oops</html>", { status: 500, statusText: "Server Error" }));

        await expect(renameDatasource("upload.csv", "Samples")).rejects.toThrow("Rename failed (500 Server Error)");
    });

    test("treats a response without success as a failure", async () => {
        fetchMock.mockResolvedValueOnce(jsonResponse(200, { success: false, error: "nope" }));

        await expect(renameDatasource("upload.csv", "Samples")).rejects.toThrow("nope");
    });
});
