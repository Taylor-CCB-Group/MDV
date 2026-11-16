import { Alert, Button, Dialog, DialogActions, DialogContent, DialogTitle, TextField, Typography } from "@mui/material";
import { useState } from "react";
import { ErrorBoundary } from "react-error-boundary";
import { getProjectRoot } from "@/dataloaders/DataLoaderUtil";
import type DataStore from "@/datastore/DataStore";
import { getDataSourceDisplayName, getDataSourceNameProblem } from "@/datastore/dataSourceDisplayName";
import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "@/utilities/Dialog";
import DebugErrorComponent from "./DebugErrorComponent";

/**
 * Sets the name `datasource` is shown under (see `MDVProject.rename_datasource`), and returns
 * whether it changed. Throws with the server's reason if the rename is refused.
 */
export async function renameDatasource(datasource: string, name: string): Promise<boolean> {
    // not getPostData: that discards the body of a 400, which is where the reason is
    const resp = await fetch(`${getProjectRoot()}/rename_datasource`, {
        method: "POST",
        body: JSON.stringify({ datasource, name }),
        headers: { "Content-Type": "application/json" },
    });
    const data = await resp.json().catch(() => null);
    if (!resp.ok || !data?.success) {
        throw new Error(data?.error ?? `Rename failed (${resp.status} ${resp.statusText})`);
    }
    return Boolean(data.changed);
}

export const RenameDatasourceDialogComponent = (props: {
    dataStore: DataStore;
    allDataStores: DataStore[];
    onClose: () => void;
    onRenamed: () => void;
}) => {
    const { dataStore, allDataStores } = props;
    const currentName = getDataSourceDisplayName(dataStore);
    const [value, setValue] = useState(currentName);
    const [saving, setSaving] = useState(false);
    const [error, setError] = useState<string | null>(null);

    const name = value.trim();
    const nameProblem = getDataSourceNameProblem(value, dataStore, allDataStores);
    const isRenameDisabled = saving || name === currentName || nameProblem !== null;

    const handleSubmit = async (event: HTMLFormElement) => {
        event.preventDefault();
        if (isRenameDisabled) return;
        setSaving(true);
        setError(null);
        try {
            await renameDatasource(dataStore.name, name);
            props.onRenamed();
        } catch (e) {
            setError(e instanceof Error ? e.message : String(e));
            setSaving(false);
        }
    };

    return (
        <Dialog
            open
            onClose={saving ? undefined : props.onClose}
            fullWidth
            maxWidth="xs"
            PaperProps={{
                component: "form",
                onSubmit: handleSubmit,
            }}
        >
            <DialogTitle>Rename Datasource</DialogTitle>
            <DialogContent dividers>
                <TextField
                    label="Datasource Name"
                    variant="outlined"
                    fullWidth
                    margin="normal"
                    autoFocus
                    value={value}
                    disabled={saving}
                    onChange={(e) => setValue(e.target.value)}
                    error={nameProblem !== null}
                    helperText={nameProblem}
                />
                <Typography variant="body2" color="text.secondary">
                    This changes the name shown in MDV. Views, charts and code that use "{dataStore.name}" keep working.
                    The page reloads after renaming.
                </Typography>
                {error && (
                    <Alert severity="error" sx={{ mt: 2 }}>
                        {error}
                    </Alert>
                )}
            </DialogContent>
            <DialogActions>
                <Button type="submit" disabled={isRenameDisabled}>
                    {saving ? "Renaming..." : "Rename"}
                </Button>
                <Button color="error" onClick={props.onClose} disabled={saving}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

/**
 * Asked before renaming when the view has unsaved changes, as the rename reloads the page.
 * `saveView` resolves to whether the view was saved.
 */
export const SaveViewFirstDialogComponent = (props: {
    saveView: () => Promise<boolean>;
    onContinue: () => void;
    onClose: () => void;
}) => {
    const [saving, setSaving] = useState(false);
    const [failed, setFailed] = useState(false);

    const handleSave = async () => {
        setSaving(true);
        setFailed(false);
        const saved = await props.saveView().catch(() => false);
        setSaving(false);
        if (saved) {
            props.onContinue();
        } else {
            setFailed(true);
        }
    };

    return (
        <Dialog open onClose={saving ? undefined : props.onClose} fullWidth maxWidth="xs">
            <DialogTitle>Save the current view?</DialogTitle>
            <DialogContent dividers>
                <Typography variant="body1">
                    Renaming a datasource reloads the page, so unsaved changes to this view will be lost.
                </Typography>
                {failed && (
                    <Alert severity="error" sx={{ mt: 2 }}>
                        The view could not be saved.
                    </Alert>
                )}
            </DialogContent>
            <DialogActions>
                <Button onClick={handleSave} disabled={saving} color="success">
                    {saving ? "Saving..." : "Save"}
                </Button>
                <Button onClick={props.onContinue} disabled={saving}>
                    Don't save
                </Button>
                <Button color="error" onClick={props.onClose} disabled={saving}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

type ViewSaver = {
    viewManager: { saveView(): Promise<void>; hasUnsavedChanges(): boolean };
    getState(): { chartErrors?: unknown[] } | undefined;
};

/** Saves the current view, and resolves to whether it was saved. */
export async function saveCurrentView(cm: ViewSaver): Promise<boolean> {
    // saveView reports nothing and swallows its own errors, so check what it left behind
    await cm.viewManager.saveView();
    // saveView refuses to save a view with chart errors, and hasUnsavedChanges reports
    // false for one, so it cannot say whether the save happened
    if (cm.getState()?.chartErrors?.length) return false;
    return !cm.viewManager.hasUnsavedChanges();
}

/** Offers to save the view if it has unsaved changes, then asks for the new name. */
export const RenameDatasourceFlow = (props: {
    dataStore: DataStore;
    allDataStores: DataStore[];
    hasUnsavedChanges: boolean;
    saveView: () => Promise<boolean>;
    onClose: () => void;
    onRenamed: () => void;
}) => {
    const [step, setStep] = useState<"save" | "rename">(props.hasUnsavedChanges ? "save" : "rename");
    if (step === "save") {
        return (
            <SaveViewFirstDialogComponent
                saveView={props.saveView}
                onContinue={() => setStep("rename")}
                onClose={props.onClose}
            />
        );
    }
    return (
        <RenameDatasourceDialogComponent
            dataStore={props.dataStore}
            allDataStores={props.allDataStores}
            onClose={props.onClose}
            onRenamed={props.onRenamed}
        />
    );
};

class RenameDatasourceDialog extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;
    constructor(dataStore: DataStore) {
        super({}, null);
        const cm = window.mdv.chartManager;
        const allDataStores = cm.dataSources.map((ds) => ds.dataStore);
        const saveView = () => saveCurrentView(cm);
        this.root = createMdvPortal(
            <ErrorBoundary
                FallbackComponent={({ error }) => (
                    <DebugErrorComponent error={error} title="Unexpected Error: please report to developers." />
                )}
            >
                <RenameDatasourceFlow
                    dataStore={dataStore}
                    allDataStores={allDataStores}
                    hasUnsavedChanges={cm.viewManager.hasUnsavedChanges()}
                    saveView={saveView}
                    onClose={() => this.close()}
                    // the new name is read when the page loads
                    onRenamed={() => window.location.reload()}
                />
            </ErrorBoundary>,
            this.dialog,
            this,
        );
        this.outer.style.display = "none";
    }
    close(): void {
        super.close();
        this.root.unmount();
    }
}

export default RenameDatasourceDialog;
