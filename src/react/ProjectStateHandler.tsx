import DebugErrorComponent, { type DebugErrorComponentProps } from "@/charts/dialogs/DebugErrorComponent";
import AlertErrorComponent from "@/charts/dialogs/AlertErrorComponent";
import { getPostData } from "@/dataloaders/DataLoaderUtil";
import { useCallback, useEffect, useState } from "react";
import ReusableAlertDialog from "@/charts/dialogs/ReusableAlertDialog";

export type ProjectStateHandlerType = {
    root: string;
    data: any;
    staticFolder: boolean;
    permission: boolean;
};

const ProjectStateHandler = ({ root, data, staticFolder, permission }: ProjectStateHandlerType) => {
    const [error, setError] = useState<DebugErrorComponentProps['error'] | null>(null);
    const [errorDialogOpen, setErrorDialogOpen] = useState(false);
    const [confirmSave, setConfirmSave] = useState(false);
    const [permissionDenied, setPermissionDenied] = useState(false);
    const cm = window.mdv.chartManager;

    const saveState = useCallback(async () => {
        if (staticFolder) return;
        try {
            const resp = await getPostData(`${root}/save_state`, data);
            if (resp.success) {
                cm.createInfoAlert("State saved", { duration: 2000 });
                cm.setAllColumnsClean();
            } else {
                throw new Error("An error occurred while saving state. Please try again later.");
            }
        } catch (err) {
            setError(err instanceof Error ? {
                message: err.message,
                stack: err?.stack,
            } : {
                message: "An error occurred while saving state. Please try again later.",
                stack: `${err}`
            });
            setErrorDialogOpen(true);
        }
    }, [cm, staticFolder, data, root]);

    useEffect(() => {
        setError(null);
        setErrorDialogOpen(false);
        setConfirmSave(false);
        setPermissionDenied(false);

        // check for permission
        if (!permission) {
            setPermissionDenied(true);
            return;
        }

        if (data.chartErrors && !Array.isArray(data.chartErrors)) {
            setConfirmSave(true);
            return;
        }

        saveState();
    }, [data, permission, saveState]);

    const handleConfirmProceed = () => {
        setConfirmSave(false);
        saveState();
    };

    return (
        <>
            <ReusableAlertDialog
                open={permissionDenied}
                handleClose={() => setPermissionDenied(false)}
                component={
                    <AlertErrorComponent title="Forbidden" message="You do not have permission to save this project." />
                }
            />

            <ReusableAlertDialog
                open={confirmSave}
                handleClose={() => setConfirmSave(false)}
                component={
                    <AlertErrorComponent
                        title="Proceed with Errors"
                        message={"There are currently errors in your current view. Would you like to proceed anyway?"}
                    />
                }
                isAlertErrorComponent
                isConfirmButton
                confirmText="Proceed"
                onConfirmClick={handleConfirmProceed}
            />

            <ReusableAlertDialog
                open={errorDialogOpen}
                handleClose={() => setErrorDialogOpen(false)}
                component={<DebugErrorComponent error={{ message: error?.message as string, stack: error?.stack }} />}
            />
        </>
    );
};

const ProjectStateHandlerWrapper = ({ root, data, staticFolder, permission }: ProjectStateHandlerType) => {
    return <ProjectStateHandler root={root} data={data} staticFolder={staticFolder} permission={permission} />;
};

export default ProjectStateHandlerWrapper;
