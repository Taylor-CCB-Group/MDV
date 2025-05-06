import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "@/utilities/Dialog";
import { useState } from "react";
import { ErrorBoundary } from "react-error-boundary";
import DebugErrorComponent from "./DebugErrorComponent";
import AddViewDialogComponent from "./AddViewDialog";
import SaveViewDialogComponent from "./SaveViewDialog";
import DeleteViewDialogComponent from "./DeleteViewDialog";
import SaveAsViewDialogComponent from "./SaveAsViewDialog";

const Wrapper = (props: {
    onDone: () => void;
    dialogType: DialogType;
    saveDialogAction?: () => void;
    saveDialogContent?: string;
}) => {
    const [open, setOpen] = useState(true);

    const getDialogComponent = () => {
        switch (props.dialogType) {
            case "save":
                return (
                    <SaveViewDialogComponent
                        open={open}
                        onClose={() => {
                            setOpen(false);
                            props.onDone();
                        }}
                        action={props?.saveDialogAction}
                        content={props?.saveDialogContent}
                    />
                );
            case "add":
                return (
                    <AddViewDialogComponent
                        open={open}
                        onClose={() => {
                            setOpen(false);
                            props.onDone();
                        }}
                    />
                );
            case "delete":
                return (
                    <DeleteViewDialogComponent
                        open={open}
                        onClose={() => {
                            setOpen(false);
                            props.onDone();
                        }}
                    />
                );
            case "save_as":
                return (
                    <SaveAsViewDialogComponent
                        open={open}
                        onClose={() => {
                            setOpen(false);
                            props.onDone();
                        }}
                    />
                );
        }
    };

    return (
        <ErrorBoundary
            FallbackComponent={({ error }) => (
                <DebugErrorComponent error={error} title="Unexpected Error: please report to developers." />
            )}
        >
            {getDialogComponent()}
        </ErrorBoundary>
    );
};

type DialogType = "save" | "add" | "delete" | "save_as";

class ViewDialogWrapper extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;
    constructor(dialogType: DialogType, saveDialogAction?: () => void, saveDialogContent?: string) {
        super(
            {
                // title: `Add View'`,
                // width: 300,
                // height: 220,
            },
            null,
        );
        this.root = createMdvPortal(
            <Wrapper
                onDone={() => this.close()}
                dialogType={dialogType}
                saveDialogAction={saveDialogAction}
                saveDialogContent={saveDialogContent}
            />,
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

// mercifully, this isn't necessary
// BaseDialog.experiment["ViewDialogWrapper"] = ViewDialogWrapper;
export default ViewDialogWrapper;
