import { BaseDialog } from "@/utilities/Dialog";
import { createMdvPortal } from "@/react/react_utils";
import type { ManageGateDialogType } from "./ManageGateDialogComponent";
import ManageGateDialogContent from "./ManageGateDialogComponent";

export type ManageGateDialogWrapperProps = Omit<ManageGateDialogType, "open">;

class ManageGateDialogWrapper extends BaseDialog {
    declare root: ReturnType<typeof createMdvPortal>;
    
    constructor(props: ManageGateDialogWrapperProps) {
        const config = {
            title: "Manage Gates",
            width: 500,
            onclose: () => props.onClose(),
        };
        super(config, props);

        this.root = createMdvPortal(
            <ManageGateDialogContent
                {...props}
                open={true}
                onClose={() => {
                    props.onClose();
                    this.close();
                }}
            />,
            this.dialog,
            this,
        );
    }

    close(): void {
        this.root?.unmount();
        super.close();
    }
}

export default ManageGateDialogWrapper;
