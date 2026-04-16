import type { ColumnRemovalImpact } from "@/charts/columnRemovalUtils";
import ColumnRemovalImpactDialog from "@/react/components/ColumnRemovalImpactDialog";
import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "@/utilities/Dialog";
import { useState } from "react";

type WrapperProps = {
    columnName: string;
    impact: ColumnRemovalImpact;
    onConfirm: () => void;
    onDone: () => void;
};

function Wrapper({ columnName, impact, onConfirm, onDone }: WrapperProps) {
    const [open, setOpen] = useState(true);

    const handleClose = () => {
        setOpen(false);
        onDone();
    };

    const handleConfirm = () => {
        onConfirm();
        handleClose();
    };

    return (
        <ColumnRemovalImpactDialog
            open={open}
            columnName={columnName}
            impact={impact}
            onClose={handleClose}
            onConfirm={handleConfirm}
        />
    );
}

export default class ColumnRemovalImpactDialogWrapper extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

    constructor(columnName: string, impact: ColumnRemovalImpact, onConfirm: () => void) {
        super({}, null);
        this.root = createMdvPortal(
            <Wrapper
                columnName={columnName}
                impact={impact}
                onConfirm={onConfirm}
                onDone={() => this.close()}
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
