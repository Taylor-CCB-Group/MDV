import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "../../../utilities/Dialog";
import ZarrMetadataDialogComponent from "./ZarrMetadataDialog";

class ZarrMetadataDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

    constructor() {
        super(
            {
                title: "Zarr Dataset Metadata Viewer",
                width: 600,
                height: 400,
            },
            null,
        );
        this.outer.classList.add("zarrMetadataDialog");
        if (!this.dialog) throw new Error("no dialog??");
        
        this.root = createMdvPortal(
            <ZarrMetadataDialogComponent
                onClose={() => this.close()}
                onResize={(width: number, height: number) => this.resizeDialog(width, height)}
            />,
            this.dialog,
        );
        
        if (this.dialog.parentElement) {
            this.dialog.parentElement.style.display = "none";
        }
    }

    resizeDialog(width: number, height: number): void {
        if (this.outer) {
            this.outer.style.width = `${width}px`;
            this.outer.style.height = `${height}px`;
        }
    }

    close(): void {
        super.close();
        if (this.root) {
            this.root.unmount();
        }
    }
}

BaseDialog.experiment["ZarrMetadataDialogReact"] = ZarrMetadataDialogReact;

export default ZarrMetadataDialogReact;