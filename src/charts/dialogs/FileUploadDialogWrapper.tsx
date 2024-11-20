import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "../../utilities/Dialog";
import FileUploadDialogComponent from "./FileUploadDialog";
import { VivProvider, createVivStores } from '../../react/components/avivatorish/state';
// import { Dialog } from "@mui/material";

class FileUploadDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

    constructor() {
        super(
            {
                title: "File Upload",
                width: 455,
                height: 300,
            },
            null,
        );
        this.outer.classList.add("fileUploadDialog");
        if (!this.dialog) throw new Error("no dialog??");
        const vivStores = createVivStores();
        this.root = createMdvPortal(
            <VivProvider vivStores={vivStores}>
                {/* <Dialog open={true} className="w-full h-full"> */}
                <FileUploadDialogComponent
                    onClose={() => this.close()}
                    onResize={(width: number, height: number) => {}}
                    // onResize={(width: number, height: number) => this.resizeDialog(width, height)}
                />
                {/* </Dialog> */}
            </VivProvider>,
            this.dialog
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


BaseDialog.experiment["FileUploadDialogReact"] = FileUploadDialogReact;

export default FileUploadDialogReact;