import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "../../utilities/Dialog";
import FileUploadDialogComponent from "./FileUploadDialog";

class FileUploadDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

    constructor() {
        super(
            {
                title: "File Upload",
                width: 450,
                height: 320,
            },
            null
        );
        this.outer.classList.add("fileUploadDialog");
        if (this.dialog) {
            this.root = createMdvPortal(
                <FileUploadDialogComponent
                    onClose={() => this.close()}
                    onResize={(width: number, height: number) => this.resizeDialog(width, height)} // Pass the resize callback
                />,
                this.dialog
            );
        } else {
            console.error("Dialog element not found");
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
