import { createRoot } from "react-dom/client";
import { BaseDialog } from "../../utilities/Dialog";
import FileUploadDialogComponent from "./FileUploadDialog";

class FileUploadDialogReact extends BaseDialog {
    root: ReturnType<typeof createRoot>;

    constructor() {
        super(
            {
                title: "File Upload",
                width: 450,
                height: 295,
            },
            null
        );
        this.outer.classList.add("fileUploadDialog");
        if (this.dialog) {
            this.root = createRoot(this.dialog);
            this.root.render(<FileUploadDialogComponent onClose={() => this.close()} />);
        } else {
            console.error("Dialog element not found");
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

export default 42;