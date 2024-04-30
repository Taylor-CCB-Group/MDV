import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "../../utilities/Dialog";
import FileUploadDialogComponent from "./FileUploadDialog";
import { ProjectProvider} from "../../modules/ProjectContext";
class FileUploadDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

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
            this.root = createMdvPortal(<ProjectProvider><FileUploadDialogComponent onClose={() => this.close()} /></ProjectProvider>, this.dialog);
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

export default FileUploadDialogReact;