import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "../../utilities/Dialog";
import ChatLogDialogComponent from "./ChatLogDialogComponent";
import { ProjectProvider } from "../../modules/ProjectContext";


/** todo reduce need for repetitive wrapper code like this... */
class ChatLogDialog extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

    constructor() {
        super(
            {
                title: "Chat logs",
                width: 800,
                height: 500,
            },
            null
        );
        if (this.dialog) {
            this.root = createMdvPortal(
                <ProjectProvider>
                    <ChatLogDialogComponent />
                </ProjectProvider>,
                this.dialog
            );
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

// BaseDialog.experiment["ChatLogDialog"] = ChatLogDialog; //we don't need this anymore

export default ChatLogDialog;
