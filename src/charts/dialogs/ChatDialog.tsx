import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "../../utilities/Dialog";
import ChatDialogComponent from "./ChatDialogComponent";
import { ProjectProvider } from "../../modules/ProjectContext";


class ChatDialog extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

    constructor() {
        super(
            {
                title: "ChatMDV",
                width: 450,
                height: 295,
            },
            null
        );
        if (this.dialog) {
            this.root = createMdvPortal(
                <ProjectProvider>
                    <ChatDialogComponent />
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

BaseDialog.experiment["ChatDialog"] = ChatDialog;

export default ChatDialog;
