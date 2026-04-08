import type DataStore from "../../datastore/DataStore.js";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import { AnnotationDialogComponent } from "./AnnotationDialogComponent";

export {
    AnnotationDialogComponent,
    getAnnotationColumnChoiceState,
    getEditableAnnotationColumnNames,
} from "./AnnotationDialogComponent";

class AnnotationDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

    constructor(dataStore: DataStore) {
        super(
            {
                title: `Annotate '${dataStore.name}'`,
                width: 560,
                height: 740,
            },
            null,
        );
        this.root = createMdvPortal(
            <AnnotationDialogComponent dataStore={dataStore} />,
            this.dialog,
            this,
        );
    }

    close(): void {
        super.close();
        this.root.unmount();
    }
}

BaseDialog.experiment["AnnotationDialogReact"] = AnnotationDialogReact;
export default "AnnotationDialogReact loaded";
