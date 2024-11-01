import type DataStore from "@/datastore/DataStore";
import { createMdvPortal } from "@/react/react_utils";
import { BaseDialog } from "@/utilities/Dialog";
import { observer } from "mobx-react-lite";
import type ChartManager from "../ChartManager";
import type { RowsAsColumnsLink } from "../charts";

type WrapperProps = {
    ds: DataStore;
    dsTo: DataStore;
    link: RowsAsColumnsLink;
    cm: ChartManager;
};

// todo some component that will show the currently active rows from the linked data store
// "dsTo" is one of the things that's confusing about the existing stuff... in some sense it's where
// we link "to" in terms of how it's defined in `datasources` - but in terms of data flow, it's the opposite.
// const ActiveRows = 
const Test = () => {
    return <>Oh, happy day!</>
}

const Wrapper = observer(({ ds, dsTo, link, cm }: WrapperProps) => {
    return (
        <div>
            <h2>Add Chart</h2>
            <p>Choose a chart type to add to the data store '{ds.name}'</p>
            <Test />
        </div>
    );
});


class AddColumnsFromRowsDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;
    constructor(ds: DataStore, dsTo: DataStore, link: RowsAsColumnsLink, cm: ChartManager) {
        super(
            {
                title: `Add Chart in '${ds.name}'`,
                width: 600,
                height: 420,
            },
            null,
        );
        const props = { ds, dsTo, link, cm };
        this.root = createMdvPortal(
            <Wrapper {...props} />,
            this.dialog,
            this,
        );
    }
    close(): void {
        super.close();
        this.root.unmount();
    }
}
// https://github.com/Taylor-CCB-Group/MDV/discussions/44
BaseDialog.experiment["AddColumnsFromRows"] = AddColumnsFromRowsDialogReact;
// export default AddColumnsFromRowsDialogReact;
