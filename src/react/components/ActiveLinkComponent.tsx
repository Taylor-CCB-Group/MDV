import type { CTypes } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import { RAComponent, type RowsAsColsProps } from "./LinksComponent";
import { useHighlightedForeignRows } from "../chartLinkHooks";

const ActiveLinkComponent = observer(<T extends CTypes, M extends boolean>(props: RowsAsColsProps<T, M>) => {
    const { linkedDs, link } = props;
    const { setSelectedColumn } = props;
    const rowNames = useHighlightedForeignRows().map(r => r.fieldName);
    const maxItems = props.multiple ? 10 : 1;

    return (
        <div>
            <RAComponent {...props} />
            {link.observableFields[0]?.fieldName}
        </div>
    );
})

export default ActiveLinkComponent;