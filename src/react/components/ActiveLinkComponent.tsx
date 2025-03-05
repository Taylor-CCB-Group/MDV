import { CTypes } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import { RAComponent, RowsAsColsProps } from "./LinksComponent";
import { useHighlightedForeignRows } from "../chartLinkHooks";

const ActiveLinkComponent = observer(<T extends CTypes,>(props: RowsAsColsProps<T>) => {
    const { linkedDs, link } = props;
    const { setSelectedColumn } = props;
    const rowNames = useHighlightedForeignRows().map(r => r.fieldName);
    const maxItems = props.multiple ? 10 : 1;

    //! this was causing an infinite loop, need to fix that!!!
    // useEffect(() => {
    //     //@ts-expect-error probably need isMulti logic here - may be able to refactor that into a hook
    //     setSelectedColumn(new RowsAsColsQuery(link, linkedDs.name, maxItems))
    // }, [link, linkedDs, maxItems, setSelectedColumn]);

    return (
        <div>
            <RAComponent {...props} />
            {link.observableFields[0].fieldName}
        </div>
    );
})

export default ActiveLinkComponent;