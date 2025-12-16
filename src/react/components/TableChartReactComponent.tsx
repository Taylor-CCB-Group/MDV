import { observer } from "mobx-react-lite";
import { useCallback } from "react";
import {
    SlickgridReact,
} from "slickgrid-react";
import FindAndReplaceDialog from "./FindAndReplaceDialog";
import useFindReplace from "../hooks/useFindReplace";
import useSlickGridReact from "../hooks/useSlickGridReact";
import useEditCell from "../hooks/useEditCell";

const TableChartReactComponent = observer(() => {
    const {
        chartId,
        config,
        dataStore,
        gridRef,
        isSelectingRef,
        isFindReplaceOpen,
        orderedParamColumns,
        searchColumn,
        orderedParamColumnsRef,
        sortedIndices,
        sortedIndicesRef,
        options,
        columnDefs,
        handleGridCreated,
        isColumnEditable,
        onDialogClose,
    } = useSlickGridReact();

    const {
        matchCount,
        disableFindPrev,
        disableFindNext,
        handleFind,
        handleFindNext,
        handleFindPrev,
        handleReplace,
        handleReplaceAll,
        onReset,
    } = useFindReplace(orderedParamColumns,
        sortedIndices,
        dataStore,
        searchColumn,
        config,
        gridRef,
        isSelectingRef);

    const {
        handleBeforeEditCell,
        handleCellChange,
    } = useEditCell(
        orderedParamColumnsRef, sortedIndicesRef, dataStore, gridRef
    );

    const onClose = useCallback(() => {
        onDialogClose();
        onReset();
    }, [onReset, onDialogClose]);

    return (
        <div id={`react-grid-${chartId}`} className="absolute w-[100%] h-[100%] slickgrid-react-container">
            <SlickgridReact
                gridId={`table-${chartId}`}
                columns={columnDefs}
                dataset={[]}
                options={options}
                onReactGridCreated={handleGridCreated}
                onBeforeEditCell={handleBeforeEditCell}
                onCellChange={handleCellChange}
            />

            <FindAndReplaceDialog
                open={isFindReplaceOpen}
                onClose={onClose}
                handleFind={handleFind}
                handleFindPrev={handleFindPrev}
                handleFindNext={handleFindNext}
                handleReplace={handleReplace}
                handleReplaceAll={handleReplaceAll}
                foundMatches={matchCount}
                columnName={searchColumn}
                disableFindPrev={disableFindPrev}
                disableFindNext={disableFindNext}
                isColumnEditable={isColumnEditable}
            />
        </div>
    );
});

export default TableChartReactComponent;
