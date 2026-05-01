import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useState } from "react";
import { SlickgridReact } from "slickgrid-react";
import FindAndReplaceDialog from "./FindAndReplaceDialog";
import useFindReplace from "../hooks/useFindReplace";
import useSlickGridReact from "../hooks/useSlickGridReact";
import useEditCell from "../hooks/useEditCell";
import ReusableAlertDialog from "@/charts/dialogs/ReusableAlertDialog";
import FeedbackAlertComponent, { type FeedbackAlert, isDebugError } from "./FeedbackAlertComponent";
import AddTableColumnDialog from "./AddTableColumnDialog";
import BulkEditColumnDialog from "./BulkEditColumnDialog";
import ColumnRemovalImpactDialog from "./ColumnRemovalImpactDialog";
import RenameTableColumnDialog from "./RenameTableColumnDialog";

/**
 * Main component for the react table chart
 * 
 * Uses:
 * - useSlickGridReact: Core grid states, refs and events
 * - useFindReplace: Find and replace logic
 * - useEditCell: Cell editing logic
 * 
 * Wrapped in an observer to react to Mobx changes
 */
// todo: add integration tests using playwright to test the working of this component
const TableChartReactComponent = observer(() => {
    const [alertDialogOpen, setAlertDialogOpen] = useState(false);

    const {
        chartId,
        config,
        dataStore,
        gridRef,
        selectionSourceRef,
        isFindReplaceOpen,
        orderedParamColumns,
        searchColumn,
        orderedParamColumnsRef,
        sortedFilteredIndices,
        sortedFilteredIndicesRef,
        options,
        columnDefs,
        handleGridCreated,
        isEditMode,
        isColumnEditable,
        onDialogClose,
        feedbackAlert,
        setFeedbackAlert,
        isAddColumnDialogOpen,
        cloneableColumns,
        addColumnDefaultPosition,
        closeAddColumnDialog,
        handleAddColumn,
        isBulkEditDialogOpen,
        bulkEditColumn,
        closeBulkEditDialog,
        handleBulkEdit,
        renameColumnState,
        closeRenameColumnDialog,
        handleRenameColumn,
        pendingColumnRemoval,
        closeColumnRemovalDialog,
        confirmColumnRemoval,
        openColumnRemovalView,
    } = useSlickGridReact();

    const handleFeedbackAlert = useCallback((alert: FeedbackAlert) => {
        setFeedbackAlert(alert);
    }, [setFeedbackAlert]);

    useEffect(() => {
        setAlertDialogOpen(Boolean(feedbackAlert));
    }, [feedbackAlert]);

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
    } = useFindReplace(
        orderedParamColumns,
        sortedFilteredIndices,
        dataStore,
        searchColumn,
        config,
        gridRef,
        selectionSourceRef,
        handleFeedbackAlert,
        isEditMode,
    );

    const { handleBeforeEditCell, handleCellChange } = useEditCell(
        orderedParamColumnsRef,
        sortedFilteredIndicesRef,
        dataStore,
        gridRef,
        handleFeedbackAlert,
        isEditMode,
    );

    const onClose = useCallback(() => {
        onDialogClose();
        onReset();
    }, [onReset, onDialogClose]);

    const onFeedbackDialogClose = useCallback(() => {
        setAlertDialogOpen(false);
        setFeedbackAlert(null);
    }, [setFeedbackAlert]);

    return (
        <div
            id={`react-grid-${chartId}`}
            data-testid="slickgrid-react-container"
            className="absolute w-[100%] h-[100%] slickgrid-react-container"
        >
            <SlickgridReact
                gridId={`table-${chartId}`}
                columns={columnDefs}
                dataset={[]}
                options={options}
                onReactGridCreated={handleGridCreated}
                onBeforeEditCell={handleBeforeEditCell}
                onCellChange={handleCellChange}
            />

            <div data-testid="find-replace-dialog-wrapper">
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
                    isColumnEditable={isColumnEditable && isEditMode}
                />
            </div>

            <AddTableColumnDialog
                open={isAddColumnDialogOpen}
                cloneableColumns={cloneableColumns}
                defaultPosition={addColumnDefaultPosition}
                onClose={closeAddColumnDialog}
                onSubmit={handleAddColumn}
            />

            <BulkEditColumnDialog
                open={isBulkEditDialogOpen}
                columnName={bulkEditColumn}
                onClose={closeBulkEditDialog}
                onSubmit={handleBulkEdit}
            />

            <RenameTableColumnDialog
                open={Boolean(renameColumnState)}
                columnField={renameColumnState?.field ?? null}
                initialName={renameColumnState?.initialName ?? ""}
                onClose={closeRenameColumnDialog}
                onSubmit={handleRenameColumn}
            />

            <ColumnRemovalImpactDialog
                impact={pendingColumnRemoval?.impact ?? null}
                onClose={closeColumnRemovalDialog}
                onConfirm={confirmColumnRemoval}
                onOpenView={openColumnRemovalView}
                open={Boolean(pendingColumnRemoval)}
            />

            {feedbackAlert && (
                <div data-testid="feedback-alert-dialog">
                    <ReusableAlertDialog
                        open={alertDialogOpen}
                        handleClose={onFeedbackDialogClose}
                        component={<FeedbackAlertComponent feedbackAlert={feedbackAlert} />}
                        isAlertErrorComponent={
                            !isDebugError(feedbackAlert)
                        }
                    />
                </div>
            )}
        </div>
    );
});

export default TableChartReactComponent;
