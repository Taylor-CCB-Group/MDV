import { observer } from "mobx-react-lite";
import { useCallback, useState } from "react";
import { SlickgridReact } from "slickgrid-react";
import FindAndReplaceDialog from "./FindAndReplaceDialog";
import useFindReplace from "../hooks/useFindReplace";
import useSlickGridReact from "../hooks/useSlickGridReact";
import useEditCell from "../hooks/useEditCell";
import ReusableAlertDialog from "@/charts/dialogs/ReusableAlertDialog";
import FeedbackAlertComponent, { type FeedbackAlert, isDebugError } from "./FeedbackAlertComponent";

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
    const [feedbackAlert, setFeedbackAlert] = useState<FeedbackAlert>(null);

    const handleFeedbackAlert = useCallback((alert: FeedbackAlert) => {
        setFeedbackAlert(alert);
        if (alert) {
            setAlertDialogOpen(true);
        }
    }, []);

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
    } = useFindReplace(
        orderedParamColumns,
        sortedFilteredIndices,
        dataStore,
        searchColumn,
        config,
        gridRef,
        selectionSourceRef,
        handleFeedbackAlert,
    );

    const { handleBeforeEditCell, handleCellChange } = useEditCell(
        orderedParamColumnsRef,
        sortedFilteredIndicesRef,
        dataStore,
        gridRef,
        handleFeedbackAlert,
    );

    const onClose = useCallback(() => {
        onDialogClose();
        onReset();
    }, [onReset, onDialogClose]);

    const onFeedbackDialogClose = useCallback(() => {
        setAlertDialogOpen(false);
        setFeedbackAlert(null);
    }, []);

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
                    isColumnEditable={isColumnEditable}
                />
            </div>

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
