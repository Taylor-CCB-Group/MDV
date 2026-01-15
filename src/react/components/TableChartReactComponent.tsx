import { observer } from "mobx-react-lite";
import { useCallback, useState } from "react";
import { SlickgridReact } from "slickgrid-react";
import FindAndReplaceDialog from "./FindAndReplaceDialog";
import useFindReplace from "../hooks/useFindReplace";
import useSlickGridReact from "../hooks/useSlickGridReact";
import useEditCell from "../hooks/useEditCell";
import DebugErrorComponent from "@/charts/dialogs/DebugErrorComponent";
import AlertErrorComponent, { type AlertType } from "@/charts/dialogs/AlertErrorComponent";
import ReusableDialog from "@/charts/dialogs/ReusableDialog";

export type FeedbackAlert = {
    type: AlertType;
    message: string;
    stack?: string;
    traceback?: string;
    title?: string;
    metadata?: object;
} | null;

export type FeedbackAlertComponentType = {
    feedbackAlert: FeedbackAlert;
};

const isDebugError = (feedbackAlert: FeedbackAlert) => {
    if (feedbackAlert)
        return (
            feedbackAlert.type === "error" && (feedbackAlert.stack || feedbackAlert.traceback || feedbackAlert.metadata)
        );
    else return false;
};

const FeedbackAlertComponent = ({ feedbackAlert }: FeedbackAlertComponentType) => {
    if (!feedbackAlert) return <></>;

    if (isDebugError(feedbackAlert)) {
        return (
            <DebugErrorComponent
                error={{
                    message: feedbackAlert.message,
                    stack: feedbackAlert?.stack,
                    traceback: feedbackAlert?.traceback,
                }}
                title={feedbackAlert.title || "Error"}
                extraMetadata={feedbackAlert.metadata}
            />
        );
    } else {
        return (
            <AlertErrorComponent
                message={feedbackAlert.message}
                title={feedbackAlert.title || feedbackAlert.type}
                alertType={feedbackAlert.type}
            />
        );
    }
};

// todo: add integration tests using playwright to test the working of this component
// playwright would be more suitable to test this component rather than vitest which would require a lot of mocks
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
        isSelectingRef,
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
        isSelectingRef,
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
                    <ReusableDialog
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
