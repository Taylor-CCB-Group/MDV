import AlertErrorComponent, { type AlertType } from "@/charts/dialogs/AlertErrorComponent";
import DebugErrorComponent from "@/charts/dialogs/DebugErrorComponent";

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

export const isDebugError = (feedbackAlert: FeedbackAlert) => {
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

export default FeedbackAlertComponent;