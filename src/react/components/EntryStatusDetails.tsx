import { Button } from "@mui/material";
import type { EntryStatus } from "@spatialdata/vis";

import { elementIssueId } from "@/react/spatialdata/data_quality";
import { describeElementIssue, remediesFor } from "@/react/spatialdata/entry_status";

function DismissButton({ label, onClick }: { label: string; onClick: () => void }) {
    return (
        <Button
            size="small"
            sx={{ flexShrink: 0, minWidth: 0, py: 0 }}
            onClick={onClick}
            aria-label={label}
            title={label}
        >
            Dismiss
        </Button>
    );
}

/**
 * One layer's status as sd.js reports it: why it could not load, with a retry where sd.js
 * says one could help, and what is degraded but still drawing. An element stored in a way
 * that makes it slow lists each issue and, once, how to fix it.
 *
 * With `onDismiss`, each notice (each issue, for an element's storage) can be dismissed;
 * it is called with the part's id (`noticeParts`). Errors cannot be dismissed here.
 */
export default function EntryStatusDetails({
    status,
    onRetry,
    onDismiss,
}: {
    status: EntryStatus;
    /** Absent while the canvas is not up to load the layer again. */
    onRetry?: () => void;
    onDismiss?: (id: string) => void;
}) {
    return (
        <div className="grid gap-1">
            {status.errors.map(({ resource, error }) => (
                <div key={`error:${resource}:${error.kind}`} role="alert" className="flex items-start gap-2">
                    <span className="grow break-words">{error.message}</span>
                    {error.retryable && onRetry && (
                        <Button size="small" variant="outlined" onClick={onRetry}>
                            Retry
                        </Button>
                    )}
                </div>
            ))}
            {status.notices.map((notice) =>
                notice.kind === "element-issues" ? (
                    <div key={notice.kind} role="status">
                        <ul className="ml-4 list-disc">
                            {notice.issues.map((issue) => {
                                const id = elementIssueId(issue);
                                return (
                                    <li key={id}>
                                        <span className="flex items-start gap-2">
                                            <span className="grow">{describeElementIssue(issue)}</span>
                                            {onDismiss && (
                                                <DismissButton
                                                    label={`Dismiss: ${describeElementIssue(issue)}`}
                                                    onClick={() => onDismiss(id)}
                                                />
                                            )}
                                        </span>
                                    </li>
                                );
                            })}
                        </ul>
                        {remediesFor(notice.issues).map((remedy) => (
                            <p key={remedy} className="opacity-80">
                                {remedy}
                            </p>
                        ))}
                    </div>
                ) : (
                    <div key={notice.kind} role="status" className="flex items-start gap-2">
                        <span className="grow break-words">{notice.message}</span>
                        {onDismiss && (
                            <DismissButton
                                label={`Dismiss: ${notice.message}`}
                                onClick={() => onDismiss(notice.kind)}
                            />
                        )}
                    </div>
                ),
            )}
        </div>
    );
}
