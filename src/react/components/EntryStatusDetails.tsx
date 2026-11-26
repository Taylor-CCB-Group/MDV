import { Button } from "@mui/material";
import type { EntryStatus } from "@spatialdata/vis";

import { elementIssueId } from "@/react/spatialdata/data_quality";
import { describeElementIssue, remediesFor } from "@/react/spatialdata/entry_status";

/**
 * One layer's status as sd.js reports it: why it could not load, with a retry where sd.js
 * says one could help, and what is degraded but still drawing. An element stored in a way
 * that makes it slow lists each issue and, once, how to fix it.
 */
export default function EntryStatusDetails({
    status,
    onRetry,
}: {
    status: EntryStatus;
    /** Absent while the canvas is not up to load the layer again. */
    onRetry?: () => void;
}) {
    return (
        <div className="grid gap-1">
            {status.errors.map(({ resource, error }) => (
                <div key={`error:${resource}:${error.kind}`} role="alert" className="flex items-start gap-2">
                    <span className="grow wrap-break-word">{error.message}</span>
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
                            {notice.issues.map((issue) => (
                                <li key={elementIssueId(issue)}>{describeElementIssue(issue)}</li>
                            ))}
                        </ul>
                        {remediesFor(notice.issues).map((remedy) => (
                            <p key={remedy} className="opacity-80">
                                {remedy}
                            </p>
                        ))}
                    </div>
                ) : (
                    <p key={notice.kind} role="status" className="wrap-break-word">
                        {notice.message}
                    </p>
                ),
            )}
        </div>
    );
}
