import type { ReactNode } from "react";

export type ChartArrayEmptyStateProps = {
    waitingForViewState?: boolean;
    waitingMessage?: string;
    configuredCellCount: number;
    loadedCellCount: number;
    noCellsConfiguredMessage: string;
    loadingCellsMessage: string;
};

export function ChartArrayEmptyState({
    waitingForViewState,
    waitingMessage = "Loading...",
    configuredCellCount,
    loadedCellCount,
    noCellsConfiguredMessage,
    loadingCellsMessage,
}: ChartArrayEmptyStateProps): ReactNode {
    if (waitingForViewState) {
        return (
            <div className="flex h-full items-center justify-center text-sm">{waitingMessage}</div>
        );
    }
    if (configuredCellCount === 0) {
        return (
            <div className="flex h-full items-center justify-center text-sm">
                {noCellsConfiguredMessage}
            </div>
        );
    }
    if (loadedCellCount === 0) {
        return (
            <div className="flex h-full items-center justify-center text-sm">{loadingCellsMessage}</div>
        );
    }
    return null;
}
