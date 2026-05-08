/// <reference types="vite/client" />
import "./wdyr";
import "./all_css";
import { createElement, useEffect, useMemo, useRef, useState } from "react";
import { createRoot, type Root } from "react-dom/client";
import { getProjectName } from "./ProjectContext";
import { getProjectBootstrapContext, loadProjectRuntime, type LoadedProjectRuntime } from "./projectRuntime";
import type { DataSource } from "@/charts/charts";
import type ChartManager from "@/charts/ChartManager";

type CreateMdvPortalFn = typeof import("@/react/react_utils").createMdvPortal;
type ProjectStateHandlerWrapperType = typeof import("@/react/ProjectStateHandler").default;

declare global {
    interface Window {
        mdv: {
            ChartManager: typeof ChartManager | undefined;
            chartManager: ChartManager;
            chartTypes?: any;
            debugChart?: any;
        };
    }
}

type BootstrapDeps = {
    ChartManager: typeof ChartManager;
    chartTypes: any;
    changeURLParam: (key: string, value: string) => void;
    createMdvPortal: CreateMdvPortalFn;
    ProjectStateHandlerWrapper: ProjectStateHandlerWrapperType;
};
let bootstrapDepsPromise: Promise<BootstrapDeps> | null = null;

function loadBootstrapDeps(): Promise<BootstrapDeps> {
    if (!bootstrapDepsPromise) {
        bootstrapDepsPromise = Promise.all([
            import("../charts/ChartManager.js"),
            import("../charts/BaseChart"),
            import("../react/HmrHack"),
            import("./desktop_index"),
            import("@/react/react_utils"),
            import("@/react/ProjectStateHandler"),
        ]).then(
            ([chartManagerModule, baseChartModule, _hmrHackModule, desktopIndexModule, reactUtilsModule, projectStateModule]) => ({
                ChartManager: chartManagerModule.default,
                chartTypes: baseChartModule.default.types,
                changeURLParam: desktopIndexModule.changeURLParam,
                createMdvPortal: reactUtilsModule.createMdvPortal,
                ProjectStateHandlerWrapper: projectStateModule.default,
            }),
        );
    }
    return bootstrapDepsPromise;
}

type BootstrapStatus =
    | { kind: "loading" }
    | { kind: "ready"; runtime: LoadedProjectRuntime }
    | { kind: "error"; error: unknown };

let stateHandlerContainer: HTMLElement | null = null;
let stateHandlerRoot: Root | null = null;

function LoadState({
    status,
    detail,
    isClosing,
    onTransitionEnd,
}: {
    status: BootstrapStatus;
    detail: string;
    isClosing: boolean;
    onTransitionEnd: () => void;
}) {
    const appTheme = window.mdv?.chartManager?.theme;
    const prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
    const isDark = appTheme ? appTheme === "dark" : prefersDark;
    const progressText =
        status.kind === "loading"
            ? "Fetching datasources, state, and views."
            : "Preparing chart runtime and rendering initial view.";
    const overlayClassName = [
        "fixed inset-0 z-[4000] flex items-center justify-center p-6 text-center transition-opacity duration-200",
        isDark ? "bg-slate-900/95" : "bg-slate-50/95",
        isClosing ? "opacity-0 pointer-events-none" : "opacity-100",
    ].join(" ");
    const panelClassName = [
        "w-full max-w-[520px] rounded-xl border px-6 py-5 shadow-2xl",
        isDark ? "border-slate-400/25 bg-slate-800/85 text-slate-200" : "border-slate-900/10 bg-white/90 text-slate-900",
    ].join(" ");
    const bodyTextClassName = isDark ? "mb-3 text-slate-400" : "mb-3 text-slate-600";
    const trackClassName = [
        "relative h-1.5 overflow-hidden rounded-full",
        isDark ? "bg-slate-700/80" : "bg-slate-300/80",
    ].join(" ");
    const indicatorClassName = [
        "absolute left-0 top-0 h-full w-2/5 rounded-full animate-pulse",
        isDark ? "bg-blue-400" : "bg-blue-600",
    ].join(" ");
    const spinnerClassName = [
        "mx-auto mb-3 h-6 w-6 animate-spin rounded-full border-2 border-slate-400/40",
        isDark ? "border-t-blue-400" : "border-t-blue-600",
    ].join(" ");
    const titleClassName = "mb-2 text-[20px] font-semibold leading-6 tracking-[-0.01em]";
    const progressClassName = `${bodyTextClassName} min-h-[24px] text-[14px] leading-6`;

    return (
        <div
            onTransitionEnd={onTransitionEnd}
            className={overlayClassName}
            style={{ fontFamily: "Roboto, 'Helvetica Neue', Arial, sans-serif" }}
        >
            <div className={panelClassName}>
                {status.kind === "error" ? (
                    <>
                        <h2 className={titleClassName}>Error loading project</h2>
                        <p className={progressClassName}>{detail}</p>
                        <button
                            type="button"
                            onClick={() => window.location.reload()}
                            className="rounded-md border border-slate-400/50 px-3 py-1.5 text-sm font-medium hover:bg-slate-500/15"
                        >
                            Retry
                        </button>
                    </>
                ) : (
                    <>
                        <div className={spinnerClassName} />
                        <h2 className={titleClassName}>Loading project...</h2>
                        <p className={progressClassName}>{progressText}</p>
                        <div className={trackClassName}>
                            <div className={indicatorClassName} />
                        </div>
                    </>
                )}
            </div>
        </div>
    );
}

function BootstrapApp({ onComplete }: { onComplete: () => void }) {
    const context = useMemo(() => getProjectBootstrapContext(), []);
    const [status, setStatus] = useState<BootstrapStatus>({ kind: "loading" });
    const [initialised, setInitialised] = useState(false);
    const [initialViewLoaded, setInitialViewLoaded] = useState(false);
    const [isClosing, setIsClosing] = useState(false);
    const completeCalledRef = useRef(false);
    const detail = useMemo(() => {
        if (status.kind !== "error") return "";
        if (status.error instanceof Error) return status.error.message;
        return `${status.error}`;
    }, [status]);

    useEffect(() => {
        if (context.isPopout) {
            document.title = "MDV popout";
            onComplete();
            return;
        }
        void getProjectName(Number(context.projectId), context.apiRoot).then((projectName) => {
            document.title = `MDV - ${projectName}`;
        });
    }, [context, onComplete]);

    useEffect(() => {
        if (context.isPopout) return;
        let cancelled = false;
        setStatus({ kind: "loading" });
        // Start heavy module loading in parallel with initial data fetch.
        void loadBootstrapDeps();
        void loadProjectRuntime(context)
            .then((runtime) => {
                if (!cancelled) setStatus({ kind: "ready", runtime });
            })
            .catch((error) => {
                if (!cancelled) setStatus({ kind: "error", error });
            });
        return () => {
            cancelled = true;
        };
    }, [context]);

    useEffect(() => {
        if (status.kind !== "ready" || initialised) return;
        const { runtime } = status;
        let cancelled = false;

        void loadBootstrapDeps()
            .then(
                ({
                    ChartManager,
                    chartTypes,
                    changeURLParam,
                    createMdvPortal,
                    ProjectStateHandlerWrapper,
                }) => {
                    window.mdv = {
                        ...(window.mdv ?? {}),
                        ChartManager,
                        chartTypes,
                    };

                    const listener = async (type: string, cm: ChartManager, data: any) => {
                        if (type === "state_saved") {
                            if (stateHandlerRoot) {
                                stateHandlerRoot.unmount();
                            }

                            if (!stateHandlerContainer) {
                                stateHandlerContainer = document.createElement("div");
                                stateHandlerContainer.id = "mdv_state_handler";
                                document.body.appendChild(stateHandlerContainer);
                            }
                            stateHandlerRoot = createMdvPortal(
                                createElement(ProjectStateHandlerWrapper, {
                                    root: runtime.root,
                                    data,
                                    staticFolder: runtime.staticFolder,
                                    permission: runtime.permission,
                                }),
                                stateHandlerContainer,
                                undefined,
                                cm,
                            );
                        }
                        if (type === "view_loaded") {
                            changeURLParam("view", cm.viewManager.current_view);
                            setInitialViewLoaded(true);
                        }
                    };

                    if (cancelled) return;
                    new ChartManager(
                        "holder",
                        runtime.datasources as DataSource[],
                        runtime.dataLoader,
                        runtime.config,
                        listener as any,
                    );
                    setInitialised(true);
                },
            )
            .catch((error) => {
                if (!cancelled) {
                    setStatus({ kind: "error", error });
                }
            });

        return () => {
            cancelled = true;
        };
    }, [status, initialised]);

    useEffect(() => {
        if (status.kind !== "error" && initialised && initialViewLoaded) {
            setIsClosing(true);
        }
    }, [status, initialised, initialViewLoaded]);

    useEffect(() => {
        if (!isClosing) return;
        const timeout = window.setTimeout(() => {
            if (!completeCalledRef.current) {
                completeCalledRef.current = true;
                onComplete();
            }
        }, 320);
        return () => {
            window.clearTimeout(timeout);
        };
    }, [isClosing, onComplete]);

    const handleTransitionEnd = () => {
        if (isClosing && !completeCalledRef.current) {
            completeCalledRef.current = true;
            onComplete();
        }
    };

    return (
        <LoadState
            status={status}
            detail={detail}
            isClosing={isClosing}
            onTransitionEnd={handleTransitionEnd}
        />
    );
}

const holder = document.getElementById("holder");
if (!holder) {
    throw new Error("Could not find #holder bootstrap container.");
}

const bootstrapOverlayContainer = document.createElement("div");
bootstrapOverlayContainer.id = "mdv_bootstrap_overlay";
document.body.appendChild(bootstrapOverlayContainer);

const bootstrapRoot = createRoot(bootstrapOverlayContainer);
const closeBootstrapOverlay = () => {
    bootstrapRoot.unmount();
    bootstrapOverlayContainer.remove();
};

bootstrapRoot.render(<BootstrapApp onComplete={closeBootstrapOverlay} />);
