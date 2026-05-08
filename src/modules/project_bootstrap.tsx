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

if (import.meta.env.DEV) {
    void import("../react/HmrHack");
}

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
            import("./desktop_index"),
            import("@/react/react_utils"),
            import("@/react/ProjectStateHandler"),
        ]).then(
            ([chartManagerModule, baseChartModule, desktopIndexModule, reactUtilsModule, projectStateModule]) => ({
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
    const colors = isDark
        ? {
              bg: "rgba(15, 23, 42, 0.94)",
              panel: "rgba(30, 41, 59, 0.86)",
              text: "#e2e8f0",
              subtext: "#94a3b8",
              accent: "#60a5fa",
          }
        : {
              bg: "rgba(248, 250, 252, 0.96)",
              panel: "rgba(255, 255, 255, 0.9)",
              text: "#0f172a",
              subtext: "#475569",
              accent: "#2563eb",
          };
    const progressText =
        status.kind === "loading"
            ? "Fetching datasources, state, and views."
            : "Preparing chart runtime and rendering initial view.";

    return (
        <div
            onTransitionEnd={onTransitionEnd}
            style={{
                position: "fixed",
                inset: 0,
                display: "flex",
                justifyContent: "center",
                alignItems: "center",
                padding: "24px",
                fontFamily: "sans-serif",
                textAlign: "center",
                background: colors.bg,
                transition: "opacity 220ms ease",
                opacity: isClosing ? 0 : 1,
                pointerEvents: isClosing ? "none" : "auto",
                zIndex: 4000,
            }}
        >
            <style>{`
                @keyframes mdvBootstrapSpin {
                    from { transform: rotate(0deg); }
                    to { transform: rotate(360deg); }
                }
                @keyframes mdvBootstrapShimmer {
                    0% { transform: translateX(-120%); }
                    100% { transform: translateX(320%); }
                }
            `}</style>
            <div
                style={{
                    minWidth: "320px",
                    maxWidth: "520px",
                    padding: "22px 24px",
                    borderRadius: "12px",
                    background: colors.panel,
                    border: isDark ? "1px solid rgba(148,163,184,0.24)" : "1px solid rgba(15,23,42,0.1)",
                    boxShadow: isDark
                        ? "0 10px 30px rgba(2,6,23,0.35)"
                        : "0 10px 30px rgba(15,23,42,0.12)",
                }}
            >
                {status.kind === "error" ? (
                    <>
                        <h2 style={{ color: colors.text, margin: "0 0 10px" }}>Error loading project</h2>
                        <p style={{ color: colors.subtext, margin: "0 0 12px" }}>{detail}</p>
                        <button type="button" onClick={() => window.location.reload()}>
                            Retry
                        </button>
                    </>
                ) : (
                    <>
                        <div
                            style={{
                                width: "24px",
                                height: "24px",
                                borderRadius: "50%",
                                border: `2px solid ${isDark ? "rgba(148,163,184,0.35)" : "rgba(148,163,184,0.45)"}`,
                                borderTopColor: colors.accent,
                                animation: "mdvBootstrapSpin 0.9s linear infinite",
                                margin: "0 auto 12px",
                            }}
                        />
                        <h2 style={{ color: colors.text, margin: "0 0 8px" }}>Loading project...</h2>
                        <p style={{ color: colors.subtext, margin: "0 0 10px" }}>{progressText}</p>
                        <div
                            style={{
                                position: "relative",
                                height: "6px",
                                borderRadius: "999px",
                                overflow: "hidden",
                                background: isDark ? "rgba(51,65,85,0.7)" : "rgba(203,213,225,0.8)",
                            }}
                        >
                            <div
                                style={{
                                    position: "absolute",
                                    inset: 0,
                                    width: "40%",
                                    borderRadius: "999px",
                                    background: colors.accent,
                                    animation: "mdvBootstrapShimmer 1.2s ease-in-out infinite",
                                }}
                            />
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
