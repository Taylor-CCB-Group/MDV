/// <reference types="vite/client" />
import "./wdyr";
import "./all_css";
import "../react/HmrHack";
import ChartManager from "../charts/ChartManager.js";
import BaseChart from "../charts/BaseChart";
import { changeURLParam } from "./desktop_index";
import { createElement, useEffect, useMemo, useState } from "react";
import { createRoot, type Root } from "react-dom/client";
import { createMdvPortal } from "@/react/react_utils";
import ProjectStateHandlerWrapper from "@/react/ProjectStateHandler";
import { getProjectName } from "./ProjectContext";
import { getProjectBootstrapContext, loadProjectRuntime, type LoadedProjectRuntime } from "./projectRuntime";
import type { DataSource } from "@/charts/charts";

declare global {
    interface Window {
        mdv: {
            ChartManager: typeof ChartManager;
            chartManager: ChartManager;
            chartTypes?: any;
            debugChart?: any;
        };
    }
}

//@ts-expect-error maybe we'll be less hacky about this later
window.mdv = {
    ChartManager,
    chartTypes: BaseChart.types,
};

type BootstrapStatus =
    | { kind: "loading" }
    | { kind: "ready"; runtime: LoadedProjectRuntime }
    | { kind: "error"; error: unknown };

let stateHandlerContainer: HTMLElement | null = null;
let stateHandlerRoot: Root | null = null;

function LoadState({ status }: { status: BootstrapStatus }) {
    const detail = useMemo(() => {
        if (status.kind !== "error") return "";
        if (status.error instanceof Error) return status.error.message;
        return `${status.error}`;
    }, [status]);

    if (status.kind === "loading") {
        return (
            <div style={{ padding: "24px", fontFamily: "sans-serif", textAlign: "center" }}>
                <h2>Loading project...</h2>
                <p>Fetching datasources, state, and views.</p>
            </div>
        );
    }
    if (status.kind === "error") {
        return (
            <div style={{ padding: "24px", fontFamily: "sans-serif", textAlign: "center" }}>
                <h2>Error loading project</h2>
                <p>{detail}</p>
                <button type="button" onClick={() => window.location.reload()}>
                    Retry
                </button>
            </div>
        );
    }
    return null;
}

function BootstrapApp() {
    const context = useMemo(() => getProjectBootstrapContext(), []);
    const [status, setStatus] = useState<BootstrapStatus>({ kind: "loading" });
    const [initialised, setInitialised] = useState(false);

    useEffect(() => {
        if (context.isPopout) {
            document.title = "MDV popout";
            return;
        }
        void getProjectName(Number(context.projectId), context.apiRoot).then((projectName) => {
            document.title = `MDV - ${projectName}`;
        });
    }, [context]);

    useEffect(() => {
        if (context.isPopout) return;
        let cancelled = false;
        setStatus({ kind: "loading" });
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
            }
        };

        new ChartManager(
            "holder",
            runtime.datasources as DataSource[],
            runtime.dataLoader,
            runtime.config,
            listener as any,
        );
        setInitialised(true);
    }, [status, initialised]);

    if (context.isPopout || (status.kind === "ready" && initialised)) {
        return null;
    }
    return <LoadState status={status} />;
}

const holder = document.getElementById("holder");
if (!holder) {
    throw new Error("Could not find #holder bootstrap container.");
}

createRoot(holder).render(<BootstrapApp />);
