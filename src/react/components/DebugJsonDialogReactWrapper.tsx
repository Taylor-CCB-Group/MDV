import { ReactQueryDevtools } from "@tanstack/react-query-devtools";
import { observer } from "mobx-react-lite";
import type { ReactNode } from "react";
import { createMdvPortal } from "@/react/react_utils";
import type BaseChart from "../../charts/BaseChart";
import { BaseDialog } from "../../utilities/Dialog";
import { createEl } from "../../utilities/ElementsTyped";
import Gui from "./DebugJsonDialogComponent";

type DebugChartReactWrapperOptions = {
    showValidationSection?: boolean;
    /** Top-level keys of the json to show as groups of sections. */
    groupKeys?: readonly string[];
};

type DebugDialogContent = {
    json: any;
    showValidationSection?: boolean;
    chartPanel?: ReactNode;
    groupKeys?: readonly string[];
    jsonTitle?: string;
};

const DebugChart = observer(
    ({
        chart,
        header,
        showValidationSection = false,
        chartPanel,
        groupKeys,
        jsonTitle,
    }: {
        chart: any;
        header?: string;
        showValidationSection?: boolean;
        chartPanel?: ReactNode;
        groupKeys?: readonly string[];
        jsonTitle?: string;
    }) => {
        return (
            <>
                <Gui
                    json={chart}
                    header={header}
                    showValidationSection={showValidationSection}
                    chartPanel={chartPanel}
                    groupKeys={groupKeys}
                    jsonTitle={jsonTitle}
                />
                <ReactQueryDevtools initialIsOpen={false} />
            </>
        );
    },
);

// don't necessarily want to inherit from BaseDialog, could consider different approach.
// this will be more consistent / less work in short-term, and a basis for refactoring later.
class DebugChartReactWrapper extends BaseDialog {
    _root?: ReturnType<typeof createMdvPortal>;
    _resizeObserver?: ResizeObserver;
    get root() {
        return this._root;
    }
    set root(v) {
        this._root = v;
    }
    constructor(json: any, chart?: BaseChart<any>, options: DebugChartReactWrapperOptions = {}) {
        const name = chart ? chart.config.title || `${chart.config.type} ${chart.config.id}` : "";
        const doc = chart ? chart.__doc__ : document;
        const config = {
            //TODO review popout behavior, use `__doc` or whatever here instead of `document` when appropriate
            width: 500,
            title: chart ? `Debug / Report ${name}` : "Settings",
            doc,
            onclose: () => {
                chart?.dialogs.splice(chart.dialogs.indexOf(this), 1);
            },
        };
        super(config, {
            json,
            showValidationSection: options.showValidationSection,
            chartPanel: chart?.getDebugPanel?.(),
            groupKeys: options.groupKeys,
            // a chart's config reads best whole, not split by key
            jsonTitle: chart ? "Chart config" : undefined,
        });
    }
    init(content: DebugDialogContent) {
        const div = createEl("div", {}, this.dialog);
        this.root = createMdvPortal(
            <DebugChart
                chart={content.json}
                showValidationSection={content.showValidationSection}
                chartPanel={content.chartPanel}
                groupKeys={content.groupKeys}
                jsonTitle={content.jsonTitle}
            />,
            div,
            this,
        );
        this.repositionWhileContentSettles(div);
    }
    /**
     * BaseDialog positions itself once, in its constructor, before the React content has rendered; a tall dialog
     * then ends up hanging off the bottom of the window, so the page scrolls as well as the dialog.
     * Re-position as the content first takes shape (not afterwards, so we don't fight the user dragging it).
     */
    private repositionWhileContentSettles(content: HTMLElement) {
        if (typeof ResizeObserver === "undefined") return;
        let settleTimer: ReturnType<typeof setTimeout> | undefined;
        const observer = new ResizeObserver(() => {
            if (content.offsetHeight === 0) return;
            this._setDialogPosition();
            clearTimeout(settleTimer);
            settleTimer = setTimeout(() => observer.disconnect(), 500);
        });
        observer.observe(content);
        this._resizeObserver = observer;
    }
    close() {
        this._resizeObserver?.disconnect();
        super.close();
        // this.root is undefined; wtf?
        // tried adding getter/setter for it - it doesn't got re-assigned.
        // I don't think there's any 'this' weirdness going on here...
        // some kind of react voodoo?
        if (this.root) this.root.unmount();
        else console.warn("not unmounting react root for dialog because of weirdness");
    }
}

export default DebugChartReactWrapper;
