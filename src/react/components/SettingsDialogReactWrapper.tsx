import { observer } from "mobx-react-lite";
import { BaseDialog } from "../../utilities/Dialog";
import { createEl } from "../../utilities/ElementsTyped";
import { createMdvPortal } from "@/react/react_utils";
import Gui from "./SettingsDialogComponent";
import type BaseChart from "@/charts/BaseChart";

const SettingsDialog = observer(<T,>({ chart }: { chart: BaseChart<T> }) => {
    // const config = chart.getConfig(); //instrument with mobx etc
    return <Gui chart={chart} />;
});

// don't necessarily want to inherit from BaseDialog, could consider different approach.
// this will be more consistent / less work in short-term, and a basis for refactoring later.
class SettingsDialogReactWrapper<T> extends BaseDialog {
    _root?: ReturnType<typeof createMdvPortal>;
    get root() {
        return this._root;
    }
    set root(v) {
        this._root = v;
    }
    constructor(chart: BaseChart<T>, position?: [number, number]) {
        // if this is intended to be a drop-in replacement for existing SettingsDialog,
        // it isn't only used by 'charts', but e.g. tracks.
        const name =
            chart.config.title || `${chart.config.type} ${chart.config.id}`;
        const config = {
            width: 500,
            title: `Settings (${name})`,
            doc: chart.__doc__ || document,
            position,
            onclose: () => {
                chart.dialogs.splice(chart.dialogs.indexOf(this), 1);
                if (chart.settingsDialog === this) chart.settingsDialog = undefined;
            },
        };
        super(config, chart);
    }
    init(chart: BaseChart<T>) {
        const div = createEl("div", {}, this.dialog);
        this.root = createMdvPortal(
            <SettingsDialog chart={chart} />,
            div,
            this,
        );
    }
    close() {
        super.close();
        // this.root is undefined; wtf?
        // tried adding getter/setter for it - it doesn't got re-assigned.
        // I don't think there's any 'this' weirdness going on here...
        // some kind of react voodoo?
        if (this.root) this.root.unmount();
        else
            console.warn(
                "not unmounting react root for dialog because of weirdness",
            );
    }
}

export default SettingsDialogReactWrapper;
