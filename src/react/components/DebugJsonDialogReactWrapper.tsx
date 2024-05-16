import { observer } from "mobx-react-lite";
import { BaseDialog } from "../../utilities/Dialog";
import { createEl } from "../../utilities/ElementsTyped";
import { createMdvPortal } from "@/react/react_utils";
import Gui from "./DebugJsonDialogComponent";
import BaseChart from "../../charts/BaseChart";

const DebugChart = observer(({chart, header}: {chart: any, header?: string}) => {
    return (<Gui json={chart} header={header} />)
});


// don't necessarily want to inherit from BaseDialog, could consider different approach.
// this will be more consistent / less work in short-term, and a basis for refactoring later.
class DebugChartReactWrapper extends BaseDialog {
    _root: ReturnType<typeof createMdvPortal>;
    get root() { return this._root; }
    set root(v) {
        this._root = v;
    }
    constructor(json: any, chart?: BaseChart) {
        const name = chart ? (chart.config.title || chart.config.type + ' ' + chart.config.id) : '';
        const config = { //TODO review popout behavior, use `__doc` or whatever here instead of `document` when appropriate
            width: 500, title: `Debug ${name}`, doc: document,
            onclose: () => { chart?.dialogs.splice(chart.dialogs.indexOf(this), 1) }
        };
        super(config, json);
    }
    init(parent: any) {
        const div = createEl('div', {}, this.dialog);
        this.root = createMdvPortal((
            <DebugChart chart={parent} />
        ), div);
    }
    close() {
        super.close();
        // this.root is undefined; wtf?
        // tried adding getter/setter for it - it doesn't got re-assigned.
        // I don't think there's any 'this' weirdness going on here...
        // some kind of react voodoo?
        if (this.root) this.root.unmount();
        else console.warn('not unmounting react root for dialog because of weirdness');
    }
}

export default DebugChartReactWrapper;