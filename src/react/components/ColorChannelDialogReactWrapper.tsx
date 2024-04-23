import { observer } from "mobx-react-lite";
import { BaseDialog } from "../../utilities/Dialog";
import type { VivMDVReact } from "./VivMDVReact";
import { createEl } from "../../utilities/ElementsTyped";
import { createRoot } from "react-dom/client";
import Gui from "./ColorChannelComponents";
import { ChartProvider } from "../context";

const ColorChannelDialogReact = observer(function ColorChannelDialogReact() {
    return (<Gui />)
});


// don't necessarily want to inherit from BaseDialog, could consider different approach.
// this will be more consistent / less work in short-term, and a basis for refactoring later.
class ColorDialogReactWrapper extends BaseDialog {
    _root: ReturnType<typeof createRoot>;
    get root() { return this._root; }
    set root(v) { 
        this._root = v;
    }
    constructor(parent: VivMDVReact) {
        const config = {
            width: 500, title: `Color Channels (${parent.config.title})`, doc: document,
            onclose: () => { (parent as any).colorDialog = null} //todo fix types
        };
        super(config, parent)
    }
    init(parent: VivMDVReact) {
        const div = createEl('div', {}, this.dialog);
        this.root = createRoot(div);
        // todo use a portal, share a root & chart context with the parent?
        this.root.render(
            <ChartProvider chart={parent}>
                <ColorChannelDialogReact />
            </ChartProvider>
        );
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

BaseDialog.experiment['ColorDialogReact'] = ColorDialogReactWrapper;
export default ColorDialogReactWrapper;