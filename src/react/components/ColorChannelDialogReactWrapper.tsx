import { observer } from "mobx-react-lite";
import { BaseDialog } from "../../utilities/Dialog";
import type { VivMDVReact } from "./VivMDVReact";
import { createEl } from "../../utilities/ElementsTyped";
import { createMdvPortal } from "@/react/react_utils";
import Gui from "./ColorChannelComponents";
import { ChartProvider } from "../context";
import type { VivContextType } from "./avivatorish/state";

const ColorChannelDialogReact = observer(function ColorChannelDialogReact({
    vivStores,
}: { vivStores: VivContextType }) {
    return <Gui vivStores={vivStores} />;
});

// don't necessarily want to inherit from BaseDialog, could consider different approach.
// this will be more consistent / less work in short-term, and a basis for refactoring later.
class ColorDialogReactWrapper extends BaseDialog {
    _root?: ReturnType<typeof createMdvPortal>;
    get root() {
        return this._root;
    }
    set root(v) {
        this._root = v;
    }
    constructor(parent: VivMDVReact) {
        if (parent.__doc__ !== document)
            console.warn(
                "ColorDialogReactWrapper may have styling issues in popouts...",
            );
        const config = {
            width: 500,
            title: `Color Channels (${parent.config.title})`,
            doc: parent.__doc__ || document,
            onclose: () => {
                parent.colorDialog = undefined;
                parent.dialogs.splice(parent.dialogs.indexOf(this), 1);
            },
        };
        super(config, parent);
    }
    init(parent: VivMDVReact) {
        const div = createEl("div", {}, this.dialog);
        this.root = createMdvPortal(
            <ChartProvider chart={parent}>
                <ColorChannelDialogReact vivStores={parent.vivStores} />
            </ChartProvider>,
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
        else {
            //! this seems to be more of a problem than I thought.
            console.warn("not unmounting react root for dialog because of weirdness");
        }
    }
}

BaseDialog.experiment["ColorDialogReact"] = ColorDialogReactWrapper;
export default ColorDialogReactWrapper;
