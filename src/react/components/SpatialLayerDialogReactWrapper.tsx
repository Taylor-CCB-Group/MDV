import { observer } from "mobx-react-lite";
import { BaseDialog } from "../../utilities/Dialog";
import { createEl } from "../../utilities/ElementsTyped";
import { createMdvPortal } from "@/react/react_utils";
import { ChartProvider } from "../context";
import { SpatialCanvasRendererProvider } from "@/react/spatial_canvas_renderer_context";
import SpatialLayerDialogComponent from "./SpatialLayerDialogComponent";
import type { SpatialDataMdvReact } from "./SpatialDataMDVReact";

const SpatialLayerDialogReact = observer(function SpatialLayerDialogReact({
    parent,
}: {
    parent: SpatialDataMdvReact;
}) {
    const rendererContext = parent.spatialRendererContext;
    if (!rendererContext) {
        return <SpatialLayerDialogComponent />;
    }
    return (
        <SpatialCanvasRendererProvider value={rendererContext}>
            <SpatialLayerDialogComponent />
        </SpatialCanvasRendererProvider>
    );
});

class SpatialLayerDialogReactWrapper extends BaseDialog {
    _root?: ReturnType<typeof createMdvPortal>;
    get root() {
        return this._root;
    }
    set root(v) {
        this._root = v;
    }

    constructor(parent: SpatialDataMdvReact) {
        if (parent.__doc__ !== document) {
            console.warn(
                "SpatialLayerDialogReactWrapper may have styling issues in popouts...",
            );
        }
        const config = {
            width: 560,
            maxHeight: 720,
            title: `Layers (${parent.config.title})`,
            doc: parent.__doc__ || document,
            onclose: () => {
                parent.layerDialog = undefined;
                parent.dialogs.splice(parent.dialogs.indexOf(this), 1);
            },
        };
        super(config, parent);
    }

    init(parent: SpatialDataMdvReact) {
        const div = createEl(
            "div",
            {
                styles: {
                    display: "block",
                    width: "100%",
                    height: "100%",
                },
            },
            this.dialog,
        );
        this.root = createMdvPortal(
            <ChartProvider chart={parent}>
                <SpatialLayerDialogReact parent={parent} />
            </ChartProvider>,
            div,
            this,
        );
    }

    close() {
        super.close();
        if (this.root) this.root.unmount();
        else console.warn("not unmounting react root for spatial layer dialog");
    }
}

BaseDialog.experiment["SpatialLayerDialogReact"] = SpatialLayerDialogReactWrapper;
export default SpatialLayerDialogReactWrapper;
