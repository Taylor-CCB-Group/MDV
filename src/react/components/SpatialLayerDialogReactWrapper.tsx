import { useEffect } from "react";
import { observer } from "mobx-react-lite";
import { SpatialDataProvider } from "@spatialdata/react";
import { BaseDialog } from "../../utilities/Dialog";
import { createEl } from "../../utilities/ElementsTyped";
import { createMdvPortal } from "@/react/react_utils";
import { getProjectURL } from "@/dataloaders/DataLoaderUtil";
import { ensureChunkWorker } from "@/react/spatialdata/ensureChunkWorker";
import { ChartProvider } from "../context";
import { useRegion } from "../hooks";
import SpatialLayerDialogComponent from "./SpatialLayerDialogComponent";
import type { SpatialDataMdvReact } from "./SpatialDataMDVReact";

function getSpatialDataUrl(region: unknown): string | undefined {
    if (!region || typeof region !== "object" || !("spatial" in region)) return undefined;
    const spatial = region.spatial;
    if (!spatial || typeof spatial !== "object" || !("file" in spatial)) return undefined;
    const file = spatial.file;
    if (typeof file !== "string" || !file) return undefined;
    return getProjectURL(`spatial/${file}`);
}

const SpatialLayerDialogReact = observer(function SpatialLayerDialogReact() {
    const rawRegion = useRegion();
    const spatialDataUrl = getSpatialDataUrl(rawRegion);
    useEffect(() => {
        ensureChunkWorker();
    }, []);
    return (
        <SpatialDataProvider source={spatialDataUrl}>
            <SpatialLayerDialogComponent />
        </SpatialDataProvider>
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
            height: 700,
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
                <SpatialLayerDialogReact />
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
