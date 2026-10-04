import { useEffect, useMemo, useState } from "react";
import { VivChannelList } from "@/react/components/ColorChannelComponents";
import {
    SpatialImagePanelContext,
    type SpatialImagePanelContextValue,
} from "@/react/components/spatialLayers/ImageLayerPanel";
import {
    VivProvider,
    createVivStores,
} from "@/react/components/avivatorish/state";
import type { SpatialThumbnailScatter } from "./SpatialThumbnailScatter";

const EMPTY_RASTER = { width: 0, height: 0, data: new Float32Array() };

/** Shared Viv channel list. One mix is applied to every thumbnail. */
export function ThumbnailChannelPanel({ chart }: { chart: SpatialThumbnailScatter }) {
    const stores = useMemo(() => createVivStores(), []);
    const [version, setVersion] = useState(0);
    useEffect(() => chart.subscribe(() => setVersion((n) => n + 1)), [chart]);

    const channels = chart.getChannels();
    const tone = chart.getTone();
    const names = chart.getChannelNames();

    useEffect(() => {
        const slots = chart.histogramSlots();
        stores.channelsStore.setState({
            raster: slots.map((slot) => slot.raster ?? EMPTY_RASTER),
            domains: slots.map((slot) => slot.domain),
            ids: channels.ids,
            selections: channels.selections,
        });
        stores.viewerStore.setState({
            isChannelLoading: channels.ids.map(() => false),
            isViewerLoading: false,
        });
    }, [channels, stores, version]);

    const value = useMemo<SpatialImagePanelContextValue>(() => {
        return {
            layerId: "spatial-thumbnail-scatter",
            loader: null,
            channelNames: names.length ? names : channels.ids.map((_, i) => `Channel ${i + 1}`),
            channelIds: channels.ids,
            colors: channels.colors,
            contrastLimits: channels.contrastLimits,
            channelsVisible: channels.channelsVisible,
            brightness: tone.brightness,
            contrast: tone.contrast,
            selections: channels.selections,
            setChannels: (patch) => {
                chart.setChannels(patch as Parameters<SpatialThumbnailScatter["setChannels"]>[0]);
            },
            addChannel: () => chart.addChannel(),
            addChannelWithSelection: () => chart.addChannel(),
            removeChannel: (index) => chart.removeChannel(index),
            patchVivLayerProps: (patch) => chart.patchTone(patch),
            patchToneAtIndex: (index, key, value) => chart.patchToneAtIndex(index, key, value),
        };
        // version refreshes arrays after chart edits
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [chart, version, names, channels, tone]);

    return (
        <VivProvider vivStores={stores}>
            <SpatialImagePanelContext.Provider value={value}>
                <div className="max-h-[50vh] overflow-auto border-b border-[hsl(var(--border))]">
                    <p className="px-2 pt-2 text-xs text-[hsl(var(--muted-foreground))]">
                        Channels apply to every image.
                    </p>
                    <VivChannelList />
                </div>
            </SpatialImagePanelContext.Provider>
        </VivProvider>
    );
}
