import BaseChart from "../../charts/BaseChart";
import DeckGL from "deck.gl/typed";
import { ColorPaletteExtension, DetailView, getDefaultInitialViewState, ScaleBarLayer, VivViewer } from "@hms-dbmi/viv";
import { useChannelStats, useChartID, useChartSize, useConfig, useImgUrl } from "../hooks";
import { useScatterplotLayer } from "../scatter_state";
import { BaseReactChart } from "./BaseReactChart";
import { observer } from "mobx-react-lite";
import { action, makeObservable, observable } from "mobx";
import { BaseDialog } from "../../utilities/Dialog";
import { ChannelsState, DEFAUlT_CHANNEL_STATE, ROI, VivConfig, VivProvider, useChannelsState, useViewerStoreApi, useVivLayerConfig } from "./avivatorish/state";
import "../../charts/VivScatterPlot"; //because we use the BaseChart.types object, make sure it's loaded.
import { OmeTiffProvider, useOmeTiff } from "../context"; 
import { useEffect, useMemo, useState } from "react";
import MDVivViewer, { getVivId } from "./avivatorish/MDVivViewer";
import SelectionOverlay from "./SelectionOverlay";

function ReactTest() {
    // to make this look more like Avivator...
    // we probably don't want OmeTiffProvider to be a thing...
    // we should use a VivProvider, with hooks that look more like Avivator's
    // so VivProvider should have whatever is necessary to adapt our config to that
    // and we'd useLoader() as opposed to useOmeTiff()
    // ... and hopefully our version of Avivator hooks will have better types ...
    return (
    <OmeTiffProvider>
        <VivProvider>
            <MainChart />
        </VivProvider>
    </OmeTiffProvider>
    )
}

const DeckImpl = observer(() => {
    const imgUrl = useImgUrl();
    const viewerStore = useViewerStoreApi();
    useEffect(() => {
        if (!imgUrl) return;
        viewerStore.setState({ source: imgUrl });
    }, [imgUrl, viewerStore]);
    const ome = useOmeTiff();
    const config = useConfig<VivMdvReactConfig>();
    const [width, height] = useChartSize();
    const id = useChartID();
    // const detailId = getVivId(id + 'detail-react');
    const detailId = id + 'detail-react';

    const channelX = config.channel; //don't do this... use the proper image_properties
    // and make it populate the channel state if necessary.
    const stats = useChannelStats(ome, channelX);
    const channelsState = useChannelsState();
    const userSet = channelsState.contrastLimits && channelsState.contrastLimits.length > 0; //hack, for now.
    const contrastLimits = userSet ? channelsState.contrastLimits : [stats ? stats.contrastLimits : [0, 1]];

    const [scatterplotLayer, hoverInfo] = useScatterplotLayer();

    const layerConfig = useVivLayerConfig();
    const [viewState, setViewState] = useState<ReturnType<typeof getDefaultInitialViewState>>();
    useEffect(() => {
        if (!ome) return;
        if (!viewState) {
            //WIP
            setViewState(getDefaultInitialViewState(ome.data, { width, height }));
        }
    }, [ome]);
    const extensions = useMemo(() => [new ColorPaletteExtension()], []);
    const detailView = useMemo(() => new DetailView({
        id: detailId,
        snapScaleBar: true,
        width, height
    }), [id, width, height]);
    // TODO get viv working in popouts (not a react thing - happens elsewhere
    // - probably need to handle lost gl context)
    const layerConfigX = {
        loader: ome.data,
        selections: [{ z: 0, c: channelX, t: 0 }],
        contrastLimits,
        extensions,
        colors: [[80, 255, 255]],
        channelsVisible: [true],
    }
    // pending proper channel state handling... show that we can set contrast limits.
    if (userSet) layerConfigX.contrastLimits = contrastLimits;
    const deckProps = {
        // initialViewState: viewState,
        // controller: true,
        getTooltip: ({ object }) => hoverInfo && hoverInfo.index !== -1 && 'i: '+hoverInfo.index,
        style: {
            zIndex: '-1',
        },
        // for now I disabled layerFilter in VivViewer, because this wasn't passing
        // but it should work if I make sure it has an id that matches expectations.
        layers: [scatterplotLayer],
        id: id + 'deck',
    }
    return (
        <>
            <SelectionOverlay scatterplotLayer={scatterplotLayer}/>
            <VivViewer 
            views={[detailView]}
            layerProps={[layerConfigX]}
            viewStates={[{...viewState, id: detailId}]}
            onViewStateChange={e => {
                // Need to clarify what is idiomatic to Viv and what is appropriate for us.
                // viewerStore is not the place for keeping track of transforms, not sure if Avivator has an equivalent.
                //xxx --- viewState in avivator is an object [id]: viewState, so we may need to do something similar.
                viewerStore.setState({ viewState: {...e.viewState, id: detailId } });
            }}
            
            deckProps={deckProps} 
            />
        </>
    );
});

const MainChart = () => {
    const ome = useOmeTiff();
    if (!ome) return <div>Loading...</div>; // todo suspense.
    return <DeckImpl />;
};

/// some type stuff... common bits should probably move elsewhere...
// even if the type isn't constrained beyond 'string' (which potentially in some cases it could be)
// this kind of thing can potentially help DX and encourage appropriate values to be passed.
type ColumnName = string; 
//could we infer or something to avoid having to repeat this?
export type ScatterPlotConfig = {
    radius: number,
    opacity: number,
    color_by: ColumnName,
    color_legend: {
        display: boolean,
        // todo: add more options here...
    },
};
const scatterDefaults: ScatterPlotConfig = {
    radius: 1,
    opacity: 1,
    color_by: null,
    color_legend: {
        display: false,
    }
};
export type VivRoiConfig = {
    // making this 'type' very specific will let us infer the rest of the type, i.e.
    // `if (config.type === 'VivMdvRegionReact')` will narrow the type to VivScatterConfig
    // ... except that we also need to check the other condition, because 'string' could also be that.
    // so it's not completely ideal.
    type: "VivMdvRegionReact" | "viv_scatter_plot",
    background_filter: {
        column: ColumnName,
        category: string,
    },
    roi: ROI,
    viv: VivConfig,
    //image_properties: ChannelsState,
} & ScatterPlotConfig;

export type VivMdvReactConfig = ScatterPlotConfig & (
    { type: 'VivMdvReact', imageURL: string, overviewOn: boolean, image_properties: ChannelsState } 
    | VivRoiConfig
) & { channel: number };
export type VivMDVReact = VivMdvReact;
class VivMdvReact extends BaseReactChart<VivMdvReactConfig> {
    colorDialog: any;
    constructor(dataStore, div, config) {
        // todo better default config
        config = {...scatterDefaults, ...config};
        if (!config.channel) config.channel = 0;
        if (config.type === 'VivMdvRegionReact') {
            if (!config.viv.image_properties) config.viv.image_properties = DEFAUlT_CHANNEL_STATE;
        } else if (config.type === 'VivMdvReact') {
            if (config.overviewOn === undefined) config.overviewOn = false;
            if (config.image_properties === undefined) config.image_properties = DEFAUlT_CHANNEL_STATE;
        }
        super(dataStore, div, config, ReactTest);
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
        });
        this.addMenuIcon("fas fa-palette", "Alter Channels").addEventListener("click", (e) => {
            // return;
            if (!this.colorDialog) {
                //this.colorDialog = new ColorChannelDialog(this);
                // 🙄 HMR hack
                this.colorDialog = new BaseDialog.experiment['ColorDialogReact'](this);
            }
        });
    }
    remove() {
        super.remove();
        if (this.colorDialog) this.colorDialog.close();
    }
    colorBy?: (i: number) => [r: number, g: number, b: number];
    colorByColumn(col?: ColumnName) {
        if (!col) return this.colorByDefault();
        this.config.color_by = col;
        this.colorBy = this.getColorFunction(col, true);
    }
    colorByDefault() {
        this.config.color_by = null;
        this.colorBy = null;
    }
    getColorOptions() {
        return {
            colorby: "all"
        }
    }

    getSettings(): { type: string; label: string; current_value: any; func: (v: any) => void; }[] {
        const c = this.config;
        const settings = super.getSettings();
        return settings.concat([
            {
                // very crude placeholder for testing mechanism...
                type: "slider",
                label: "channel test",
                current_value: c.channel || 0,
                min: 0,
                max: 10,
                step: 1,
                continuous: true,
                func: x => {
                    c.channel = x;
                }
            },
            {
                type: "slider",
                label: "radius",
                current_value: c.radius || 1,
                min: 0,
                max: 50,
                continuous: true,
                func: x => {
                    c.radius = x;
                }
            },
            {
                type: "slider",
                label: "opacity",
                current_value: c.opacity || 1,
                min: 0,
                max: 1,
                continuous: true,
                func: x => {
                    c.opacity = x;
                }
            },
            // no longer using PictureInPictureViewer - would need to re-implement to some extent
            // in order for combined viewer to work (also add overviewOn to every config)
            // {
            //     type: "check",
            //     label: "overview",
            //     current_value: c.overviewOn || false,
            //     func: x => {
            //         c.overviewOn = x;
            //     }
            // }
        ]);
    }
}

BaseChart.types["VivMdvRegionReact"] = {
    ...BaseChart.types["viv_scatter_plot"],
    "class": VivMdvReact,
    name: "Viv Scatter Plot (react)",
}

export type VivMdvReactType = typeof VivMdvReact;
// we rely on the side-effect of this import to register the chart type
///- will register otherwise, but this tricks HMR into working...
export default 42;