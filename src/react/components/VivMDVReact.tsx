import BaseChart from "../../charts/BaseChart";
import DeckGL, { ScatterplotLayer } from "deck.gl/typed";
import { ColorPaletteExtension, DetailView, getDefaultInitialViewState } from "@hms-dbmi/viv";
import { useChannelStats, useChartID, useChartSize, useConfig, useFilteredIndices, useOmeTiff, useParamColumns } from "../hooks";
import { BaseConfig, BaseReactChart } from "./BaseReactChart";
import { observer } from "mobx-react-lite";
import { action, makeObservable, observable } from "mobx";
import { BaseDialog } from "../../utilities/Dialog";
import { ChannelsState, DEFAUlT_CHANNEL_STATE, ROI, VivConfig } from "../viv_state";
import "../../charts/VivScatterPlot"; //because we use the BaseChart.types object, make sure it's loaded.
import { useChart } from "../context"; 


// we need to add observer here, not just in BaseReactChart, for HMR to work.
const ReactTest = observer(() => {
    const config = useConfig<VivMdvReactConfig>();
    const [width, height] = useChartSize();
    const ome = useOmeTiff();
    // const {id} = config; //there should always be a persistent, reliable id here.
    const id = useChartID(); // hook in case we want to change something about implementation later.
    // (in which case the hook signature will probably also change)
    
    const { colorBy } = useChart() as VivMdvReact; //<<< todo better types
    
    const channelX = config.channel; //don't do this... use the proper image_properties
    // and make it populate the channel state if necessary.
    const stats = useChannelStats(ome, channelX);
    const contrastLimits = [stats ? stats.contrastLimits : [0, 1]];

    const data = useFilteredIndices(); // consider 'gray out' vs 'hide' for filtered points...
    const [cx, cy] = useParamColumns();
    const scatterplotLayer = new ScatterplotLayer({
        id: id+'scatter-react',
        data,
        opacity: 0.1,
        radiusScale: 1,
        getFillColor: colorBy ?? [0, 200, 200],
        getRadius: 1,
        getPosition: (i, {target}) => {
            // how slow is a callback function here?
            // given that we need to draw in data from multiple sources...
            // ... unless we want to do the work on a worker thread,
            // I don't think there's a significantly more efficient way
            target[0] = cx.data[i];
            target[1] = cy.data[i];
            target[2] = 0;
            return target as unknown as Float32Array; // 🤮 deck.gl types are wrong AFAICT
        },
        updateTriggers: {
            getFillColor: colorBy,
        }
    });
    if (!ome) return <div>Loading...</div>; // todo suspense.
    // note: higher-level VivViewer components remove some control over how we use layers
    // and some other props (e.g. if we wanted to override cursor style for lasso mode, we'd be SOL).
    // ... that may mean that the specific benefits of React over other frameworks may be less relevant.
    // TODO get viv working in popouts (not a react thing - happens elsewhere
    // - probably need to handle lost gl context)
    const detailView = new DetailView({
        id: id+'detail-react',
        snapScaleBar: true,
        width, height
    });
    const layerConfig = {
        loader: ome.data,
        selections: [{ z: 0, c: channelX, t: 0 }],
        contrastLimits,
        extensions: [new ColorPaletteExtension()], //garbage collection?
        colors: [[255, 255, 255]],
        channelsVisible: [true],
    }
    const detailLayers = detailView.getLayers({
        viewStates: [],
        props: layerConfig
    }); //includes a ScaleBarLayer... but I'm not seeing it. Also... two xr-layer passes? (looking at spector draw calls)
    return (
        <> 
            {`view id '${id}', channel #${channelX}: '${ome?.metadata?.Pixels?.Channels[channelX]?.Name}'`}
            <br />
            {contrastLimits[0][0]}-{contrastLimits[0][1]}
            <br />
            {data.length} points
        <DeckGL id={id + 'deck'} 
        views={[detailView.getDeckGlView()]}
        initialViewState={getDefaultInitialViewState(ome.data, {width, height})}
        style={{
            zIndex: '-1',
        }}
        // viewState={viewState}
        controller={true}
        layers={[detailLayers, scatterplotLayer]}>
        </DeckGL>
        </>
    );
});

/// some type stuff... common bits should probably move elsewhere...
// even if the type isn't constrained beyond 'string' (which potentially in some cases it could be)
// this kind of thing can potentially help DX and encourage appropriate values to be passed.
type ColumnName = string; 
export type ScatterPlotConfig = {
    radius: number,
    color_by: ColumnName,
    color_legend: {
        display: boolean,
        // todo: add more options here...
    },
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
} & BaseConfig & ScatterPlotConfig;

export type VivMdvReactConfig = 
BaseConfig & ScatterPlotConfig & (
    { type: 'VivMdvReact', imageURL: string, overviewOn: boolean, image_properties: ChannelsState } 
    | VivRoiConfig
) & { channel: number };
export type VivMDVReact = VivMdvReact;
class VivMdvReact extends BaseReactChart<VivMdvReactConfig> {
    colorDialog: any;
    constructor(dataStore, div, config: VivMdvReactConfig) {
        // todo better default config
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
    // next up... colorBy (for the scatterplot)...
    // I'd like to be able to have more of a layer-based thing...
    // I'd also quite like not to have these classes have non-trivial amounts of code in them...
    // For now, what I should do is a thing that bridges the existing colorBy stuff
    colorBy?: (i: number) => [r: number, g: number, b: number];
    colorByColumn(col?: ColumnName) {
        if (!col) return this.colorByDefault();
        this.config.color_by = col; //test this...
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
            // no longer using PictureInPictureViewer - would need to re-implement to some extent
            // in order for combined viewer to work.
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

BaseChart.types["VivMdvReact"] = {
    "class": VivMdvReact,
    name: "VivMdvReact",
    methodsUsingColumns: ["colorByColumn"], //consider annotations / different approach...
    params: [
        {
            type: "number",
            name: "x",
        },
        {
            type: "number",
            name: "y",
        }
    ],
    extra_controls: (ds) => {
        return [
            {
                //todo use `regions` config... or make a different version that uses that...
                type: "string",
                name: "imageURL",
                label: "ome.tiff URL",
                defaultVal: ""
            }
        ]
    },
};

BaseChart.types["VivMdvRegionReact"] = {
    ...BaseChart.types["viv_scatter_plot"],
    "class": VivMdvReact,
    name: "Viv Scatter Plot (react)",
}

export type VivMdvReactType = typeof VivMdvReact;
// we rely on the side-effect of this import to register the chart type
///- will register otherwise, but this tricks HMR into working...
export default 42;