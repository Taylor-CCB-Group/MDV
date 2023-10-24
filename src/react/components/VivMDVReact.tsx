import BaseChart from "../../charts/BaseChart";
import DeckGL, { OrthographicView, ScatterplotLayer } from "deck.gl/typed";
import { ColorPaletteExtension, DetailView, getDefaultInitialViewState } from "@hms-dbmi/viv";
import { useChannelStats, useChartID, useChartSize, useFilteredIndices, useOmeTiff, useParamColumns } from "../hooks";
import { BaseConfig, BaseReactChart } from "./BaseReactChart";
import { observer } from "mobx-react-lite";

// we need to add observer here, not just in BaseReactChart, for HMR to work.
const ReactTest = observer(({ parent }: { parent: VivMdvReact }) => {
    const [width, height] = useChartSize(parent);
    const ome = useOmeTiff(parent.config.imageURL);
    // const {id} = parent.config; //there should always be a persistent, reliable id here.
    const id = useChartID(parent); // hook in case we want to change something about implementation later.
    // (in which case the hook signature will probably also change)
    
    const channelX = parent.config.channel;//useConfigItem(parent, 'channel');
    
    const stats = useChannelStats(ome, channelX);
    const contrastLimits = [stats ? stats.contrastLimits : [0, 1]];

    const data = useFilteredIndices(parent); //<< subject to revision
    const [cx, cy] = useParamColumns(parent);
    const scatterplotLayer = new ScatterplotLayer({
        id: id+'scatter-react',
        data,
        opacity: 0.5,
        radiusScale: 10,
        getFillColor: [0, 200, 200],
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
        }
    });
    if (!ome) return <div>Loading...</div>; // todo suspense.
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
    });
    const { SizeX, SizeY } = ome.metadata.Pixels;
    return (
        <> 
            {`view id '${id}', channel #${channelX}: '${ome?.metadata?.Pixels?.Channels[channelX]?.Name}'`}
            <br />
            {contrastLimits[0]}-{contrastLimits[1]}
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

type VivMdvReactConfig = { channel: number, imageURL: string, overviewOn: boolean } & BaseConfig;

class VivMdvReact extends BaseReactChart<VivMdvReactConfig> {
    constructor(dataStore, div, config: VivMdvReactConfig) {
        // todo better default config
        if (!config.channel) config.channel = 0;
        if (config.overviewOn === undefined) config.overviewOn = false;
        super(dataStore, div, config, ReactTest);
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
                type: "string",
                name: "imageURL",
                label: "ome.tiff URL",
                defaultVal: ""
            }
        ]
    },
};

// we rely on the side-effect of this import to register the chart type
///- will register otherwise, but this tricks HMR into working...
export default 42;