import BaseChart from "../../charts/BaseChart";
import DeckGL, { OrthographicView, ScatterplotLayer } from "deck.gl/typed";
import { PictureInPictureViewer } from "@hms-dbmi/viv";
import { useChannelStats, useChartID, useChartSize, useDataModel, useOmeTiff, useParamColumns } from "../hooks";
import { useMemo, useState } from "react";
import { BaseConfig, BaseReactChart } from "./BaseReactChart";
import { observer } from "mobx-react-lite";

// we need to add observer here, not just in BaseReactChart, for HMR to work.
const ReactTest = observer(({ parent }: { parent: VivMdvReact }) => {
    const [width, height] = useChartSize(parent);
    const ome = useOmeTiff(parent.config.imageURL);
    const view = useMemo(()=>new OrthographicView({}), []);
    // const {id} = parent.config; //there should always be a persistent, reliable id here.
    const id = useChartID(parent); // hook in case we want to change something about implementation later.
    // (in which case the hook signature will probably also change)
    
    const channelX = parent.config.channel;//useConfigItem(parent, 'channel');
    
    const stats = useChannelStats(ome, channelX);
    const contrastLimits = stats ? stats.contrastLimits : [0, 1];

    const [viewState, setViewState] = useState<any>();
    
    const { data } = useDataModel(parent);
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
            target[0] = 2*cx.data[i];
            target[1] = 2*cy.data[i];
            target[2] = 0;
            return target as unknown as Float32Array; // 🤮 deck.gl types are wrong AFAICT
        }
    });
    if (!ome) return <div>Loading...</div>; // todo suspense.

    return (
        <> 
            {true && <PictureInPictureViewer 
                contrastLimits={[contrastLimits]}
                colors={[[255, 255, 255]]}
                channelsVisible={[true]}
                loader={ome?.data}
                selections={[{ z: 0, c: channelX, t: 0 }]}
                snapScaleBar={true}
                overview={undefined}
                overviewOn={false}
                height={height} width={width}
                onViewStateChange={(viewState) => {
                    // setViewState(viewState.viewState);//infinitely recursive
                }}
                deckProps={{
                    // how to get the appropriate view so that scatterplotLayer is visible?
                    // can't get a ref to PIPViewer?
                    views: [view],
                    layers: [scatterplotLayer],
                    style: {
                        opacity: 1,
                        zIndex: -10,
                        // display: 'none'
                    }
                }}
            />}
            {`view id '${id}', channel #${channelX}: '${ome?.metadata?.Pixels?.Channels[channelX]?.Name}'`}
            <br />
            {contrastLimits[0]}-{contrastLimits[1]}
        {/* <DeckGL id={id + 'deck'} 
        views={[view]}
        initialViewState={{
            target: [0, 0, 0],
            zoom: -10,
        }}
        // viewState={viewState}
        // controller={true}
        layers={[scatterplotLayer]}>
        </DeckGL> */}
        </>
    );
});

type VivMdvReactConfig = { channel: number, imageURL: string } & BaseConfig;

class VivMdvReact extends BaseReactChart<VivMdvReactConfig> {
    constructor(dataStore, div, config: VivMdvReactConfig) {
        if (!config.channel) config.channel = 0;
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
            }
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