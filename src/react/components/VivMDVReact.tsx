import { useEffect, useMemo, useState } from "react";
import BaseChart from "../../charts/BaseChart";
import { createRoot } from "react-dom/client";
import { ScatterplotLayer } from "deck.gl/typed";
import { PictureInPictureViewer } from "@hms-dbmi/viv";
import { useChannelStats, useChartSize, useConfigItem, useDataModel, useOmeTiff, useParamColumns } from "../hooks";

let ID = 0;
function ReactTest({ parent }: { parent: VivMdvReact }) {
    const [width, height] = useChartSize(parent);
    const ome = useOmeTiff(parent.config.imageURL);
    
    const [id] = useState(ID++);
    const channelX = useConfigItem(parent, 'channel');
    
    const stats = useChannelStats(ome, channelX);
    const contrastLimits = stats ? stats.contrastLimits : [0, 1];
    
    const { data } = useDataModel(parent);
    const [cx, cy] = useParamColumns(parent);
    const scatterplotLayer = new ScatterplotLayer({
        id: id+'scatter-react',
        data,
        opacity: 0.5,
        radiusScale: 1000,
        getFillColor: [255, 200, 200],
        getRadius: 1,
        getPosition: (i, {target}) => {
            // how slow is a callback function here?
            // given that we need to draw in data from multiple sources...
            // ... unless we want to do the work on a worker thread,
            // I don't think there's a more efficient way
            target[0] = 2*cx.data[i];
            target[1] = 2*cy.data[i];
            target[2] = 0;
            return target as unknown as Float32Array; // 🤮 deck.gl types are wrong AFAICT
        }
    });
    if (!ome) return <div>Loading...</div>; // todo suspense.

    return (
        <> 
            <PictureInPictureViewer 
                contrastLimits={[contrastLimits]}
                colors={[[255, 255, 255]]}
                channelsVisible={[true]}
                loader={ome?.data}
                selections={[{ z: 0, c: channelX, t: 0 }]}
                snapScaleBar={true}
                overview={undefined}
                overviewOn={false}
                // well and good, but why is aspect ratio wrong?
                height={height} width={width}
                deckProps={{
                    layers: [scatterplotLayer]
                }}
            />
            {ome?.metadata?.Pixels?.Channels[channelX]?.Name}
            <br />
            {contrastLimits[0]}-{contrastLimits[1]}
        {/* <DeckGL id={id + 'deck'} layers={[scatterplotLayer]}>
        </DeckGL> */}
        </>
    );
}

class VivMdvReact extends BaseChart {
    constructor(dataStore, div, config) {
        if (!config.channel) config.channel = 0;
        super(dataStore, div, config);
        createRoot(this.contentDiv).render(<ReactTest parent={this} />);
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
    remove(): void {
        super.remove();
        // make sure dim and anything else relevant is removed...
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