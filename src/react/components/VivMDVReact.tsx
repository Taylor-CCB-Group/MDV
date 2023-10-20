import { useEffect, useMemo, useState } from "react";
import BaseChart from "../../charts/BaseChart";
import { createRoot } from "react-dom/client";
import DeckGL from '@deck.gl/react';
import { ScatterplotLayer } from "deck.gl/typed";
import { DataModel } from "../../table/DataModel";
import { PictureInPictureViewer, VivViewer } from "@hms-dbmi/viv";
import {
    loadOmeTiff,
    DetailView,
    DETAIL_VIEW_ID,
    getChannelStats
} from '@hms-dbmi/viv';
import { useChartSize, useDataModel, useOmeTiff } from "../hooks";

type TiffLoader = Awaited<ReturnType<typeof loadOmeTiff>>;

let ID = 0;
function ReactTest({ parent }: { parent: VivMdvReact }) {
    const { dataStore } = parent;
    
    const [width, height] = useChartSize(parent);
    const dataModel = useDataModel(parent);
    const loader = useOmeTiff(parent.config.imageURL);

    const [id] = useState(ID++);
    const [channel, setChannel] = useState(0); //not really more usable than editing the code...


    const { data } = dataModel;
    const { columnIndex } = dataStore;
    const cx = columnIndex[parent.config.param[0]];
    const cy = columnIndex[parent.config.param[1]];

    const scatterplotLayer = new ScatterplotLayer({
        id: id+'scatter',
        data,
        getPosition: (i, {target}) => {
            // how slow is a callback function here?
            // given that we need to draw in data from multiple sources...
            // ... unless we want to do the work on a worker thread,
            // I don't think there's a more efficient way
            target[0] = cx.data[i];
            target[1] = cy.data[i];
            target[2] = 0;
            return target as unknown as Float32Array; // 🤮 deck.gl types are wrong AFAICT
        }
    });
    if (!loader) return <div>Loading...</div>;
    console.log('loader', loader);
    // xxx: logging deckrenderer.js:55 opts.layers.length 
    //it's growing out of control with Cillian's image (which it doesn't in Avivator)
    // ^^^ I'd passed SizeX & SizeY where I should be passing size of component...


    //channel 20, 0-10000 is legible in test_viv
    //cillian /project/cillian/images/avivator/BKV3Yc.ome.tiff
    return (
        <> 
            <PictureInPictureViewer 
                contrastLimits={[0, 1000000]} 
                colors={[[100, 245, 200]]}
                channelsVisible={[true]}
                loader={loader?.data}
                selections={[{ z: 0, c: channel, t: 0 }]}
                overview={undefined}
                overviewOn={false}
                height={height} width={width}
                deckProps={{
                    layers: [scatterplotLayer]
                }}
            />
            <input type="number" value={channel} onChange={(e) => { setChannel(+e.target.value); }}
                min={0} max={loader?.metadata?.Pixels?.SizeC - 1}
            />
            {loader?.metadata?.Pixels?.Channels[channel]?.Name}
        {/* <DeckGL id={id + 'deck'} layers={[scatterplotLayer]}>
        </DeckGL> */}
        </>
    );
}

class VivMdvReact extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        createRoot(this.contentDiv).render(<ReactTest parent={this} />);
    }
    drawChart() {
        this._callListeners("text", "Hello World");
    }
    // getSettings(): { type: string; label: string; current_value: any; func: (v: any) => void; }[] {
    //     const c = this.config;
    //     const settings = super.getSettings();
    //     return settings.concat([
    //         {
    //             type: "slider",
    //             label: "Word Size",
    //             current_value: c.wordSize || 100,
    //             min: 10,
    //             max: 100,
    //             func: x => {
    //                 c.wordSize = x;
    //                 this.drawChart();
    //             }
    //         }
    //     ]);
    // }
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