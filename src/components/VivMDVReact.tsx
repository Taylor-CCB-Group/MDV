import { useMemo, useState } from "react";
import BaseChart from "../charts/BaseChart";
import { createRoot } from "react-dom/client";
import DeckGL from '@deck.gl/react';
import { ScatterplotLayer } from "deck.gl/typed";
import { DataModel } from "../table/DataModel";


let ID = 0;
function ReactTest({ parent }: { parent: VivMdvReact }) {
    const { dataStore } = parent;
    const [id] = useState(ID++);

    const dataModel = useMemo(() => {
        const dataModel = new DataModel(dataStore, { autoUpdate: true });
        dataModel.setColumns(parent.config.param);
        dataModel.updateModel();
        return dataModel;
    }, [dataStore]);
    const { data } = dataModel;
    const { columnIndex } = dataStore;
    const cx = columnIndex[parent.config.param[0]];
    const cy = columnIndex[parent.config.param[1]];

    const scatterplotLayer = new ScatterplotLayer({
        id: id+'scatter',
        data,
        getPosition: (i, {target}) => {
            target[0] = cx.data[i];
            target[1] = cy.data[i];
            target[2] = 0;
            return target as unknown as Float32Array; // 🤮 deck.gl types are wrong AFAICT
        }
    })
    return (
        <>
        <DeckGL id={id + 'deck'} layers={[scatterplotLayer]}>
        </DeckGL>
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