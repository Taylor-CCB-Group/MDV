import { useEffect, useState } from "react";
import BaseChart from "../../charts/BaseChart";
import { createRoot } from "react-dom/client";

let ID = 0;
function ReactTest({parent}) {
    const {dataStore} = parent;
    const [id] = useState(ID++);
    const [filterSize, setFilterSize] = useState(dataStore.filterSize);
    const [dim] = useState(dataStore.getDimension("catcol_dimension"));
    const [words, setWords] = useState(['hello', 'world']);
    const [text, setText] = useState(parent.config.text);
    useEffect(() => {
        dataStore.addListener(id, () => setFilterSize(dataStore.filterSize));
        parent.addListener("text", (type, data) => setText(data + ' ' + parent.config.wordSize));
        const colNameWords = parent.config.param[0];
        const colNameSize = parent.config.param[1];

        const wordsVals = dataStore.getColumnValues(colNameWords);
        const sizeVals = dataStore.getColumnValues(colNameSize);
        const wordCol = dataStore.columnIndex[colNameWords];
        const sizeCol = dataStore.columnIndex[colNameSize];

        // for (let i = 0; i < wordCol.data.length; i++) {
        //     const word = wordsVals[wordCol.data[i]];
        //     const size = sizeVals[sizeCol.data[i]];
        //     console.log(word, size || `no size for ${i}`);
        // }
        dataStore.getRowAsObject()
        setWords(wordsVals);
        console.log(sizeVals);

        if (dim && dim.getAverages) {
            dim.getAverages(data => {
                console.log(data);
            }, [colNameSize], {})
        } else {
            console.log("no averages");
        }

        // dim.getCategories(data => {
        //     const d = new Array(data.length);
        //     let maxCount = 1;
        //     for (let i = 0; i < data.length; i++) {
        //         d[i] = [data[i], i];
        //         maxCount = Math.max(maxCount, data[i]);
        //     }
        //     const vals = dataStore.getColumnValues(colNameWords);
        //     const col = d.map(v => vals[v[1]]);
        //     setWords(col);
        // }, colNameWords, {});
        return () => {
            dataStore.removeListener(id);
            parent.removeListener("text");
        }
    }, [dataStore, dim, text]);
    return (
        <div style={{padding: '0.3em', overflow: 'auto'}}>
        <h2>Words:</h2>
        <pre>{JSON.stringify(words, null, 2)}</pre>
        </div>
    );
}

class ReactWordCloudChart extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        createRoot(this.contentDiv).render(<ReactTest parent={this} />);
    }
    drawChart() {
        this._callListeners("text", "Hello World");
    }
    getSettings(): { type: string; label: string; current_value: any; func: (v: any) => void; }[] {
        const c = this.config;
        const settings = super.getSettings();
        return settings.concat([
            {
                type: "slider",
                label: "Word Size",
                current_value: c.wordSize || 100,
                min: 10,
                max: 100,
                func: x => {
                    c.wordSize = x;
                    this.drawChart();
                }
            }
        ]);
    }
    remove(): void {
        super.remove();
        // make sure dim and anything else relevant is removed...
    }
}

BaseChart.types["WordCloud2"] = {
    "class": ReactWordCloudChart,
    name: "WordCloud Chart",
    params: [
        {
            type: ["text", "multitext"],
            name: "text",
        },
        {
            type: ["text", "multitext"],
            name: "word size",
        }
    ],
    init: (config, dataSource, extraControls) => {
        config.wordSize = 20;
    }
};

// we rely on the side-effect of this import to register the chart type
export default 42;