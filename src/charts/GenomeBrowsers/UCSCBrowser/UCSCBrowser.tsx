import BaseChart  from "../../BaseChart";
import { g } from "../../../lib/utils";
import type { BaseConfig } from "../../BaseChart";
import UCSCBrowserComponent from "./UCSCBrowserComponent";
import { BaseReactChart } from "../../../react/components/BaseReactChart";
import type DataStore from "../../../datastore/DataStore";
import { applyViewMargins } from "../genomicLocationUtils";


interface UCSCBrowserLocation {
    chr: string;
    start: number;
    end: number;
}
interface UCSCBrowserViewMargins{
    type:"fixed_length" | "absolute" | "percentage";
    value:number;
}

//all should be required but have default values
//will zod take care of this
interface UCSCBrowserConfig extends BaseConfig {
    src?: string;
    location?: UCSCBrowserLocation;
    view_margins?: UCSCBrowserViewMargins;
    highlight_selected_region?: boolean;

}

function getLocation(location:UCSCBrowserLocation, vm:UCSCBrowserViewMargins) : UCSCBrowserLocation  {
    return applyViewMargins(location, vm);
}


class UCSCBrowser extends BaseReactChart<UCSCBrowserConfig> {
    constructor(dataStore: DataStore, div: string | HTMLDivElement, config: UCSCBrowserConfig) {
        //set defaults - needed to be done before super call
        //is there a better way to this - can you set defaults in Zod schema?
        config.view_margins = config.view_margins || {type:"fixed_length", value:1000};
        config.src = config.src || "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38";
        config.location = config.location || {chr:"chr1", start:1000000, end:1001000};
        config.highlight_selected_region = config.highlight_selected_region ?? false;
        super(dataStore, div, config, UCSCBrowserComponent);
        this.contentDiv.style.overflowY = "scroll";
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        // Initialize view_margins if not present
        if (!c.view_margins) {
            c.view_margins = {type:"fixed_length", value:1000};
        }
        return [
            ...settings,
            g({
                type: "text",
                current_value: c.src || "",
                label: "Session URL",
                func: (x: string) => {
                    c.src = x;
                }
            }),
            g({
                label: "Highlight selected region",
                type: "check",
                current_value: c.highlight_selected_region || false,
                func: (x) => {
                    c.highlight_selected_region = x;
                },
            }),
        g(
            {
            label: "View margin length",
            type: "text",
            current_value: c.view_margins.value.toString(),
            func: (x) => {
                let n = Number.parseInt(x);
                n = Number.isNaN(n)
                    ? c.view_margins.type === "percentage"
                        ? 20
                        : 1000
                    : n;
                c.view_margins = { type: c.view_margins.type, value: n };
            },
        }

        ),
        g({
            label: "View margin type",
            type: "radiobuttons",
            choices: [
                ["% of feature length", "percentage"],
                ["absolute (bp)", "absolute"],
                ["Fixed Length (bp)", "fixed_length"],
            ],
            current_value: c.view_margins.type,
            func: (x) => {
                c.view_margins = { type: x as "percentage" | "absolute" | "fixed_length", value: c.view_margins.value };
            },
        })
        ];
    }

}

BaseChart.types["ucsc_browser"] = {
    name: "UCSC Browser",
    class: UCSCBrowser,
    //params are defined by the genome's genomic_location columns
    params: [],
    required:(ds)=>{
        return (ds.genome?.genomic_location || ds.genome?.svs) ;
    },
    extra_controls: (ds) => {
        return [
            {
                type: "text",
                name: "url",
                label: "URL",
            },
        ];
    },
    init(config, ds, extraControls) {
        config.src = extraControls.url;
        if (ds.genome?.genomic_location){
            const cols = ds.genome?.genomic_location?.columns;
            config.param = [cols["chr"], cols["start"], cols["end"]];
        //svs start and end may be on different chromosome
        } else if (ds.genome?.svs){
            const cols = ds.genome?.svs?.sv_columns;
            config.param = [cols["chr1"], cols["pos1"], cols["pos2"], cols["chr2"]];
        }
    }
};
export {getLocation, UCSCBrowserConfig, UCSCBrowserLocation, UCSCBrowserViewMargins};
export default UCSCBrowser;