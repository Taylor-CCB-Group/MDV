import BaseChart  from "../../BaseChart";
import { g } from "../../../lib/utils";
import type { BaseConfig } from "../../BaseChart";
import UCSCBrowserComponent from "./UCSCBrowserComponent";
import { BaseReactChart } from "../../../react/components/BaseReactChart";
import type DataStore from "../../../datastore/DataStore";
import { applyViewMargins } from "../genomicLocationUtils";
import {type GenomeLocation,GenomeViewMargins} from "../genomicLocationUtils";


//all should be required but have default values
//will zod take care of this
interface UCSCBrowserConfig extends BaseConfig {
    src?: string;
    location?: GenomeLocation;
    view_margins?: GenomeViewMargins;
    highlight_selected_region?: boolean;

}


class UCSCBrowser extends BaseReactChart<UCSCBrowserConfig> {
    constructor(dataStore: DataStore, div: string | HTMLDivElement, config: UCSCBrowserConfig) {
        //set defaults - needed to be done before super call
        //is there a better way to this - can you set defaults in Zod schema?
        config.view_margins = config.view_margins || {type:"fixed_length", value:1000};
        config.location = config.location || {chr:"chr1", start:1000000, end:1001000};
        config.highlight_selected_region = config.highlight_selected_region ?? false;
        config.src = config.src || "";
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
        const vm = c.view_margins;
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
            current_value: vm.value.toString(),
            func: (x) => {
                const cur = c.view_margins ?? { type: "fixed_length" as const, value: 1000 };
                let n = Number.parseInt(x);
                n = Number.isNaN(n)
                    ? cur.type === "percentage"
                        ? 20
                        : 1000
                    : n;
                c.view_margins = { type: cur.type, value: n };
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
            current_value: vm.type,
            func: (x) => {
                const cur = c.view_margins ?? { type: "fixed_length" as const, value: 1000 };
                c.view_margins = { type: x as "percentage" | "absolute" | "fixed_length", value: cur.value };
            },
        })
        ];
    }

}

BaseChart.types["ucsc_browser"] = {
    name: "UCSC Browser",
    class: UCSCBrowser,
    required:(ds)=>{
        return (ds.genome) ;
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
    }
};
export {UCSCBrowserConfig};
export default UCSCBrowser;
