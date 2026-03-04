//import BaseChart from "@/charts/BaseChart";
import JsonView from "react18-json-view";
import "react18-json-view/src/style.css";
// If dark mode is needed, import `dark.css`.
import "react18-json-view/src/dark.css";
import { useMemo, useState } from "react";
import { useDebounce } from "use-debounce";
import "../../utilities/css/JsonDialogStyles.css";
import { Accordion, AccordionDetails, AccordionSummary, Divider, Link, Typography } from "@mui/material";
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import { DEBUG_JSON_REPORT_MESSAGE, MDV_EMAIL } from "@/utilities/constants";

type JSONObject = { [key: string]: any };

/**
 * given a filter string and a json object, return a new json object that is a subset of the input json
 * matching the filter should mean that a given node is included in the output
 * if all (case insensitive) elements of filterArray appear somewhere in the node or its path.
 * credit ChatGPT for the code (copilot failed) 
 * Adendum: the ChatGPT code did not meet the spec, so I had to fix it.
 * This had caused me to spend a lot longer debugging when I was trying to use it before I realised this was at fault.
 * */
function filterJSON(obj: JSONObject, filter: string): JSONObject {
    if (!filter) return obj;
    const filterArray = filter.toLowerCase().split(" ");

    function matchesFilter(value: string, path: string[]): boolean {
        const combinedString = path.concat(value).join(" ").toLowerCase();
        return filterArray.every((f) => combinedString.includes(f));
    }

    function recursiveFilter(
        current: JSONObject,
        path: string[],
    ): JSONObject | null {
        // we should return the whole subtree if the current node matches the filter
        if (matchesFilter("", path)) {
            return current;
        }
        if (typeof current === "string") {
            return matchesFilter(current, path) ? current : null;
        }
        //`typeof current === "object"` probably no longer needed
        if (typeof current === "object" && current !== null) {
            const filteredObj: JSONObject = {};
            let childMatch = false;

            for (const key in current) {
                const result = recursiveFilter(current[key], path.concat(key));
                if (result !== null) {
                    filteredObj[key] = result;
                    childMatch = true;
                }
            }
            if (childMatch || matchesFilter("", path)) {
                return filteredObj;
            }
        }
        return null;
    }

    return recursiveFilter(obj, []) || {};
}

export default function ({ json, header }: { json: any; header?: string }) {
    const [filter, setFilter] = useState("");
    const [debouncedFilter] = useDebounce(filter, 300); //why is this returning a tuple?
    const filteredJson = useMemo(
        () => filterJSON(json, debouncedFilter),
        [json, debouncedFilter]
    );

    return (
        <div className="overflow-auto p-4" style={{ userSelect: "text" }}>
            {header && <h2 className="text-lg font-semibold mb-4">{header}</h2>}
            <Accordion sx={{
                border: "1px solid var(--border_menu_bar_color)",
                mb: 2,
            }}>
                <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                    <Typography variant="h6">
                        Report an issue or a bug
                    </Typography>
                </AccordionSummary>
                <Divider />
                <AccordionDetails sx={{ mt: 1 }}>
                    <Typography>
                        {DEBUG_JSON_REPORT_MESSAGE[0]}
                        {' '}
                        <Link 
                            href={`mailto:${MDV_EMAIL}?subject=MDV%20Bug%20Report`} 
                            sx={{ 
                                textDecoration: "none",
                                "&.MuiLink-root": { color: "info.main" }
                            }}
                        >
                            {MDV_EMAIL}
                        </Link>
                    </Typography>
                    <Typography sx={{mt: 1}}>
                        {DEBUG_JSON_REPORT_MESSAGE[1]}
                    </Typography>
                    <Typography sx={{mt: 1}}>
                        {DEBUG_JSON_REPORT_MESSAGE[2]}
                    </Typography>
                </AccordionDetails>
            </Accordion>
            <input
                type="text"
                value={filter}
                onChange={(e) => setFilter(e.target.value)}
                placeholder="Filter..."
                className="w-full mb-4 mt-4 p-2 border rounded-md dark:bg-gray-700 dark:border-gray-600"
            />
            <JsonView 
                src={filteredJson} 
                style={{
                    fontSize: '0.875rem',
                    fontFamily: 'ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace',
                    padding: "10px"
                }}
                collapsed={1}
            />
        </div>
    );
}
