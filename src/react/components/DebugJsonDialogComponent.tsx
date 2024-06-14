//import BaseChart from "@/charts/BaseChart";
import JsonView from 'react18-json-view'
import 'react18-json-view/src/style.css'
// If dark mode is needed, import `dark.css`.
import 'react18-json-view/src/dark.css'
import { useState } from 'react';

type JSONObject = { [key: string]: any };

/**
 * given a filter string and a json object, return a new json object that is a subset of the input json
 * matching the filter should mean that a given node is included in the output
 * if all (case insensitive) elements of filterArray appear somewhere in the node or its path.
 * credit ChatGPT for the code (copilot failed) */
function filterJSON(obj: JSONObject, filter: string): JSONObject {
    if (!filter) return obj;
    const filterArray = filter.toLowerCase().split(' ');

    function matchesFilter(value: string, path: string[]): boolean {
        const combinedString = path.concat(value).join(' ').toLowerCase();
        return filterArray.every(f => combinedString.includes(f));
    }

    function recursiveFilter(current: JSONObject, path: string[]): JSONObject | null {
        if (typeof current === 'string') {
            return matchesFilter(current, path) ? current : null;
        } else if (typeof current === 'object' && current !== null) {
            const filteredObj: JSONObject = {};
            let childMatch = false;

            for (const key in current) {
                const result = recursiveFilter(current[key], path.concat(key));
                if (result !== null) {
                    filteredObj[key] = result;
                    childMatch = true;
                }
            }
            if (childMatch || matchesFilter('', path)) {
                return filteredObj;
            }
        }
        return null;
    }

    return recursiveFilter(obj, []) || {};
}


export default function ({json, header}: {json: any, header?: string}) {
    const [filter, setFilter] = useState('');
    const filteredJson = filterJSON(json, filter);
    return (
        <div className='max-h-[90vh] overflow-auto'>
            {header && <h2>{header}</h2>}
            <input type='text' value={filter} onChange={e => setFilter(e.target.value)} placeholder='Filter...'></input>
            <JsonView src={filteredJson} />
        </div>
    );
}