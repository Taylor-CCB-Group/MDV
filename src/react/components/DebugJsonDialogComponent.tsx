//import BaseChart from "@/charts/BaseChart";
import JsonView from 'react18-json-view'
import 'react18-json-view/src/style.css'
// If dark mode is needed, import `dark.css`.
import 'react18-json-view/src/dark.css'
import { useState } from 'react';

export default function ({json, header}: {json: any, header?: string}) {
    const [filter, setFilter] = useState('');
    const filterArray = filter.toLowerCase().split(' ');
    // todo something better than this
    const filteredJson = JSON.parse(JSON.stringify(json, (key, value) => {    
        if (typeof value === 'string') {
            for (const filter of filterArray) {
                if (value.toLowerCase().includes(filter) || key.toLowerCase().includes(filter)) {
                    return value;
                }
            }
            return undefined;
        }
        return value;
    }));
    return (
        <div className='max-h-[90vh] overflow-auto'>
            {header && <h2>{header}</h2>}
            <input type='text' value={filter} onChange={e => setFilter(e.target.value)} placeholder='Filter...'></input>
            <JsonView src={filteredJson} collapsed={2}/>
        </div>
    );
}