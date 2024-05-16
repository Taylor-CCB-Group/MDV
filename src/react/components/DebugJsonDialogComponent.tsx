//import BaseChart from "@/charts/BaseChart";
import JsonView from 'react18-json-view'
import 'react18-json-view/src/style.css'
// If dark mode is needed, import `dark.css`.
import 'react18-json-view/src/dark.css'

export default function ({json, header}: {json: any, header?: string}) {
    return (
        <div className='max-h-[90vh] overflow-auto'>
            {header && <h2>{header}</h2>}
            <JsonView src={json} />
        </div>
    );
}