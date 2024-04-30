import BaseChart from "@/charts/BaseChart";
import JsonView from 'react18-json-view'
import 'react18-json-view/src/style.css'
// If dark mode is needed, import `dark.css`.
import 'react18-json-view/src/dark.css'

export default function ({chart}: {chart: BaseChart}) {
    return (
        <div>
            <h1>Debug Chart</h1>
            {/* <pre>
                {JSON.stringify(chart.config, null, 2)}
            </pre> */}
            <JsonView src={chart.config} />
        </div>
    );
}