import type { Chart } from "@/charts/charts";


export default function decorateColumnMethod<T extends Chart>(method: string, chart: T, dataSource: string) {
    const cm = window.mdv.chartManager;
    const newMethod = `_${method}`;
    chart[newMethod] = chart[method];
    console.log("decorating", method, chart.config);
    //if original method is called check whether column has data
    //first argument must be column(s) needed
    // biome-ignore lint/complexity/useArrowFunction: we are using `arguments` here, think we need `function` syntax
    chart[method] = function () {
        //column not needed
        if (arguments[0] == null) {
            chart[newMethod](...arguments);
        } else {
            const cols = Array.isArray(arguments[0])
                ? arguments[0]
                : [arguments[0]];
            // if we have a special value indicating live data, we can do something, perhaps with chart.mobxAutorun(), here...
            // we want to make sure that
            // - the result of the live column change will be applied
            // - the value stored in the config will be the special value, rather than the result of the live column change
            // - reaction disposes when the live column is no longer needed
            // there could be anything happening inside the method... but perhaps if we know that the config is always a
            // proxied mobx object, we can somehow generically detect changes and re-patch???
            // perhaps safer to have a copy of the config that will be used in chart.getConfig(), which has the 'special value',
            // but then when the decorated method does things to its config, it won't be reflected in the 'special value' config that we save
            // this is really a 'config vs state' thing.
            cm._getColumnsThen(dataSource, cols, () =>
                chart[newMethod](...arguments),
            );
        }
    };
}