import type BaseChart from "@/charts/BaseChart";
import type { FieldSpec } from "@/react/components/ColumnSelectionComponent";
import type { IReactionDisposer } from "mobx";

/**
 * Apply the {@link loadColumnData} decorator to a method on a chart class manually at runtime, rather than using the `@loadColumnData` syntax
 * which as of this writing is only working with vanilla ts code - not js or react.
 */
export default function decorateColumnMethod<T extends BaseChart>(method: string, chart: T, dataSource: string) {
    // could try to make `method: keyof T`, but that would require a more type juggling...
    chart[method] = loadColumnData(chart[method], { target: chart });
}

// also consider configEntriesUsingColumns...
// perhaps we can have computed getters and custom setters for the columns in the config

/**
 * Decorator for a method that requires column data to be loaded.
 * 
 * Methods decorated with this will use the first argument as a column id or array of column ids.
 * 
 * Anything specified as a virtual column that may change at runtime should cause the original method to be called again
 * as a `reaction`.
 */
export function loadColumnData<This extends BaseChart, Args extends any[], Return>(
    target: (this: This, ...args: Args) => void, // we could probably type this to have first argument specified...
    context: ClassMethodDecoratorContext<This, (this: This, ...args: Args) => void> | { target: This },
) {
    let disposer: IReactionDisposer | null = null;
    function replacementMethod(this: This, ...args: any[]) {
        const dataSource = this.dataSource.name;
        // if we have a special value indicating live data, we can do something, perhaps with chart.mobxAutorun(), here...
        // we want to make sure that
        // - the result of the live column change will be applied ✅
        // - the value stored in the config will be the special value, rather than the result of the live column change
        // - reaction disposes when the live column is no longer needed ✅ (to be tested)
        // there could be anything happening inside the method... but perhaps if we know that the config is always a
        // proxied mobx object, we can somehow generically detect changes and re-patch???
        // perhaps safer to have a copy of the config that will be used in chart.getConfig(), which has the 'special value',
        // but then when the decorated method does things to its config, it won't be reflected in the 'special value' config that we save
        // this is really a 'config vs state' thing.
        if (args[0] == null) {
            // columns not needed
            target.call(this, ...args);
        } else {
            // this should only happen if there are any live columns...
            // and we should dispose if replacementMethod is called again
            if (disposer) {
                disposer();
                this.reactionDisposers = this.reactionDisposers.filter(v => v !== disposer);
            }
            const multiCol = Array.isArray(args[0]);
            const colsOriginal = multiCol ? args[0] : [args[0]];
            args[0] = colsOriginal;
            disposer = this.mobxAutorun(() => {
                const cols = args[0].flatMap(
                    (v: FieldSpec) => (typeof v === "string" ? v : v.fields)
                ).filter(v => v);
                const a = multiCol ? cols : cols[0];
                const newArgs = [a, ...args.slice(1)];
                const cm = window.mdv.chartManager;
                console.log('loading columns:', cols);
                cm._getColumnsThen(dataSource, cols, () => {
                    console.log('calling', target.name, 'with', newArgs);
                    target.call(this, ...newArgs);
                });
            });
        }
    }

    return replacementMethod;
}