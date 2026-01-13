import { createMdvPortal } from "@/react/react_utils";
import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import type DataStore from "../../datastore/DataStore";
import { createEl } from "../../utilities/ElementsTyped";
import { ChartProvider } from "../context";

function Fallback() {
    return (
        <>
            <h2>BaseReactChart</h2>
            <p>BaseReactChart is a base class for charts that use React.</p>
            <p>
                As of this writing, in order to implement a new React-based
                chart, make a class that extends BaseReactChart, and make it
                call the super constructor with an argument for the (mobx
                observable) function that does the react rendering.
            </p>
            <p>See `react.md` for more detail...</p>
        </>
    );
}

export type TComponent<T extends BaseConfig> = () => React.JSX.Element;

/**
 * Base class for charts that use React.
 *
 * Components should be functional, but will require a very thin wrapper class extending this, such
 * that any BaseChart methods called from outside will appropriately affect the component state.
 *
 * Internally, this class uses mobx to make the config object observable, and creates a React root with the given
 * `ReactComponentFunction`. Use of MobX was originally intended to be hidden by this abstraction - but if we change to a
 * different state management system, it will have implications for downstream components that expect to react to
 * changes in the config object. Any component that needs to react to changes in the config object should be wrapped in
 * an `observer()` call.
 *
 * We may also want to consider a different approach to the React root, i.e. a single root with portals for each chart, in
 * which case it should be handled in this class and should not (hopefully) require child classes/components to change.
 */
export abstract class BaseReactChart<T extends BaseConfig> extends BaseChart<T> {
    // declare config: T & BaseConfig; // would be good to review this T & BaseConfig thing...
    declare popoutIcon: HTMLElement;
    // get dataSource(): DataSource {
    //     return window.mdv.chartManager.charts[this.config.id].dataSource;
    // }
    useMobx = true;
    root?: ReturnType<typeof createMdvPortal>;
    reactEl: HTMLDivElement;
    ComponentFn: TComponent<T>;
    protected constructor(
        dataStore: DataStore,
        div: string | HTMLDivElement,
        config: T,
        ReactComponentFunction: TComponent<T> = Fallback,
    ) {
        super(dataStore, div, config);
        //! for review - can we bypass this in react charts after a bit more evolution of the state management?
        //delete this.activeQueries.methodToConfigMap['setParams'];
        config = this.config; //original config will be copied by super, before doing things like adding id to it...
        // note: a previous version of this used makeObservable for keeping track of onDataFiltered...
        // that worked, with extra extraneous number that changed to be observed by the hook...
        // What I have now done is change DataStore to be observable, and added a method for getting filtered indices
        // in a way that can be shared by different charts (react or otherwise).

        // note: observer() needs to be applied at all levels of the component tree that need to react to changes in
        // any mobx state, so we can't just hide it in the base class.
        // (although maybe we could design a hook that hides it?)
        // const Observed = observer(ReactComponentFunction);
        this.reactEl = createEl(
            "div",
            { className: "react-chart" },
            this.contentDiv,
        ); //other things may still be added to contentDiv outside react (e.g. legend)
        this.ComponentFn = ReactComponentFunction;
        this.mountReact();
    }
    private mountReact() {
        const ReactComponentFunction = this.ComponentFn;
        this.root = createMdvPortal(
            //! todo sort out this chart type thing...
            //(apparently ok again now)
            <ChartProvider chart={this}>
                <ReactComponentFunction />
            </ChartProvider>,
            this.reactEl,
            this,
        );
    }
    //todo: implement this
    // getImage(callback: (imgCa: string) => void, type: "svg" | "png" = "png"): void {
    //     // "Charts should have a `getImage()` method that takes a function and a type, either 'svg' or 'png'.
    //     // The callback should pass a serialized svg string in the case of svg, or
    //     // a canvas element in the case of a png. (??)
    //     // if (type === "png") toPng(this.reactEl).then(callback);
    //     // if (type === "svg") toSvg(this.reactEl).then(callback);
    //     if (type === "png") {

    //     }
    // }
    remove(): void {
        // make sure dim and anything else relevant is removed...
        // **is there any React teardown we should be considering?**
        this.root?.unmount();
        super.remove();
        // this is handled by super.remove()
        // for (const disposer of this.reactionDisposers) {
        //     disposer();
        // }
    }
}
