import { makeAutoObservable } from "mobx";
import BaseChart from "../../charts/BaseChart";
import { createRoot } from "react-dom/client";
import type DataStore from '../../datastore/DataStore'
import { ChartProvider } from "../context";

function Fallback() {
    return <>
        <h2>BaseReactChart</h2>
        <p>BaseReactChart is a base class for charts that use React.</p>
        <p>As of this writing, in order to implement a new React-based chart, 
            make a class that extends BaseReactChart, and make it call the super constructor
            with an argument for the (mobx observable) function that does the react rendering.</p>
        <p>See `react.md` for more detail...</p>
    </>
}

export type BaseConfig = { 
    id: string, 
    size: [x: number, y: number], 
    title: string,
    legend: string,
    type: string,
    param: string[] | string,
};
type TComponent<T extends BaseConfig> = () => JSX.Element;

/**
 * Base class for charts that use React.
 * 
 * Components should be functional, but will require a very thin wrapper class extending this, such
 * that any BaseChart methods called from outside will appropriately affect the component state.
 * 
 * Internally, this class uses mobx to make the config object observable, and creates a React root with the given
 * `ReactComponentFunction`. Use of MobX is intended to be hidden by this abstraction - but if we change to a
 * different state management system, it will have implications for downstream components that expect to react to
 * changes in the config object.
 * 
 * We may also want to consider a different approach to the React root, i.e. a single root with portals for each chart, in
 * which case it should be handled in this class and should not (hopefully) require child classes/components to change.
 */
export abstract class BaseReactChart<T> extends BaseChart {
    declare config: T & BaseConfig;
    useMobx = true;
    root: ReturnType<typeof createRoot>;
    protected constructor(dataStore: DataStore, div: string | HTMLDivElement, config: T & BaseConfig, ReactComponentFunction: TComponent<T & BaseConfig> = Fallback) {
        super(dataStore, div, config);
        config = this.config; //original config will be copied by super, before doing things like adding id to it...
        makeAutoObservable(config);
        Object.defineProperty(this, 'config', {
            get: () => config,
            set: (v) => {
                config = v;
                makeAutoObservable(config);
            }
        });
        // note: a previous version of this used makeObservable for keeping track of onDataFiltered...
        // that worked, with extra extraneous number that changed to be observed by the hook...
        // What I have now done is change DataStore to be observable, and added a method for getting filtered indices
        // in a way that can be shared by different charts (react or otherwise).
        
        // note: observer() needs to be applied at all levels of the component tree that need to react to changes in 
        // any mobx state, so we can't just hide it in the base class.
        // (although maybe we could design a hook that hides it?)
        // const Observed = observer(ReactComponentFunction);
        this.root = createRoot(this.contentDiv);
        this.root.render((
            <ChartProvider chart={this}>
                <ReactComponentFunction />
            </ChartProvider>
        ));
    }
    remove(): void {
        // make sure dim and anything else relevant is removed...
        // **is there any React teardown we should be considering?**
        this.root.unmount();
        super.remove();
    }
}
