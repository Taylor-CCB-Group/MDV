import { action, makeAutoObservable, makeObservable, observable } from "mobx";
import { observer } from "mobx-react-lite";
import BaseChart from "../../charts/BaseChart";
import { createRoot } from "react-dom/client";
import Dimension from "../../datastore/Dimension";

function Fallback() {
    return <>
        <h2>BaseReactChart</h2>
        <p>BaseReactChart is a base class for charts that use React.</p>
        <p>As of this writing, in order to implement a new React-based chart, 
            make a class that extends BaseReactChart, and make it call the super constructor
            with an argument for the (mobx observable) function that does the react rendering.</p>
        <p>TODO: documentation...</p>
    </>
}

export type BaseConfig = { id: string, size: [x: number, y: number] };
type TComponent<T extends BaseConfig> = (props: {parent: BaseReactChart<T>}) => JSX.Element;

/**
 * Base class for charts that use React.
 * 
 * Components should be functional, but will require a very thin wrapper class extending this, such
 * that any BaseChart methods called from outside will appropriately effect the component state.
 * 
 * Internally, this class uses mobx to make the config object observable, and creates a React root with the given
 * `ReactComponentFunction`. Use of MobX is hidden by this abstraction - but if we change to a different state management
 * system, it will have implications for downstream components that expect to react to changes in the config object.
 * 
 * We could probably abstract that away more fully; if we always access the config object through a hook (which might internally
 * use context API), we can probably save ourselves from even needing to pass 'parent' as a prop.
 * 
 * We may also want to consider a different approach to the React root, i.e. a single root with portals for each chart, in
 * which case it should be handled in this class and should not (hopefully) require child classes/components to change.
 */
export abstract class BaseReactChart<TConfig extends BaseConfig> extends BaseChart {
    declare config: TConfig;
    useMobx = true;
    dataFilter: Dimension | null = null;
    dataFilterSignal = 0;
    constructor(dataStore, div, config: TConfig, ReactComponentFunction: TComponent<TConfig> = Fallback) {
        super(dataStore, div, config);
        makeAutoObservable(config);
        Object.defineProperty(this, 'config', {
            get: () => config,
            set: (v) => {
                config = v;
                makeAutoObservable(config);
            }
        });
        makeObservable(this, {
            onDataFiltered: action,
            dataFilter: observable,
            dataFilterSignal: observable,
        })
        
        // note: applying observer here seems nice in that it means we can hide that implementation detail from child components,
        // but it stops HMR from working.
        // const Observed = observer(ReactComponentFunction);
        createRoot(this.contentDiv).render(<ReactComponentFunction parent={this} />);
    }
    onDataFiltered(dim: Dimension): void {
        super.onDataFiltered(dim);
        this.dataFilter = dim;
        this.dataFilterSignal++;
    }
    remove(): void {
        // make sure dim and anything else relevant is removed...
        // **is there any React teardown we should be considering?**
        super.remove();
    }
}
