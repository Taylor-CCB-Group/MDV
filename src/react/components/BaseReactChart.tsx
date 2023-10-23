import { makeAutoObservable } from "mobx";
import { observer } from "mobx-react-lite";
import BaseChart from "../../charts/BaseChart";
import { createRoot } from "react-dom/client";

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

export type BaseConfig = { id: string, width: number, height: number };
type TComponent<T extends BaseConfig> = (props: {parent: BaseReactChart<T>}) => JSX.Element;

/**
 * Base class for charts that use React.
 * 
 * Internally, this class uses mobx to make the config object observable, and creates a React root with the given
 * `ReactComponentFunction`. Use of MobX is hidden by this abstraction - but if we change to a different state management
 * system, it will have implications for downstream components that expect to react to changes in the config object.
 */
export abstract class BaseReactChart<TConfig extends BaseConfig> extends BaseChart {
    declare config: TConfig;
    reactComponentFunction: TComponent<TConfig> = Fallback;
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
        const Observed = observer(ReactComponentFunction);
        createRoot(this.contentDiv).render(<Observed parent={this} />);
    }
}