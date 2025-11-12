import type { FeatureCollection } from '@turf/helpers';
import { useViewerStore } from "./components/avivatorish/state";
import type { DeckScatterConfig } from "./components/DeckScatterReactWrapper";
import { useConfig } from "./hooks";


/**
 * we have different ways of managing view state in viv vs deck - should be unified...
 * this is an agnostic shim for mobx vs zustand pending refactor.
 */
export function useViewState() {
    // we should do something nicer than try/catch here,
    // but longer term we'll probably have a single store/mobx state.
    try {
        return useViewerStore(state => state.viewState);
    } catch (e) {
        const config = useConfig<DeckScatterConfig>();
        if (!config.viewState) {
            throw new Error("No view state found in config, and no viewer store found.");
        }
        return config.viewState;
    }
}

export const getEmptyFeatureCollection = () => ({
    type: "FeatureCollection",
    features: []
} as FeatureCollection);
