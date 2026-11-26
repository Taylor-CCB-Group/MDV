import type { SimpleFeatureCollection } from "@deck.gl-community/editable-layers";
import { useContext } from "react";
import { useIsInVivContext, useOptionalViewerStore } from "./components/avivatorish/state";
import type { DeckScatterConfig } from "./components/DeckScatterReactWrapper";
import { ChartContext } from "./context";

/**
 * we have different ways of managing view state in viv vs deck - should be unified...
 * this is an agnostic shim for mobx vs zustand pending refactor.
 */
export function useViewState() {
    // longer term we'll probably have a single store/mobx state.
    const inViv = useIsInVivContext();
    const vs = useOptionalViewerStore((state) => state.viewState);
    // read directly rather than via useConfig, which would throw for a Viv viewer outside a chart context
    const config: DeckScatterConfig | undefined = useContext(ChartContext)?.config;
    if (inViv) {
        // Return a default viewState when it hasn't been initialized yet
        // This can happen when a new chart is initialized before the image loader is ready
        // The viewState should be initialized from config in MainChart, but this provides a safety net
        return vs ?? { zoom: 0, target: [0, 0] };
    }
    if (!config?.viewState) {
        throw new Error("No view state found in config, and no viewer store found.");
    }
    return config.viewState;
}

export const getEmptyFeatureCollection = () =>
    ({
        type: "FeatureCollection",
        features: [],
    }) as SimpleFeatureCollection;
