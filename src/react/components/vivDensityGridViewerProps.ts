import type { OrbitViewState, OrthographicViewState } from "deck.gl";
import type { VivViewerWrapperProps } from "./avivatorish/MDVivViewer";

/** Props bundle for {@link MDVivViewer} in Viv scatter density-grid mode. */
export type VivDensityGridViewerProps = VivViewerWrapperProps & {
    onViewStateChange: (e: {
        viewId: string;
        viewState: OrthographicViewState | OrbitViewState;
    }) => void;
};
