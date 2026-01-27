import * as React from "react";
import DeckGL from "@deck.gl/react";
import type { MjolnirEvent } from 'mjolnir.js';
// import { getVivId } from '@hms-dbmi/views';
// No need to use the ES6 or React variants.
import equal from "fast-deep-equal";
import { ScaleBarLayer } from "@hms-dbmi/viv";
import type { OrthographicViewState, OrbitViewState, DeckGLProps, PickingInfo } from "deck.gl";
import { rebindMouseEvents } from "@/lib/deckMonkeypatch";
import type { EditableGeoJsonLayer } from "@deck.gl-community/editable-layers";
export function getVivId(id: string) {
    return `-#${id}#`;
}
export type ViewState = (OrthographicViewState | OrbitViewState) & { id: string }; //questionable
// type ViewStates = { [id: string]: ViewState } & { find: any }; //highly questionable - should look at runtime/docs...
export type ViewStates = ViewState[];
export type View = {id: string} & any; //placeholder
export type VivPickInfo = PickingInfo<any, any> & { tile: any }; //placeholder
const areViewStatesEqual = (viewState: ViewState, otherViewState?: ViewState) => {
    return (
        otherViewState === viewState ||
        (viewState?.zoom === otherViewState?.zoom &&
            //@ts-ignore CBA to discriminate between Orbit and Ortho viewStates
            viewState?.rotationX === otherViewState?.rotationX &&
            //@ts-ignore
            viewState?.rotationOrbit === otherViewState?.rotationOrbit &&
            equal(viewState?.target, otherViewState?.target))
    );
};

/**
 * @typedef viewStateChangeProps
 * @type {object}
 * @property {string} args.viewId
 * @property {object} args.viewState
 * @property {object} args.oldViewState
 * @ignore
 */

/**
 * @callback ViewStateChange
 * @param {viewStateChangeProps} args
 * @ignore
 */

/**
 * @callback Hover
 * @param {Object} info
 * @param {Object} event
 * @ignore
 */

/**
 * @callback HandleValue
 * @param {Array.<number>} valueArray pixel values for the image under the hover location
 * @ignore
 */

/**
 * @callback HandleCoordinate
 * @param {Object} coordnate The coordinate in the image from which the values are picked.
 * @ignore
 */

export type VivViewerWrapperProps = {
    views: View[];
    // viewStates: { [id: string]: ViewState & { id: string } }; //todo
    viewStates: ViewStates; //Map<string, ViewState & { id: string }>;
    onViewStateChange?: any;
    onHover?: (info: PickingInfo, event: MjolnirEvent) => void;
    hoverHooks?: any;
    layerProps: any;
    deckProps?: DeckGLProps;
    randomize?: boolean;
    useDevicePixels?: boolean;
    outerContainer?: HTMLElement;
    selectionLayer?: EditableGeoJsonLayer;
};
export type VivViewerWrapperState = {
    viewStates: any;
    //@ts-expect-error --- issue with new deck.gl - check up on this / review viv PR? ---
    deckRef?: React.MutableRefObject<DeckGL>;
    outerContainer?: Element;
};
/**
 * @typedef HoverHooks
 * @type {object}
 * @property {HandleValue} handleValue
 * @property {HandleCoordinate} handleCoordinate
 * @ignore
 */
class MDVivViewerWrapper extends React.PureComponent<
    VivViewerWrapperProps,
    VivViewerWrapperState
> {
    constructor(props: VivViewerWrapperProps) {
        console.log(
            "using custom VivViewerWrapper via MDVivViewer - now necessary for fixing mouse events on popouts",
        );
        super(props);
        this.state = {
            viewStates: {},
            deckRef: React.createRef(),
            outerContainer: props.outerContainer,
        };
        const { viewStates } = this.state;
        const { views, viewStates: initialViewStates } = this.props;
        views.forEach((view: any) => {
            viewStates[view.id] = view.filterViewState({
                viewState: initialViewStates.find((v: any) => v.id === view.id),
            });
        });
        this._onViewStateChange = this._onViewStateChange.bind(this);
        this.layerFilter = this.layerFilter.bind(this);
        this.onHover = this.onHover.bind(this);
    }

    /**
     * This prevents only the `draw` call of a layer from firing,
     * but not other layer lifecycle methods.  Nonetheless, it is
     * still useful.
     * @param {object} args
     * @param {object} args.layer Layer being updated.
     * @param {object} args.viewport Viewport being updated.
     * @returns {boolean} Whether or not this layer should be drawn in this viewport.
     */
    // eslint-disable-next-line class-methods-use-this
    layerFilter({ layer, viewport }: any) {
        //return true; // for testing whether viv id matching is an issue
        return layer.id.includes(getVivId(viewport.id));
    }

    /**
     * This updates the viewState as a callback to the viewport changing in DeckGL
     * (hence the need for storing viewState in state).
     */
    _onViewStateChange({ viewId, viewState, interactionState, oldViewState }: any) {
        // Save the view state and trigger rerender.
        const { views, onViewStateChange } = this.props;
        // eslint-disable-next-line no-param-reassign
        viewState =
            onViewStateChange?.({
                viewId,
                viewState,
                interactionState,
                oldViewState,
            }) || viewState;
        this.setState((prevState) => {
            const viewStates: any = {};
            views.forEach((view) => {
                const currentViewState = prevState.viewStates[view.id];
                viewStates[view.id] = view.filterViewState({
                    viewState: { ...viewState, id: viewId },
                    oldViewState,
                    currentViewState,
                });
            });
            return { viewStates };
        });
        return viewState;
    }


    /**
     * prevent the deckMonkeypatch from double-binding etc.
     * as of anything related to that, this is dubious and liable to need attention in the future.
     */
    _cleanupMouseEvents?: () => void;
    componentDidUpdate(prevProps: VivViewerWrapperProps) {
        const { props } = this;
        const { views, outerContainer, selectionLayer } = props;

        if (
            outerContainer !== this.state.outerContainer &&
            this.state.deckRef?.current
        ) {
            try {
                this.setState({ outerContainer });
                const deck = this.state.deckRef.current.deck;
                this._cleanupMouseEvents?.();
                //this should be common with DeckScatterComponent - make a helper/hook...
                this._cleanupMouseEvents = rebindMouseEvents(deck, selectionLayer);
            } catch (e) {
                console.error(
                    "attempt to reset deck eventManager element failed",
                    e,
                );
            }
        }
        

        // Only update state if the previous viewState prop does not match the current one
        // so that people can update viewState
        // eslint-disable-next-line react/destructuring-assignment
        const viewStates = { ...this.state.viewStates };
        let anyChanged = false;
        views.forEach((view) => {
            const currViewState = props.viewStates?.find(
                (viewState) => viewState.id === view.id,
            );
            if (!currViewState) {
                return;
            }
            const prevViewState = prevProps.viewStates?.find(
                (viewState) => viewState.id === view.id,
            );
            if (areViewStatesEqual(currViewState, prevViewState)) {
                return;
            }
            anyChanged = true;
            const { height, width } = view;
            viewStates[view.id] = view.filterViewState({
                viewState: {
                    ...currViewState,
                    height,
                    width,
                    id: view.id,
                },
            });
        });
        if (anyChanged) {
            // eslint-disable-next-line react/no-did-update-set-state
            this.setState({ viewStates });
        }
    }
    componentWillUnmount(): void {
        this._cleanupMouseEvents?.();
        this._cleanupMouseEvents = undefined;
    }

    /**
     * This updates the viewStates' height and width with the newest height and
     * width on any call where the viewStates changes (i.e resize events),
     * using the previous state (falling back on the view's initial state) for target x and y, zoom level etc.
     */
    static getDerivedStateFromProps(props: VivViewerWrapperProps, prevState: VivViewerWrapperState) {
        const { views, viewStates: viewStatesProps } = props;
        // Update internal viewState on view changes as well as height and width changes.
        // Maybe we should add x/y too?
        if (
            views.some(
                (view) =>
                    !prevState.viewStates[view.id] ||
                    view.height !== prevState.viewStates[view.id].height ||
                    view.width !== prevState.viewStates[view.id].width,
            )
        ) {
            const viewStates: ViewStates = {} as any; //nonono
            views.forEach((view) => {
                const { height, width } = view;
                const currentViewState = prevState.viewStates[view.id];
                viewStates[view.id] = view.filterViewState({
                    viewState: {
                        ...(currentViewState ||
                            viewStatesProps.find((v: any) => v.id === view.id)),
                        height,
                        width,
                        id: view.id,
                    },
                });
            });
            return { viewStates };
        }
        return prevState;
    }

    // eslint-disable-next-line consistent-return
    onHover(info: VivPickInfo, event: MjolnirEvent) {
        const { tile, coordinate, sourceLayer: layer } = info;
        const { onHover, hoverHooks } = this.props;
        if (onHover) {
            onHover(info, event);
        }
        if (!hoverHooks || !coordinate || !layer) {
            return null;
        }
        const { handleValue = () => {}, handleCoordnate = () => {} } =
            hoverHooks;
        let hoverData;
        // Tiled layer needs a custom layerZoomScale.
        if (layer.id.includes("Tiled")) {
            if (!tile?.content) {
                return null;
            }
            const {
                content,
                bbox,
                index: { z },
            } = tile;
            if (!content.data || !bbox) {
                return null;
            }
            const { data, width, height } = content;
            const { left, right, top, bottom } = bbox;
            const bounds = [
                left,
                data.height < layer.tileSize ? height : bottom,
                data.width < layer.tileSize ? width : right,
                top,
            ];
            if (!data) {
                return null;
            }
            // The zoomed out layer needs to use the fixed zoom at which it is rendered.
            const layerZoomScale = Math.max(1, 2 ** Math.round(-z));
            const dataCoords = [
                Math.floor((coordinate[0] - bounds[0]) / layerZoomScale),
                Math.floor((coordinate[1] - bounds[3]) / layerZoomScale),
            ];
            const coords = dataCoords[1] * width + dataCoords[0];
            hoverData = data.map((d: any) => d[coords]);
        } else {
            const { channelData } = layer.props;
            if (!channelData) {
                return null;
            }
            const { data, width, height } = channelData;
            if (!data || !width || !height) {
                return null;
            }
            const bounds = [0, height, width, 0];
            // Using floor means that as we zoom out, we are scaling by the zoom just passed, not the one coming.
            const { zoom } = layer.context.viewport;
            const layerZoomScale = Math.max(1, 2 ** Math.floor(-zoom));
            const dataCoords = [
                Math.floor((coordinate[0] - bounds[0]) / layerZoomScale),
                Math.floor((coordinate[1] - bounds[3]) / layerZoomScale),
            ];
            const coords = dataCoords[1] * width + dataCoords[0];
            hoverData = data.map((d: any) => d[coords]);
        }
        handleValue(hoverData);
        handleCoordnate(coordinate);
    }

    /**
     * This renders the layers in the DeckGL context.
     */
    _renderLayers() {
        const { onHover } = this;
        const { viewStates } = this.state;
        const { views, layerProps } = this.props;
        return views.map((view, i) =>
            view.getLayers({
                viewStates,
                props: {
                    ...layerProps[i],
                    onHover,
                },
            }),
        );
    }

    render() {
        /* eslint-disable react/destructuring-assignment */
        const {
            views,
            randomize,
            useDevicePixels = true,
            deckProps,
        } = this.props;
        const { viewStates } = this.state;
        const deckGLViews = views.map((view) => view.getDeckGlView());
        // DeckGL seems to use the first view more than the second for updates
        // so this forces it to use the others more evenly.  This isn't perfect,
        // but I am not sure what else to do.  The DeckGL render hooks don't help,
        // but maybe useEffect() would help?  I couldn't work it out as
        // The issue is that I'm not sure how React would distinguish between forced updates
        // from permuting the views array and "real" updates like zoom/pan.
        // I tried keeping a counter but I couldn't figure out resetting it
        // without triggering a re-render.
        if (randomize) {
            const random = Math.random();
            const holdFirstElement = deckGLViews[0];
            // weight has to go to 1.5 because we use Math.round().
            const randomWieghted = random * 1.49;
            const randomizedIndex = Math.round(
                randomWieghted * (views.length - 1),
            );
            deckGLViews[0] = deckGLViews[randomizedIndex];
            deckGLViews[randomizedIndex] = holdFirstElement;
        }
        // MDV: make sure the scalebar is on top of other layers
        // perhaps this should be in _renderLayers(), yolo.
        const vivLayers = this._renderLayers()[0];
        const scaleBarLayer = vivLayers.find((layer: any) => layer instanceof ScaleBarLayer);
        const otherLayers = vivLayers.filter((layer: any) => layer !== scaleBarLayer);
        const layers = deckProps?.layers === undefined
            ? [vivLayers]
            : [otherLayers, ...deckProps.layers];
        // XXX: including an undefined scaleBarLayer above causes some other layers to not render
        if (scaleBarLayer) {
            layers.push(scaleBarLayer);
        }
        return (
            <DeckGL
                // MDV: we mess with deck internals via ref to get eventManager to work in popouts
                ref={this.state.deckRef}
                // eslint-disable-next-line react/jsx-props-no-spreading
                {...(deckProps ?? {})}
                layerFilter={this.layerFilter}
                layers={layers}
                onViewStateChange={this._onViewStateChange}
                views={deckGLViews}
                viewState={viewStates}
                useDevicePixels={useDevicePixels}
                getCursor={({ isDragging }) => {
                    return isDragging ? "grabbing" : "crosshair";
                }}
            />
        );
    }
}

/**
 * This is a wrapper around the VivViewer component from @hms-dbmi/viv
 * *** THIS IS NOW ACTUALLY NECESSARY ***
 * to fix issues with mouse events in popouts.
 * In future, we may handle more interesting things here to do with layer rendering.
 */
export default (props: VivViewerWrapperProps) => (
    <MDVivViewerWrapper {...props} />
);
