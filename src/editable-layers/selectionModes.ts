import {
    DrawPolygonMode,
    ModifyMode,
    // DrawRectangleMode,
    // TransformMode,
    TranslateMode,
    CompositeMode,
    DrawPolygonByDraggingMode,
    type ClickEvent,
    type FeatureCollection,
    type GeoJsonEditMode,
    type ModeProps,
    type SimpleFeatureCollection,
} from "@deck.gl-community/editable-layers";
import { DrawRectangleByDraggingMode } from "./deck-community-ish/draw-rectangle-by-dragging-mode";

/**
 * Selection geometry lives in the chart's data coordinates under an `OrthographicView`, so the
 * selection layer declares itself cartesian; editable-layers derives its edit-mode maths
 * (translate, midpoints etc.) from this.
 */
export const SELECTION_COORDINATE_SYSTEM = "cartesian";

//wtf do we need this adapter for and how can the types be so confusing/wrong?
//"Type 'ModifyMode' is not assignable to type 'GeoJsonEditMode'." etc
//even though `class ModifyMode extends GeoJsonEditMode`
function compositeModes(modes: unknown[]) {
    return modes as GeoJsonEditMode[];
}

/**
 * Upstream `CompositeMode` (as of editable-layers 9.4.5) forwards click / move / drag to its
 * modes but not double-click, so a double-click falls through to the no-op in `GeoJsonEditMode`.
 * That left the Polygon tool unable to finish on double-click: it only closed when a double-click
 * happened to be recognised as two clicks, with the second landing on the last vertex.
 */
class SelectionCompositeMode extends CompositeMode {
    handleDoubleClick(event: ClickEvent, props: ModeProps<FeatureCollection>) {
        this._coalesce((handler) => handler.handleDoubleClick(event, props));
    }
}

/**
 * Upstream finishes on any double-click, so one with fewer than three vertices (e.g. if a click
 * was lost) commits a zero-area polygon that replaces the current selection. Keep drawing instead.
 */
class DrawSelectionPolygonMode extends DrawPolygonMode {
    handleDoubleClick(event: ClickEvent, props: ModeProps<SimpleFeatureCollection>) {
        if (this.getClickSequence().length < 3) return;
        super.handleDoubleClick(event, props);
    }
}

export class EditMode extends CompositeMode {
    constructor() {
        super(compositeModes([
            new TranslateMode(),
            new ModifyMode(),
        ]));
    }

    /**
     * `ModifyMode` removes a vertex on click. A double-click is recognised as a single `dblclick`
     * rather than two clicks, so treat it the same way rather than silently ignoring it.
     */
    handleDoubleClick(event: ClickEvent, props: ModeProps<FeatureCollection>) {
        this.handleClick(event, props);
    }
}

export class RectangleMode extends SelectionCompositeMode {
    constructor() {
        super(compositeModes([
            // we pass in modeProps to the layer, not the edit mode.
            // new DrawRectangleMode(), //with modeConfig: { dragToDraw: true }
            //our version, which we shouldn't need but looks/works a bit different as of now
            new DrawRectangleByDraggingMode(),
            new TranslateMode(),
        ]));
    }
}

export class PolygonMode extends SelectionCompositeMode {
    constructor() {
        super(compositeModes([
            new DrawSelectionPolygonMode(),
            new TranslateMode(),
        ]));
    }
}

// todo - figure out weird conflict behaviour with this mode...
export class FreehandMode extends SelectionCompositeMode {
    constructor() {
        super(compositeModes([
            new DrawPolygonByDraggingMode(),
            new TranslateMode(),
        ]));
    }
}
