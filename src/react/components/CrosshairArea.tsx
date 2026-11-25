import "./pointerCrosshair.css";

/**
 * The crosshair's lines, for an area that is positioned (see `placeCrosshair`). Only grids whose
 * cells share a view draw them: there the lines mark the same point in every cell, which a
 * single view, or cells with views of their own, have no use for.
 */
export function CrosshairLines() {
    return (
        <>
            <div className="mdv-crosshair" data-axis="x" />
            <div className="mdv-crosshair" data-axis="y" />
        </>
    );
}
