import { g } from "@/lib/utils";

export type CrosshairConfig = {
    /** Mark the pointer with a horizontal and a vertical line in every cell of a grid; absent is on. */
    crosshair?: boolean | null;
    dimension?: "2d" | "3d";
};

/**
 * Whether a grid (of density fields, colours or a layer's fill) marks the pointer with a crosshair
 * in every cell: on unless turned off, and never in 3D. Single views and region grids have none.
 */
export function isCrosshairShown(config: CrosshairConfig) {
    return config.dimension !== "3d" && (config.crosshair ?? true);
}

/** The setting that turns the crosshair on and off. */
export function getCrosshairSetting(config: CrosshairConfig) {
    return g({
        type: "check",
        label: "Grid: mark the pointer in every cell",
        current_value: isCrosshairShown(config),
        func: (x) => {
            config.crosshair = x;
        },
    });
}

/**
 * Moves the crosshair to the pointer at (`x`, `y`) over `area`, or hides it when there is no
 * area. The position is set on `root` as an offset from the area's centre, so every area under
 * `root` with the same view (the cells of a grid) marks the same point. Written straight to the
 * style, as a re-render per pointer move would also set every prop of the deck again.
 */
export function placeCrosshair(root: HTMLElement, area: Element | undefined, x: number, y: number) {
    if (!area) {
        hideCrosshair(root);
        return;
    }
    const rect = area.getBoundingClientRect();
    root.style.setProperty("--mdv-crosshair-dx", `${x - rect.left - rect.width / 2}px`);
    root.style.setProperty("--mdv-crosshair-dy", `${y - rect.top - rect.height / 2}px`);
    root.dataset.crosshair = "on";
}

export function hideCrosshair(root: HTMLElement) {
    delete root.dataset.crosshair;
}
