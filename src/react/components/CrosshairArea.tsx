import type { ReactNode } from "react";
import { hideCrosshair, placeCrosshair } from "./pointerCrosshair";
import "./pointerCrosshair.css";

/** The crosshair's lines, for an area that is positioned (see `placeCrosshair`). */
export function CrosshairLines() {
    return (
        <>
            <div className="mdv-crosshair" data-axis="x" />
            <div className="mdv-crosshair" data-axis="y" />
        </>
    );
}

/** An area filling its parent that, when `enabled`, marks the pointer over it with a crosshair. */
export function CrosshairArea({ enabled, children }: { enabled: boolean; children: ReactNode }) {
    return (
        // biome-ignore lint/a11y/noStaticElementInteractions: the handlers only draw the crosshair
        <div
            className="relative h-full w-full"
            onPointerMove={
                enabled
                    ? (event) => placeCrosshair(event.currentTarget, event.currentTarget, event.clientX, event.clientY)
                    : undefined
            }
            onPointerLeave={(event) => hideCrosshair(event.currentTarget)}
        >
            {children}
            {enabled && <CrosshairLines />}
        </div>
    );
}
