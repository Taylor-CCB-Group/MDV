/** Compile-time smoke check that @spatialdata/* public API resolves. */
import type { RenderStack } from "@spatialdata/layers";
import type { SpatialCanvasViewer } from "@spatialdata/vis";

export type SpatialDataSmokeRenderStack = RenderStack;
export type SpatialDataSmokeViewer = typeof SpatialCanvasViewer;
