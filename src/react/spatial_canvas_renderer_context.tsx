import type { SpatialData } from "@spatialdata/core";
import type {
    ElementsByType,
    useSpatialCanvasRenderer,
} from "@spatialdata/vis";
import { createContext, useContext, type PropsWithChildren } from "react";
import type { TableAssociation } from "./spatial_table_association";

export type SpatialCanvasRendererValue = {
    spatialData: SpatialData | null;
    coordinateSystem: string | null;
    availableElements: ElementsByType;
    getImageLayerLoadedData: ReturnType<typeof useSpatialCanvasRenderer>["getImageLayerLoadedData"];
    getLabelsLayerLoadedData: ReturnType<typeof useSpatialCanvasRenderer>["getLabelsLayerLoadedData"];
    getLayerLoadState: ReturnType<typeof useSpatialCanvasRenderer>["getLayerLoadState"];
    isLoading: boolean;
    isBlocking: boolean;
    inferTableAssociation: (elementKey?: string | null) => TableAssociation;
};

const SpatialCanvasRendererContext = createContext<SpatialCanvasRendererValue | null>(null);

export function SpatialCanvasRendererProvider({
    value,
    children,
}: PropsWithChildren<{ value: SpatialCanvasRendererValue }>) {
    return (
        <SpatialCanvasRendererContext.Provider value={value}>
            {children}
        </SpatialCanvasRendererContext.Provider>
    );
}

export function useSpatialCanvasRendererContext() {
    const value = useContext(SpatialCanvasRendererContext);
    if (!value) {
        throw new Error("useSpatialCanvasRendererContext requires SpatialCanvasRendererProvider");
    }
    return value;
}

export function useOptionalSpatialCanvasRendererContext() {
    return useContext(SpatialCanvasRendererContext);
}
