import type { SimpleFeatureCollection } from "@deck.gl-community/editable-layers";

export interface Gate {
    id: string;
    name: string;
    geometry: SimpleFeatureCollection;
    columns: [string, string];
    createdAt: number;
    labelPosition: [number, number];
    region?: string;
    color?: [number, number, number];
}
