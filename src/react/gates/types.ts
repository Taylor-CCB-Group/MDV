import type { FeatureCollection } from '@turf/helpers';

export interface Gate {
    id: string;
    name: string;
    geometry: FeatureCollection;
    columns: [string, string];
    createdAt: number;
}

