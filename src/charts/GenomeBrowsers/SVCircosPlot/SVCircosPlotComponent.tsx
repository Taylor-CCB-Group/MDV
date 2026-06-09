import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import type DataStore from "@/datastore/DataStore";
import { rebindMouseEvents } from "@/lib/deckMonkeypatch";
import { useChartSize } from "@/react/hooks";
import { useBinaryChartColors } from "@/react/hooks/useBinaryChartColors";
import { OrthographicView, type OrthographicViewState, type ViewStateChangeParameters } from "@deck.gl/core";
import { PathLayer, TextLayer } from "@deck.gl/layers";
import DeckGL, { type DeckGLRef } from "@deck.gl/react";
import { observer } from "mobx-react-lite";
import { useEffect, useMemo, useRef, useState } from "react";
import { v4 as uuidv4 } from "uuid";
import { useChart } from "../../../react/context";
import { useOuterContainer } from "../../../react/screen_state";
import { useGenomicInfo } from "../genomicLocationUtils";
import SVLayer from "./SVLayer.js";

type RgbColor = [number, number, number];
type RgbaColor = [number, number, number, number];
type Point = [number, number];

interface FilterState {
    timeStamp?: number;
    dimension?: unknown;
}

interface ChromosomeSegment {
    chr: string;
    start: number;
    length: number;
    index: number;
}

interface LabelDatum {
    text: string;
    position: Point;
}

interface TickLineDatum {
    path: Point[];
    major: boolean;
}

interface SVBinaryData {
    attributes: {
        getSV: Float32Array;
        getColor: Uint8Array;
        getFilter: Uint8Array;
    };
}

interface CalculatedSVs {
    data: SVBinaryData;
    totalLength: number;
    chromosomeSegments: ChromosomeSegment[];
}

interface HoverInfo {
    index?: number;
}

type SVCircosChartConfig = BaseConfig & {
    color_by?: string;
    color_legend?: { display: boolean };
    param?: string[];
};

type SVCircosChart = BaseChart<SVCircosChartConfig> & {
    dataStore: DataStore;
    setColorLegend: () => void;
};

type GenomeDataStore = DataStore & {
    filterArray: Uint8Array;
    genome: {
        chromosomes: Record<string, number>;
    };
    size: number;
};

const CHROMOSOME_COLORS: RgbColor[] = [
    [31, 119, 180],
    [255, 127, 14],
    [44, 160, 44],
    [214, 39, 40],
    [148, 103, 189],
    [140, 86, 75],
    [227, 119, 194],
    [127, 127, 127],
    [188, 189, 34],
    [23, 190, 207],
];

// don't want to use filteredindexes as the filter buffer is used directly
// in webgl, just need to know when filtering happens
const useFilter = (datastore: DataStore): FilterState => {
    const [filter, setFilter] = useState<FilterState>({});

    useEffect(() => {
        const id = uuidv4();
        datastore.addListener(id, (type: string, data: unknown) => {
            if (type === "filtered") {
                setFilter({
                    timeStamp: Date.now(),
                    dimension: data,
                });
            }
        });
        return () => {
            datastore.removeListener(id);
        };
    }, [datastore]);

    return filter;
};

// convert svs into a flat array of start, length, type
// for webgl - should move this to a web worker if it gets too slow
function calculateSVs(svColumns: any, dataStore: GenomeDataStore, colorArray: Uint8Array): CalculatedSVs {
    const chrs = Object.keys(dataStore.genome.chromosomes).sort((a, b) =>
        a.localeCompare(b, undefined, { numeric: true, sensitivity: "base" }),
    );
    const chrCumLen: Record<string, number> = {};
    let cumSize = 0;
    for (const chr of chrs) {
        chrCumLen[chr] = cumSize;
        cumSize += dataStore.genome.chromosomes[chr];
    }

    const svMap: Record<string, number> = {
        TRA: 0,
        BND: 0,
        INV: 1,
        DEL: 2,
        INS: 3,
        DUP: 4,
    };

    const svBuff = new Float32Array(dataStore.size * 3);

    for (let i = 0; i < dataStore.size; i++) {
        const idx = i * 3;
        const chr1Value = String(svColumns.chr1.getValue(i) ?? "");
        const chr2Value = String(svColumns.chr2.getValue(i) ?? "");
        const pos1Value = Number(svColumns.pos1.getValue(i) ?? 0);
        const pos2Value = Number(svColumns.pos2.getValue(i) ?? 0);
        const svStart = pos1Value + (chrCumLen[chr1Value] ?? 0);
        const svEnd = pos2Value + (chrCumLen[chr2Value] ?? 0);
        const svType = String(svColumns.svtype.getValue(i) ?? "")
            .trim()
            .toUpperCase();
        const typeIndex = svMap[svType] ?? 5;
        svBuff[idx] = svStart;
        // Keep signed delta so links where end is before start are rendered correctly.
        svBuff[idx + 1] = svEnd - svStart;
        svBuff[idx + 2] = typeIndex;
    }

    const chromosomeSegments = chrs.map((chr, index) => ({
        chr,
        start: chrCumLen[chr],
        length: dataStore.genome.chromosomes[chr],
        index,
    }));

    return {
        data: {
            attributes: {
                getSV: svBuff,
                getColor: colorArray,
                getFilter: dataStore.filterArray,
            },
        },
        totalLength: cumSize,
        chromosomeSegments,
    };
}

function getChromosomeSegments(dataStore: GenomeDataStore): {
    chromosomeSegments: ChromosomeSegment[];
    totalLength: number;
} {
    const chrs = Object.keys(dataStore.genome.chromosomes).sort((a, b) =>
        a.localeCompare(b, undefined, { numeric: true, sensitivity: "base" }),
    );
    const chrCumLen: Record<string, number> = {};
    let cumSize = 0;
    for (const chr of chrs) {
        chrCumLen[chr] = cumSize;
        cumSize += dataStore.genome.chromosomes[chr];
    }
    return {
        chromosomeSegments: chrs.map((chr, index) => ({
            chr,
            start: chrCumLen[chr],
            length: dataStore.genome.chromosomes[chr],
            index,
        })),
        totalLength: cumSize,
    };
}

function getEmptyCalculatedSVs(dataStore: GenomeDataStore): CalculatedSVs {
    const { chromosomeSegments, totalLength } = getChromosomeSegments(dataStore);
    return {
        data: {
            attributes: {
                getSV: new Float32Array(0),
                getColor: new Uint8Array(0),
                getFilter: new Uint8Array(0),
            },
        },
        totalLength,
        chromosomeSegments,
    };
}

function getChromosomePath(start: number, length: number, totalLength: number, radius: number, steps = 48): Point[] {
    const startAngle = (start * 2 * Math.PI) / totalLength - Math.PI / 2;
    const endAngle = ((start + length) * 2 * Math.PI) / totalLength - Math.PI / 2;
    const stepCount = Math.max(2, steps);
    const points: Point[] = [];
    for (let i = 0; i <= stepCount; i++) {
        const t = i / stepCount;
        const angle = startAngle + (endAngle - startAngle) * t;
        points.push([Math.cos(angle) * radius, Math.sin(angle) * radius]);
    }
    return points;
}

function pointAtGenomePosition(position: number, totalLength: number, radius: number): Point {
    const angle = (position * 2 * Math.PI) / totalLength - Math.PI / 2;
    return [Math.cos(angle) * radius, Math.sin(angle) * radius];
}

function getChromosomeLabelData(
    chromosomeSegments: ChromosomeSegment[],
    totalLength: number,
    radius: number,
): LabelDatum[] {
    return chromosomeSegments.map((segment) => {
        const mid = segment.start + segment.length / 2;
        return {
            text: segment.chr,
            position: pointAtGenomePosition(mid, totalLength, radius),
        };
    });
}

function formatTickLabel(bp: number): string {
    if (bp >= 1000000) return `${Math.round(bp / 1000000)}Mb`;
    if (bp >= 1000) return `${Math.round(bp / 1000)}kb`;
    return `${bp}bp`;
}

function getTickData(
    chromosomeSegments: ChromosomeSegment[],
    totalLength: number,
    tickStep = 25000000,
    majorEvery = 2,
): { outerTickLines: TickLineDatum[]; outerTickLabels: LabelDatum[] } {
    const outerTickLines: TickLineDatum[] = [];
    const outerTickLabels: LabelDatum[] = [];

    chromosomeSegments.forEach((segment) => {
        const ticks = Math.floor(segment.length / tickStep);
        for (let i = 0; i <= ticks; i++) {
            const offset = i * tickStep;
            if (offset > segment.length) break;
            const genomePos = segment.start + offset;
            const pOuter1 = pointAtGenomePosition(genomePos, totalLength, 1052);
            const isMajor = i % majorEvery === 0;
            const pOuter2 = pointAtGenomePosition(genomePos, totalLength, isMajor ? 1080 : 1066);
            outerTickLines.push({
                path: [pOuter1, pOuter2],
                major: isMajor,
            });

            if (isMajor && i > 0) {
                const outerLabelPos = pointAtGenomePosition(genomePos, totalLength, 1095);
                outerTickLabels.push({
                    text: formatTickLabel(offset),
                    position: outerLabelPos,
                });
            }
        }
    });

    return { outerTickLines, outerTickLabels };
}

function getDisplaySpacing(zoomLevel: number): { tickStep: number; majorEvery: number } {
    if (zoomLevel < -0.4) return { tickStep: 100000000, majorEvery: 1 };
    if (zoomLevel < 0.3) return { tickStep: 50000000, majorEvery: 1 };
    if (zoomLevel < 1.0) return { tickStep: 25000000, majorEvery: 2 };
    return { tickStep: 10000000, majorEvery: 2 };
}

const SVCircosPlotComponent = observer(() => {
    const chart = useChart() as SVCircosChart;
    const { areColumnsLoaded, genomicColumns } = useGenomicInfo();

    const [colorArray, colorChange] = useBinaryChartColors(chart);

    const { data, totalLength, chromosomeSegments } = useMemo(
        () =>
            areColumnsLoaded
                ? calculateSVs(genomicColumns, chart.dataStore as GenomeDataStore, colorArray)
                : getEmptyCalculatedSVs(chart.dataStore as GenomeDataStore),
        [chart.dataStore, genomicColumns],
    );

    const outerContainer = useOuterContainer();
    const deckRef = useRef<DeckGLRef<OrthographicView | null> | null>(null);
    const [width, height] = useChartSize();
    const [viewState, setViewState] = useState<OrthographicViewState>({
        target: [0, 0],
        zoom: [0, 0],
    });
    const [hoveredIndex, setHoveredIndex] = useState(-1);
    const filter = useFilter(chart.dataStore);

    useEffect(() => {
        deckRef.current?.deck?.redraw();
    }, [colorArray]);

    useEffect(() => {
        const plotDiameter = 2 * 1120 + 100;
        const minDim = Math.min(width, height);
        const zoom = Math.log2(minDim / plotDiameter);
        setViewState({
            target: [0, 0],
            zoom: [zoom, zoom],
        });
    }, [width, height]);

    const handleHover = (info: HoverInfo) => {
        setHoveredIndex(info.index !== undefined && info.index >= 0 ? info.index : -1);
    };

    const handleClick = (info: HoverInfo) => {
        if (info.index !== undefined && info.index >= 0) {
            chart.dataStore.dataHighlighted([info.index], chart);
            setHoveredIndex(info.index);
        }
    };

    useEffect(() => {
        void outerContainer;
        if (deckRef.current?.deck) {
            try {
                return rebindMouseEvents(deckRef.current.deck);
            } catch (e) {
                console.error(
                    "attempt to reset deck eventManager element failed - could be related to brittle deck monkeypatch",
                    e,
                );
            }
        }
    }, [outerContainer]);

    const zoomLevel = Array.isArray(viewState.zoom) ? viewState.zoom[0] : (viewState.zoom ?? 0);
    const { tickStep, majorEvery } = getDisplaySpacing(zoomLevel);
    const showTicks = zoomLevel > -0.9;
    const showTickLabels = showTicks;
    const view = useMemo(() => new OrthographicView(), []);

    const chromosomeLayer = new PathLayer<ChromosomeSegment>({
        id: "chromosome_layer",
        data: chromosomeSegments,
        getPath: (d) => getChromosomePath(d.start, d.length, totalLength, 1045),
        getColor: (d) => CHROMOSOME_COLORS[d.index % CHROMOSOME_COLORS.length],
        getWidth: 8,
        widthUnits: "pixels",
        widthMinPixels: 2,
        pickable: false,
        capRounded: false,
        jointRounded: false,
    });

    const chromosomeLabels = getChromosomeLabelData(chromosomeSegments, totalLength, 1120);
    const { outerTickLines, outerTickLabels } = getTickData(chromosomeSegments, totalLength, tickStep, majorEvery);

    const chromosomeLabelLayer = new TextLayer<LabelDatum>({
        id: "chromosome_label_layer",
        data: chromosomeLabels,
        getPosition: (d) => d.position,
        getText: (d) => d.text,
        getColor: [30, 30, 30, 255] satisfies RgbaColor,
        getSize: 14,
        sizeUnits: "pixels",
        sizeMinPixels: 10,
        getTextAnchor: "middle",
        getAlignmentBaseline: "center",
        pickable: false,
    });

    const tickLayer = new PathLayer<TickLineDatum>({
        id: "chromosome_tick_layer",
        data: outerTickLines,
        getPath: (d) => d.path,
        getColor: (d) => (d.major ? [70, 70, 70, 255] : [120, 120, 120, 220]),
        getWidth: (d) => (d.major ? 2 : 1),
        widthUnits: "pixels",
        widthMinPixels: 1,
        visible: showTicks,
        pickable: false,
    });

    const outerTickLabelLayer = new TextLayer<LabelDatum>({
        id: "chromosome_tick_outer_label_layer",
        data: outerTickLabels,
        getPosition: (d) => d.position,
        getText: (d) => d.text,
        getColor: [80, 80, 80, 255] satisfies RgbaColor,
        getSize: 10,
        sizeUnits: "pixels",
        sizeMinPixels: 8,
        getTextAnchor: "middle",
        getAlignmentBaseline: "center",
        visible: showTickLabels,
        pickable: false,
    });

    const layers = [
        chromosomeLayer,
        tickLayer,
        outerTickLabelLayer,
        chromosomeLabelLayer,
        new SVLayer({
            id: "sv_layer",
            data,
            pickable: true,
            highlightColor: [0, 0, 0, 255],
            highlightedObjectIndex: hoveredIndex,
            radius: 1000,
            totalLength,
                container: outerContainer,
                center: [0, 0],
                rotation: 0,
                updateTriggers: {
                    getFilter: [filter.timeStamp],
                    getColor: [colorChange],
                },
            }),
    ];

    const handleViewStateChange = ({ viewState: nextViewState }: ViewStateChangeParameters<OrthographicViewState>) => {
        setViewState(nextViewState);
    };

    return !areColumnsLoaded ? (
        <div style={{ padding: "8px", fontSize: "12px" }}>Loading genomic columns...</div>
    ) : (
        <DeckGL
            ref={deckRef}
            layers={layers}
            viewState={viewState}
            onViewStateChange={handleViewStateChange}
            onClick={handleClick}
            onHover={handleHover}
            getCursor={() => (hoveredIndex >= 0 ? "pointer" : "default")}
            views={view}
            controller={true}
        />
    );
});

export default SVCircosPlotComponent;
