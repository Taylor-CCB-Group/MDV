import { ScatterplotLayer } from "deck.gl/typed";
import { ScatterPlotConfig } from "./components/VivMDVReact";
import { useChart } from "./context";
import { useChartID, useConfig, useFilteredIndices, useParamColumns } from "./hooks";

export function useScatterplotLayer() {
    const id = useChartID();
    const colorBy = (useChart() as any).colorBy;
    const config = useConfig<ScatterPlotConfig>();
    // seem to be reacting fine to changes, why did I think I needed to use extra autorun or reaction?
    const { opacity, radius } = config;

    const data = useFilteredIndices();
    const [cx, cy] = useParamColumns();
    const scatterplotLayer = new ScatterplotLayer({
        id: id + 'scatter-react',
        data,
        opacity,
        radiusScale: radius,
        getFillColor: colorBy ?? [0, 200, 200],
        getRadius: 1,
        getPosition: (i, { target }) => {
            // how slow is a callback function here?
            // given that we need to draw in data from multiple sources...
            // ... unless we want to do the work on a worker thread,
            // I don't think there's a significantly more efficient way
            
            /// === HACK /4 what is going on here? ===
            target[0] = cx.data[i] // 4;
            target[1] = cy.data[i] // 4;
            target[2] = 0;
            return target as unknown as Float32Array; // 🤮 deck.gl types are wrong AFAICT
        },
        updateTriggers: {
            getFillColor: colorBy,
        }
    });
    return scatterplotLayer;
}