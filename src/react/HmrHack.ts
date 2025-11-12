import vSideEffect from "../react/components/VivMDVReact";
const arr: any[] = [];
arr.push(vSideEffect);
import dSideEffect from "../react/components/DeckScatterReactWrapper";
arr.push(dSideEffect);
import SideEffect from "../charts/dialogs/AnnotationDialogReact";
arr.push(SideEffect);
// import addChart from "../charts/dialogs/AddChartDialogReact";
// arr.push(addChart);
import hfc from "./components/HighlightedFeatureChart";
arr.push(hfc);
import "./components/SelectionDialogReact";
import "../charts/dialogs/AddColumnsFromRowsDialogReact";
import tSideEffect from "../react/components/TableChartReactWrapper";
arr.push(tSideEffect);

/**
 * Charts are registered with ChartManager by being added to a BaseChart.types object,
 * generally as a side effect of importing the chart module. For some reason, I've found
 * that having some kind of explicit export from the module helps HMR to work in the MDV context.
 *
 * Something should be done with the result of the import, otherwise it gets tree-shaken away.
 *
 * I've also realised that you *must not* export the chart class from the module.
 *
 * Better documentation of this would be good.
 *
 * We may need to explicitly use the HMR API to make sure that the modules are reloaded in
 * correctly. I've used it in the past, and was pleasantly surprised that it worked without
 * too much hassle, but in principle would rather avoid it if possible.
 */
export default function DoImportSideEffects() {
    console.log(`DoImportSideEffects: ${arr.length} modules imported`);
}
