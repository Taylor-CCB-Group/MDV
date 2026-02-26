import "../react/components/VivMDVReact";
import "../react/components/DeckScatterReactWrapper";
import "../charts/dialogs/AnnotationDialogReact";
// import addChart from "../charts/dialogs/AddChartDialogReact";
// arr.push(addChart);
import "./components/HighlightedFeatureChart";
import "./components/GeneNetworkChart";
import "./components/SelectionDialogReact";
import "../charts/dialogs/AddColumnsFromRowsDialogReact";
import "../react/components/TableChartReactWrapper";

/**
 * Charts are registered with ChartManager by being added to a BaseChart.types object,
 * generally as a side effect of importing the chart module. For some reason, in the past I'd found
 * that having some kind of explicit export from the module helps HMR to work in the MDV context.
 * 
 * This no longer appears to be the case - as of December 2025, it seems that importing the modules
 * and having the BaseChart.types mutation side-effect is enough. Still generally keeping react HMR-able
 * code separate from chart classes.
 * 
 * So - the imports done in this file aren't so much of an HmrHack any more.
 *
 */
