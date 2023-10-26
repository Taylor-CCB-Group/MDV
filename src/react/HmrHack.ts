import vSideEffect from "../react/components/VivMDVReact";
const arr = [];
arr.push(vSideEffect);
import SideEffect from "../charts/dialogs/AnnotationDialogReact";
arr.push(SideEffect);

/**
 * Charts are registered with ChartManager by being added to a BaseChart.types object,
 * generally as a side effect of importing the chart module. For some reason, I've found
 * that having some kind of explicit export from the module helps HMR to work in the MDV context.
 * 
 * As of this writing... I haven't completely understood what combination of things is 
 * needed to make sure that HMR works properly for these modules; VivMDVReact HMR was working,
 * then I made some changes to get a new ColorChannelDialogReactWrapper used by it working... at some
 * point it regressed. Gradually reintroducing the changes I'd made in development to figure out what caused
 * the problem, WIP...
 * 
 * We may need to explicitly use the HMR API to make sure that the modules are reloaded in
 * correctly. I've used it in the past, and was pleasantly surprised that it worked without
 * too much hassle, but in principle would rather avoid it if possible.
 */
export default function DoImportSideEffects() {
    console.log(`DoImportSideEffects: ${arr.length} modules imported`);
}