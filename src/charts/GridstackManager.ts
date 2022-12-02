import 'gridstack/dist/gridstack.min.css';
import { GridStack } from 'gridstack';
import { debounce } from '../utilities/Utilities';
import { DataSource, Chart, ChartManager } from './charts';
import AnnotationDialog from './dialogs/AnnotationDialog';

function clearPosition(div) {
    div.style.position = '';
    div.style.left = '';
    div.style.right = '';
    div.style.top = '';
    div.style.bottom = '';
    div.style.width = '';
    div.style.height = '';
}

export default class GridStackManager {
    cellHeight = 150;
    cellApproxWidth = 300;
    grids: Map<DataSource, GridStack>;
    dragHandle = '.ciview-chart-title';
    chartManager: ChartManager;
    constructor(chartManager: ChartManager) {
        this.grids = new Map();
        this.chartManager = chartManager;
    }

    getGrid(ds: DataSource) {
        if (!this.grids.has(ds)) {
            const div = ds.contentDiv;
            div.classList.add('grid-stack');
            const rect = div.getBoundingClientRect();
            const column = Math.round(rect.width / this.cellApproxWidth);
            console.log(column);
            const rows = Math.round(rect.height / this.cellHeight);
            this.cellHeight = rect.height / rows; //hack, we might want different sizes for different views, not one 'global'
            console.log('gridstack cellHeight', this.cellHeight);
            const grid = GridStack.init({
                cellHeight: this.cellHeight, handle: this.dragHandle,
                float: true
                // these options not working as expected...
                // margin: 10
                // column
            }, div);
            grid.on('resizestop', (ev, el) => {
                if (el instanceof HTMLElement) {
                    el.style.filter = '';
                }
            });
            grid.on('resizestart', (ev, el) => {
                if (el instanceof HTMLElement) {
                    el.style.filter = 'blur(1px) opacity(0.5)';
                }
            });
            this.chartManager.addMenuIcon(ds.name, "fas fa-th", "compact layout", ()=>grid.compact());
            this.chartManager.addMenuIcon(ds.name, "fas fa-tags", "Tag annotation", () => {new AnnotationDialog(ds.dataStore)})
            this.grids.set(ds, grid);
        }
        return this.grids.get(ds)!;
    }

    manageChart(chart: Chart, ds: DataSource, autoPosition?: boolean) {
        const grid = this.getGrid(ds);
        const div = chart.getDiv();
        const rect = div.getBoundingClientRect();
        const w = Math.round(rect.width / grid.cellWidth());
        const h = Math.round(rect.height / this.cellHeight);
        const x = Math.round(rect.x / grid.cellWidth());
        const y = Math.floor(rect.y / this.cellHeight);
        //maybe better not to add to dom before getting here, and use grid.addWidget
        //current version seems to be working moderately ok.
        clearPosition(div);
        div.style.transition = 'filter 0.5s ease';
        
        // using ResizeObserver rather than gridstack callbacks
        // means we have closure on 'chart' without needing to maintain another data structure
        // or attach more properties to the div.
        const ro = new ResizeObserver(debounce(() => {
            try {
                chart.setSize();
            } catch { }
        }, 20));
        ro.observe(div);
        
        // div.remove();
        // grid.addWidget(div, {w, h, x, y, autoPosition});
        
        grid.makeWidget(div);
        //nb, autoPosition property doesn't apply in update()?
        //passing options to makeWidget() or addWidget() does not evoke joy.
        if (!autoPosition) {
            console.log({w, h, x, y});
            grid.update(div, {w, h, x, y});
        } else {
            grid.update(div, {w, h});
        }
        
        const oldChangeBase = addPopOutHandler(); //revertModifications() used here refers to values stored by subsequent calls
        const lockButton = addPositionLock();
        const oldRemove = addRemoveHandler();
        function addPositionLock() {
            let locked = false;
            const lockButton = chart.addMenuIcon("fas fa-unlock", "lock position");
            const lockIcon = lockButton.children[0];
            lockButton.addEventListener("click", () => {
                locked = !locked;
                lockIcon.classList.toggle("fa-lock", locked);
                lockIcon.classList.toggle("fa-unlock", !locked);
                grid.update(div, {locked});
                div.classList.toggle('gridLock', locked);
            });
            return lockButton;
        }
        function addRemoveHandler() {
            const {remove} = chart;
            chart.remove = () => {
                remove.apply(chart);
                //don't remove from DOM, as it happens elsewhere.
                //doesn't get rid of regl 'must not double destroy framebuffer' error.
                grid.removeWidget(div, false);
            }
            return remove;
        }
        function addPopOutHandler() {
            const {changeBaseDocument} = chart;
            chart.changeBaseDocument = (doc) => {
                changeBaseDocument.apply(chart, [doc]);
                if (doc === document) console.error('changeBaseDocument should be restored to normal by the time we get here...');
                revertModifications();
            }
            return changeBaseDocument;
        }
        // this happens when we close PopOut... we expect to add them all again straight after
        // >>> this implementation might change if we change more of the ChartManager logic for when manageChart is called <<<
        function revertModifications() {
            grid.removeWidget(div, true); //doesn't remove listeners from handle...
            //// leaving listeners as a bug for now.
            lockButton.remove();
            chart.remove = oldRemove;
            chart.changeBaseDocument = oldChangeBase;
        }
    }
}
