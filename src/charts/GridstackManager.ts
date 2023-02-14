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

declare global {
    interface HTMLDivElement {
        gridstackNode?: any;
        gridstackPopoutCallback?: () => void;
    }
}
type GridInstance = {
    grid:GridStack;
    charts:Set<Chart>;
    icon:HTMLElement;
}

export default class GridStackManager {
    cellHeight = 150;
    cellApproxWidth = 300;
    grids: Map<DataSource,GridInstance>;
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
                cellHeight: this.cellHeight,
                handle: this.dragHandle,
                float: true,
                oneColumnSize:400
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
         
            const i = this.chartManager.addMenuIcon(ds.name, "fas fa-compress-arrows-alt", "compact layout", ()=>grid.compact());
            //this.chartManager.addMenuIcon(ds.name, "fas fa-tags", "Tag annotation", () => {new AnnotationDialog(ds.dataStore)})
            this.grids.set(ds, {grid:grid,charts:new Set(),icon:i});
        }
        
        return this.grids.get(ds)!;
    }

    destroy(ds: DataSource){
        if (!this.grids.has(ds)){
            console.error(`destroying non-existent gridstack - ${ds.name}`)
        }
        const gi = this.grids.get(ds);
        const div = ds.contentDiv;
         //store sizes/positions of div elements
        const sizes = new Map();
        for (let chart of gi.charts){
            const d = chart.getDiv();
            sizes.set(d,[d.offsetWidth,d.offsetHeight,d.offsetLeft,d.offsetTop]);
            chart.removeLayout();
            delete chart.config.gssize;
            delete chart.config.gsposition;
        }
                 
        gi.grid.destroy(false);

        //convert back to absolute positioning plus other clean up on the div
        for (let chart of gi.charts){
            const  d= chart.getDiv();
            const s = sizes.get(d);
            d.style.position="absolute";
            d.style.width=(s[0]-5)+"px";
            d.style.height=(s[1]-5)+"px";
            d.style.left=s[2]+"px";
            d.style.top=s[3]+"px";
            chart.config.size=[s[0]-5,s[1]-5];
            chart.config.position=[s[2],s[3]];
            d.classList.remove("grid-stack-item")
        }
        div.classList.remove('grid-stack');
        gi.icon.remove();
        this.grids.delete(ds)
    }


    getCellDimensions(ds:DataSource){
        const gridInstance = this.getGrid(ds);
        const grid= gridInstance.grid;
        return [grid.cellWidth(),this.cellHeight]
    }


    manageChart(chart: Chart, ds: DataSource, autoPosition?: boolean) {
        const gridInstance = this.getGrid(ds);
        const grid= gridInstance.grid;
        gridInstance.charts.add(chart);
        const div = chart.getDiv();
        const rect = div.getBoundingClientRect();
        let w = Math.round(rect.width / grid.cellWidth());
        let h = Math.round(rect.height / this.cellHeight);
        let x = Math.round(div.offsetLeft / grid.cellWidth()); //relative location
        //cases where the old container was wider than the current one
        x=x>11?11:x;
        let y = Math.floor(div.offsetTop / this.cellHeight);  //relative location
        if (chart.config.gsposition){
            [x,y]= chart.config.gsposition;
            [w,h]= chart.config.gssize;
            autoPosition=false;
        }
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
                grid.removeWidget(div, true);
                gridInstance.charts.delete(chart);

                //div.gridstackNode = undefined; //doesn't help
                remove.apply(chart);
            }
            return remove;
        }
        function addPopOutHandler() {
            const {changeBaseDocument} = chart;
            div.gridstackPopoutCallback = () => {
                // console.log('removing from gridstack');
                grid.removeWidget(div, true);
            }
            chart.changeBaseDocument = (doc) => {
                //by now, it's too late... the parent element is no longer the grid.
                //grid.removeWidget(div, true);
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
        chart.removeLayout=()=>{
            grid.removeWidget(div,false); //doesn't remove listeners from handle...
            lockButton.remove();
            delete div.gridstackPopoutCallback;
            chart.remove = oldRemove;
            chart.changeBaseDocument = oldChangeBase;
            ro.disconnect();
        }

        
    }
    
}
