import 'gridstack/dist/gridstack.min.css';
import { GridStack } from 'gridstack';
import { debounce } from '../utilities/Utilities';

function clearPosition(div) {
    div.style.position = '';
    div.style.left = '';
    div.style.right = '';
    div.style.top = '';
    div.style.bottom = '';
    div.style.width = '';
    div.style.height = '';
}
type DataSource = {contentDiv: HTMLDivElement};
type Chart = {
    getDiv: ()=>HTMLDivElement;
    remove: ()=>void;
    addMenuIcon: (classes: string, info: string) => HTMLElement;
    setSize: (x?: number, y?: number) => void;
};

export default class GridStackManager {
    cellHeight = 150;
    cellApproxWidth = 300;
    grids: Map<DataSource, GridStack>;
    constructor() {
        this.grids = new Map();
    }

    getGrid(ds: DataSource) {
        //ds is undefined for PopOuts... and lo, they behave pretty poorly ATM.
        if (!this.grids.has(ds)) {
            const div = ds.contentDiv;
            div.classList.add('grid-stack');
            const rect = div.getBoundingClientRect();
            const column = Math.round(rect.width / this.cellApproxWidth);
            console.log(column);
            const rows = Math.round(rect.height / this.cellHeight);
            this.cellHeight = rect.height / rows; //hack, we could have different sizes for different ds etc.
            console.log('gridstack cellHeight', this.cellHeight);
            const grid = GridStack.init({
                cellHeight: this.cellHeight, handle: '.ciview-chart-title',
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
            this.grids.set(ds, grid);
        }
        return this.grids.get(ds);
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
        if (autoPosition) {
            grid.update(div, {w, h, x, y});
        } else {
            grid.update(div, {w, h});
        }
        
        addPositionLock();
        addRemoveHandler();
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
        }
        function addRemoveHandler() {
            const {remove} = chart;
            chart.remove = () => {
                remove.apply(chart);
                //don't remove from DOM, as it happens elsewhere.
                //doesn't get rid of regl 'must not double destroy framebuffer' error.
                grid.removeWidget(div, false);
            }
        }
    }
}
