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

export default class GridStackManager {
    cellHeight = 300;
    cellApproxWidth = 300;
    constructor(chartManager) {
        this.grids = new Map();
        this.chartManager = chartManager;
    }

    getGrid(ds) {
        if (!this.grids.has(ds)) {
            const div = ds.contentDiv;
            div.classList.add('grid-stack');
            const column = Math.round(div.getBoundingClientRect().width / this.cellApproxWidth);
            console.log(column);
            const grid = GridStack.init({
                cellHeight: this.cellHeight, handle: '.ciview-chart-title',
                // these options not working as expected...
                // margin: 10
                // column
            }, div);
            grid.on('resizestop', (ev, el) => {
                // console.log('resizestop', ev, el);
                el.style.filter = '';
            });
            grid.on('resizestart', (ev, el) => {
                el.style.filter = 'blur(1px) opacity(0.5)';
            });
            this.grids.set(ds, grid);
        }
        return this.grids.get(ds);
    }

    manageChart(chart, ds) {
        const grid = this.getGrid(ds);
        const div = chart.getDiv();
        const rect = div.getBoundingClientRect();
        const w = Math.round(rect.width / grid.cellWidth());
        const h = Math.round(rect.height / this.cellHeight);
        const x = Math.floor(rect.x / grid.cellWidth());
        const y = Math.floor(rect.y / this.cellHeight);
        // div.remove();
        //maybe better not to add to dom before getting here, and use grid.addWidget
        //current version seems to be working moderately ok.
        clearPosition(div);
        // grid.addWidget(chart.getDiv(), {autoPosition: true, w, h, x, y});
        div.style.transition = 'filter 0.5s ease';
        grid.makeWidget(div);
        //adding x, y seems to be making infinite loop
        //resize doesn't work well until items are properly integrated in grid.
        grid.update(div, {
            w, h//, x, y
        });
        // using ResizeObserver rather than gridstack callbacks
        // means we have closure on 'chart' without needing to maintain another data structure
        // or attach more properties to the div.
        const ro = new ResizeObserver(debounce(() => {
            try {
                chart.setSize();
            } catch {}
        }, 20));
        ro.observe(div);
        // chart.setSize();
    }
}
