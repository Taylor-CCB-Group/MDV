import 'gridstack/dist/gridstack.min.css';
import { GridStack } from 'gridstack';

export default class GridStackManager {
    constructor(chartManager) {
        this.grids = new Map();
        this.chartManager = chartManager;
    }

    getGrid(ds) {
        if (!this.grids.has(ds)) {
            const div = ds.contentDiv;
            div.classList.add('grid-stack');
            const grid = GridStack.init({
                cellHeight: 300, handle: '.ciview-chart-title'
            }, div);
            grid.on('resizestop', (ev, el) => {
                const {width, height} = el.gridstackNode;
                el.__chart.setSize(width, height);
            });
            this.grids.set(ds, grid);
        }
        return this.grids.get(ds);
    }

    manageChart(chart, ds) {
        const grid = this.getGrid(ds);
        const div = chart.getDiv();
        // div.remove();
        div.__chart = chart; // not really my style
        grid.addWidget(chart.getDiv(), {autoPosition: true});
        grid.makeWidget(div);
    }
}
