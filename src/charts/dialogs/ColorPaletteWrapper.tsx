import { createMdvPortal } from "@/react/react_utils";
import type { ChartManager } from "../ChartManager.js";
import type DataStore from "../../datastore/DataStore.js";
import { BaseDialog } from "../../utilities/Dialog.js";
import ColorPaletteComponent from "./ColorPaletteComponent";

type ColorPaletteChartManager = Pick<ChartManager, "_sync_colors" | "dataSources" | "getDataSource">;

class ColorPaletteWrapper extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;
    private cm: ColorPaletteChartManager;
    private dataSource: DataStore;

    constructor(cm: ColorPaletteChartManager, ds: ChartManager["dataSources"][number]) {
        super(
            {
                width: 500,
                height: 600,
                title: `Color Palette for ${ds.name}`,
            },
            null,
        );
        this.cm = cm;
        this.dataSource = cm.getDataSource(ds.dataStore.name);
        this.dialog.style.padding = "0";
        this.dialog.style.display = "flex";
        this.dialog.style.flexDirection = "column";
        this.dialog.style.overflow = "hidden";
        const container = document.createElement("div");
        Object.assign(container.style, {
            width: "100%",
            height: "100%",
            display: "flex",
            flex: "1 1 auto",
            minHeight: "0",
            padding: "0",
        });
        this.dialog.append(container);
        this.root = createMdvPortal(
            <ColorPaletteComponent
                dataSource={this.dataSource}
                dataSourceLabel={ds.name}
                onApply={(column, colors) => this.applyColors(column, colors)}
            />,
            container,
            this,
        );
    }

    applyColors(column: string, colors: string[]) {
        const { cm, dataSource } = this;
        dataSource.setColumnColors(column, colors);
        dataSource.dataChanged([column], false);

        for (const dataSourceItem of cm.dataSources) {
            const linkedStore: DataStore = dataSourceItem.dataStore;
            const syncedColumns = linkedStore.syncColumnColors
                .filter((item: { dataSource: string }) => item.dataSource === dataSource.name)
                .flatMap((item: { columns: { col: string; link_to: string }[] }) =>
                    item.columns.filter((linkedColumn) => linkedColumn.link_to === column),
                );

            if (syncedColumns.length > 0) {
                cm._sync_colors(syncedColumns, dataSource, linkedStore);
                linkedStore.dataChanged(
                    syncedColumns.map((linkedColumn: { col: string }) => linkedColumn.col),
                    false,
                );
            }

            const hasLinkedColumn = linkedStore.linkColumns.some(
                (item: { dataSource: string; columns: string[] }) =>
                    item.dataSource === dataSource.name && item.columns.includes(column),
            );
            if (hasLinkedColumn) {
                linkedStore.setColumnColors(column, colors);
                linkedStore.dataChanged([column], false);
            }
        }
    }

    close(): void {
        super.close();
        this.root.unmount();
    }
}

export default ColorPaletteWrapper;
