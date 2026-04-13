import { createMdvPortal } from "@/react/react_utils";
import type { ChartManager } from "../ChartManager.js";
import type DataStore from "../../datastore/DataStore.js";
import { createEl } from "../../utilities/ElementsTyped";
import { BaseDialog } from "../../utilities/Dialog.js";
import ColorPaletteComponent from "./ColorPaletteComponent";

type ColorPaletteChartManager = Pick<ChartManager, "_sync_colors" | "dataSources" | "getDataSource">;
type ChartManagerDataSource = ChartManager["dataSources"][number];
type ColorPaletteContent = { cm: ColorPaletteChartManager; ds: ChartManagerDataSource };
type SyncColumnColorLink = { col: string; link_to: string };
type SyncColumnColorMapping = { columns: SyncColumnColorLink[]; dataSource: string };
type LinkedColumnsMapping = { columns: string[]; dataSource: string };

class ColorPaletteWrapper extends BaseDialog {
    declare root: ReturnType<typeof createMdvPortal>;
    declare cm: ColorPaletteChartManager;
    declare ds: ChartManagerDataSource;

    constructor(cm: ColorPaletteChartManager, ds: ChartManagerDataSource) {
        super(
            {
                width: 500,
                height: 600,
                title: `Color Palette for ${ds.name}`,
            },
            { cm, ds },
        );
    }

    init(content: ColorPaletteContent) {
        this.cm = content.cm;
        this.ds = content.ds;

        this.dialog.style.padding = "0";
        this.dialog.style.display = "flex";
        this.dialog.style.flexDirection = "column";
        this.dialog.style.overflow = "hidden";

        const container = createEl(
            "div",
            {
                styles: {
                    width: "100%",
                    height: "100%",
                    display: "flex",
                    flex: "1 1 auto",
                    minHeight: "0",
                    padding: "0",
                },
            },
            this.dialog,
        );

        this.root = createMdvPortal(
            <ColorPaletteComponent
                dataSource={this.cm.getDataSource(this.ds.dataStore.name)}
                dataSourceLabel={this.ds.name}
                onApply={(column, colors) => {
                    this.applyColors(column, colors);
                }}
            />,
            container,
            this,
        );
    }

    applyColors(column: string, colors: string[]) {
        const dataSource: DataStore = this.cm.getDataSource(this.ds.dataStore.name);
        dataSource.setColumnColors(column, colors);
        dataSource.dataChanged([column], false);

        for (const dataSourceItem of this.cm.dataSources) {
            const linkedStore: DataStore = dataSourceItem.dataStore;
            const syncedColumns = linkedStore.syncColumnColors
                .filter((item: SyncColumnColorMapping) => item.dataSource === dataSource.name)
                .flatMap((item: SyncColumnColorMapping) =>
                    item.columns.filter((linkedColumn: SyncColumnColorLink) => linkedColumn.link_to === column),
                );

            if (syncedColumns.length > 0) {
                this.cm._sync_colors(syncedColumns, dataSource, linkedStore);
                linkedStore.dataChanged(syncedColumns.map((linkedColumn) => linkedColumn.col), false);
            }

            const hasLinkedColumn = linkedStore.linkColumns.some(
                (item: LinkedColumnsMapping) => item.dataSource === dataSource.name && item.columns.includes(column),
            );
            if (hasLinkedColumn) {
                linkedStore.setColumnColors(column, colors);
                linkedStore.dataChanged([column], false);
            }
        }
    }

    close(): void {
        this.root?.unmount();
        super.close();
    }
}

export default ColorPaletteWrapper;
