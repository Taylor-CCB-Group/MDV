import { createMdvPortal } from "@/react/react_utils";
import { createEl } from "../../utilities/ElementsTyped";
import { BaseDialog } from "../../utilities/Dialog.js";
import ColorPaletteComponent from "./ColorPaletteComponent";

interface ColumnColorSyncLink {
    col: string;
    link_to: string;
}

interface SyncColumnColorsConfig {
    columns: ColumnColorSyncLink[];
    dataSource: string;
}

interface LinkColumnsConfig {
    columns: string[];
    dataSource: string;
}

interface ColorDataSource {
    dataChanged: (columns: string[], shouldUpdateFilters?: boolean, shouldCalcExtents?: boolean) => void;
    getColumnColors: (column: string) => string[];
    getColumnList: (datatype: string) => Array<{ field: string; name: string }>;
    getColumnValues: (column: string) => Array<string | number | null | undefined>;
    linkColumns: LinkColumnsConfig[];
    name: string;
    setColumnColors: (column: string, colors: string[]) => void;
    syncColumnColors: SyncColumnColorsConfig[];
}

interface ChartManagerDataSourceItem {
    dataStore: ColorDataSource;
}

interface ChartManagerLike {
    _sync_colors: (
        columns: ColumnColorSyncLink[],
        sourceDataSource: ColorDataSource,
        targetDataSource: ColorDataSource,
    ) => void;
    dataSources: ChartManagerDataSourceItem[];
    getDataSource: (name: string) => ColorDataSource;
}

interface DialogSource {
    dataStore: {
        name: string;
    };
    name: string;
}

interface ColorPaletteInitContent {
    cm: ChartManagerLike;
    ds: DialogSource;
}

class ColorPaletteWrapper extends BaseDialog {
    declare root: ReturnType<typeof createMdvPortal>;
    declare cm: ChartManagerLike;
    declare ds: DialogSource;

    constructor(cm: ChartManagerLike, ds: DialogSource) {
        super(
            {
                width: 500,
                height: 600,
                title: `Color Palette for ${ds.name}`,
            },
            { cm, ds },
        );
    }

    init(content: ColorPaletteInitContent) {
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
        const dataSource = this.cm.getDataSource(this.ds.dataStore.name);
        dataSource.setColumnColors(column, colors);
        dataSource.dataChanged([column], false, false);

        for (const dataSourceItem of this.cm.dataSources) {
            const linkedStore = dataSourceItem.dataStore;
            const syncedColumnColors = linkedStore.syncColumnColors.find((item) => item.dataSource === dataSource.name);

            if (syncedColumnColors) {
                const linkedColumn = syncedColumnColors.columns.find((item) => item.link_to === column);
                if (linkedColumn) {
                    this.cm._sync_colors([linkedColumn], dataSource, linkedStore);
                    linkedStore.dataChanged([linkedColumn.col], false, false);
                }
            }

            const linkedColumns = linkedStore.linkColumns.find((item) => item.dataSource === dataSource.name);
            if (linkedColumns?.columns.includes(column)) {
                linkedStore.setColumnColors(column, colors);
                linkedStore.dataChanged([column], false, false);
            }
        }
    }

    close(): void {
        this.root?.unmount();
        super.close();
    }
}

export default ColorPaletteWrapper;
export type { ChartManagerLike, ColorDataSource, DialogSource };
