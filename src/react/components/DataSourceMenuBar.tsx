import { type MouseEvent, useCallback, useSyncExternalStore } from "react";
import type { DataSource } from "@/charts/charts";
import AddChartDialogReact from "@/charts/dialogs/AddChartDialogReact";
import { AnnotationDialogReact } from "@/charts/dialogs/AnnotationDialogReact";
import ColorPaletteWrapper from "@/charts/dialogs/ColorPaletteWrapper";
import RenameDatasourceDialog from "@/charts/dialogs/RenameDatasourceDialog";
import type DataStore from "@/datastore/DataStore";
import { getDataSourceDisplayName } from "@/datastore/dataSourceDisplayName";
import { useChartManager } from "../hooks";
import "./DataSourceMenuBar.css";

/** microtip's tooltip attributes; it only shows a tooltip on an element with role="tooltip" */
function tooltipProps(text: string, size: "small" | "medium" = "small", position = "bottom-left") {
    return {
        role: "tooltip",
        "aria-label": text,
        "data-microtip-size": size,
        "data-microtip-position": position,
    } as const;
}

type MenuIconProps = {
    icon: string;
    tooltip: string;
    tooltipPosition?: string;
    size?: string;
    onClick: (event: MouseEvent<HTMLSpanElement>) => void;
};

/**
 * The same markup as `createMenuIcon` in `utilities/Elements.js`, so these icons take the
 * legacy `.ciview-menu-icon` and microtip styles and sit alongside the ones that
 * `ChartManager.addMenuIcon` still appends to the bar.
 */
const MenuIcon = ({ icon, tooltip, tooltipPosition = "bottom-right", size = "18px", onClick }: MenuIconProps) => (
    <span
        role="tooltip"
        aria-label={tooltip}
        data-microtip-size="small"
        data-microtip-position={tooltipPosition}
        onClick={onClick}
    >
        <i className={`ciview-menu-icon ${icon}`} style={{ fontSize: size }} />
    </span>
);

let nextListenerId = 0;

/** The number of rows passing the datastore's filters, kept current as filters change. */
function useFilterSize(dataStore: DataStore) {
    const subscribe = useCallback(
        (onChange: () => void) => {
            const listenerId = `datasource-menu-bar-${nextListenerId++}`;
            dataStore.addListener(listenerId, (type: string) => {
                if (type === "filtered") onChange();
            });
            return () => dataStore.removeListener(listenerId);
        },
        [dataStore],
    );
    return useSyncExternalStore(subscribe, () => dataStore.filterSize);
}

function subscribeToFullscreen(onChange: () => void) {
    document.addEventListener("fullscreenchange", onChange);
    return () => document.removeEventListener("fullscreenchange", onChange);
}

/** Whether `getElement()` is the element currently shown full screen. */
function useIsFullscreen(getElement: () => Element | null | undefined) {
    return useSyncExternalStore(subscribeToFullscreen, () => {
        const element = getElement();
        return !!element && document.fullscreenElement === element;
    });
}

/**
 * The bar above a datasource's charts: add chart, reset filters, colours, layout, tagging and
 * full screen on the left; the datasource's name (with a rename button when the project is
 * editable) and filtered row count on the right.
 *
 * `ChartManager._setUpMenu` mounts this inside `dataSource.menuBar`, each time a view is built.
 */
const DataSourceMenuBar = ({ dataSource }: { dataSource: DataSource }) => {
    const cm = useChartManager();
    const { name, dataStore } = dataSource;
    const displayName = getDataSourceDisplayName(dataStore);
    const isEditable = cm.config.permission === "edit";

    const size = dataStore.size;
    const filterSize = useFilterSize(dataStore);

    // the split pane holding this bar and the datasource's charts
    const getPane = useCallback(() => dataSource.menuBar?.parentElement, [dataSource]);
    const isFullscreen = useIsFullscreen(getPane);

    const handleColorSchemeClick = () => {
        try {
            new ColorPaletteWrapper(cm, dataSource);
        } catch (error) {
            console.error("error making ColorPalette", error);
            cm.createInfoAlert("Error making Color Palette", { type: "warning", duration: 2000 });
        }
    };

    const handleFullscreenClick = async () => {
        try {
            if (isFullscreen) {
                await document.exitFullscreen();
            } else {
                await getPane()?.requestFullscreen();
            }
        } catch (error) {
            console.error("fullscreen error caused: ", error);
        }
    };

    return (
        <>
            <MenuIcon icon="fas fa-chart-bar" tooltip="Add Chart" onClick={() => new AddChartDialogReact(dataStore)} />
            <MenuIcon icon="fas fa-sync-alt" tooltip="Reset All Filters" onClick={() => dataStore.removeAllFilters()} />
            <MenuIcon icon="fas fa-palette" tooltip="Change Color Scheme" onClick={handleColorSchemeClick} />
            <MenuIcon
                icon="fas fa-th"
                tooltip="Change layout"
                onClick={(e) => cm.layoutMenus[name].show(e.nativeEvent)}
            />
            <MenuIcon
                icon="fas fa-tags"
                tooltip="Tag Annotation"
                tooltipPosition="bottom"
                onClick={() => new AnnotationDialogReact(dataStore)}
            />
            <div className="mdv-ds-menu-info">
                <span className={isEditable ? "mdv-ds-name mdv-ds-name--editable" : "mdv-ds-name"}>
                    {/* the tooltip carries the full name, which is truncated when space is short */}
                    <span className="mdv-ds-name-tooltip" {...tooltipProps(displayName, "medium")}>
                        <span className="mdv-ds-name-text">{displayName}</span>
                    </span>
                    {isEditable && (
                        <span {...tooltipProps("Rename datasource")}>
                            <button
                                type="button"
                                className="mdv-ds-rename"
                                aria-label="Rename datasource"
                                onClick={() => new RenameDatasourceDialog(dataStore)}
                            >
                                <i className="fas fa-pen" aria-hidden="true" />
                            </button>
                        </span>
                    )}
                </span>
                <span {...tooltipProps(`${filterSize.toLocaleString()} of ${size.toLocaleString()} rows`)}>
                    <span className="mdv-ds-count">
                        {filterSize.toLocaleString()}
                        <span className="mdv-ds-count-track">
                            <span
                                className="mdv-ds-count-fill"
                                style={{ display: "block", width: `${size ? (filterSize / size) * 100 : 0}%` }}
                            />
                        </span>
                    </span>
                </span>
            </div>
            <MenuIcon
                icon={isFullscreen ? "fas fa-compress" : "fas fa-expand"}
                tooltip={isFullscreen ? "Exit Full Screen" : "Full Screen"}
                // it can sit at the bar's right edge, where a tooltip opening rightwards would widen the page
                tooltipPosition="bottom-left"
                onClick={handleFullscreenClick}
            />
        </>
    );
};

export default DataSourceMenuBar;
