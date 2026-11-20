import { DrawPolygonByDraggingMode, type GeoJsonEditMode, TranslateMode } from "@deck.gl-community/editable-layers";
import { GridViewOutlined, TuneOutlined } from "@mui/icons-material";
import AddCircleOutlinedIcon from "@mui/icons-material/AddCircleOutlined";
import ChevronLeftIcon from "@mui/icons-material/ChevronLeft";
import ChevronRightIcon from "@mui/icons-material/ChevronRight";
import CloseOutlinedIcon from "@mui/icons-material/CloseOutlined";
import DoneOutlinedIcon from "@mui/icons-material/DoneOutlined";
import LayersOutlinedIcon from "@mui/icons-material/LayersOutlined";
import PanToolOutlinedIcon from "@mui/icons-material/PanToolOutlined";
import { ButtonGroup, Divider } from "@mui/material"; //check tree-shaking...
import { LassoIcon, SplinePointerIcon, VectorSquareIcon } from "lucide-react";
import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { EditMode, PolygonMode, RectangleMode } from "@/editable-layers/selectionModes";
import { useChart } from "../context";
import type { DualContourLegacyConfig } from "../contour_state";
import { useChartID, useConfig, useParamColumns, useTheme } from "../hooks";
import useGateActions from "../hooks/useGateActions";
import { useSpatialLayers } from "../spatial_context";
import type { DeckScatterConfigWithRegion } from "./DeckScatterReactWrapper";
import { densityGridConfig, isDensityGridShown, supportsDensityGridMode } from "./densityGridUtils";
import GateNameDialog from "./GateNameDialog";
import IconWithTooltip from "./IconWithTooltip";
import ManageGateDialogWrapper from "./ManageGateDialogWrapper";
import PolygonIcon from "./PolygonIcon";

// the conventional icon for each selection tool, from open-licence sets (lucide: ISC)
export const Tools = {
    pan: {
        name: "Pan",
        ToolIcon: PanToolOutlinedIcon,
        mode: TranslateMode,
    },
    rectangle: {
        name: "Rectangle",
        ToolIcon: VectorSquareIcon,
        mode: RectangleMode,
    },
    // todo: add these back in once we have deck EditableGeoJsonLayer etc in place...
    polygon: {
        name: "Polygon",
        ToolIcon: PolygonIcon,
        mode: PolygonMode,
    },
    freehand: {
        name: "Freehand",
        ToolIcon: LassoIcon,
        mode: DrawPolygonByDraggingMode,
    },
    transform: {
        name: "Modify",
        ToolIcon: SplinePointerIcon,
        mode: EditMode,
    },
    // measure: {
    //     name: "Measure",
    //     ToolIcon: StraightenIcon,
    //     mode: DrawLineStringMode
    // },
} as const;
type ToolIcon = (typeof Tools)[keyof typeof Tools]["ToolIcon"];
export type Tool = (typeof Tools)[keyof typeof Tools]["name"];
const ToolArray = Object.values(Tools);
const EDIT_MODE_TOOL_NAMES: Tool[] = ["Pan", "Modify"];

type ToolButtonProps = {
    name: Tool;
    ToolIcon: ToolIcon;
    selectedTool: Tool;
    setSelectedTool: (tool: Tool) => void;
};

const ToolButton = observer(({ name, ToolIcon, selectedTool, setSelectedTool }: ToolButtonProps) => {
    const theme = useTheme();
    const selectedBackgroundColor = theme === "dark" ? "rgba(255,255,255,0.2)" : "rgba(48, 46, 46, 0.18)";
    const style = useMemo(
        () => ({
            backgroundColor: selectedTool === name ? selectedBackgroundColor : "transparent",
            color: "var(--text_color)",
            borderRadius: 10,
            zIndex: 2,
        }),
        [selectedTool, name, selectedBackgroundColor],
    );
    return (
        <IconWithTooltip
            tooltipText={name}
            onClick={() => setSelectedTool(name)}
            iconButtonProps={{
                sx: style,
                "aria-label": name,
            }}
        >
            <ToolIcon size={20} />
        </IconWithTooltip>
    );
});

export default observer(function SelectionOverlay() {
    const { selectionProps } = useSpatialLayers();
    const chart = useChart<DeckScatterConfigWithRegion>();
    const chartId = useChartID();
    const config = useConfig<DeckScatterConfigWithRegion & DualContourLegacyConfig>();
    const canShowDensityGrid = supportsDensityGridMode(config.type) && config.dimension !== "3d";
    // gates are x / y outlines in the data; a 3D selection is an outline on the screen
    const canUseGates = config.dimension !== "3d";
    const densityGridShown = isDensityGridShown(config);
    const toggleDensityGrid = useCallback(() => {
        action(() => {
            config.grid = densityGridConfig(!isDensityGridShown(config));
        })();
    }, [config]);
    const {
        setSelectionMode,
        selectionFeatureCollection,
        editingGateId,
        setSelectionFeatureCollection,
        selectedTool,
        setSelectedTool: setSelectedToolX,
    } = selectionProps;
    const [cx, cy] = useParamColumns();

    const [gateDialogOpen, setGateDialogOpen] = useState(false);
    // the toolbar starts collapsed to the active tool, so it covers less of the plot
    const [expanded, setExpanded] = useState(false);
    // editing a gate needs the edit tools in view
    const showAllTools = expanded || editingGateId != null;
    const manageGateByChartRef = useRef<Map<string, ManageGateDialogWrapper<DeckScatterConfigWithRegion>>>(new Map());
    const theme = useTheme();

    const {
        onDeleteGate,
        onExportClick,
        onRenameGate,
        onSaveGate,
        onEditGate,
        onConfirmEditGate,
        onCancelEditGate,
        onColorChange,
    } = useGateActions();

    useEffect(() => {
        // close the open dialog when unmounting the component
        return () => {
            const dlg = manageGateByChartRef.current.get(chartId);
            if (dlg) {
                setTimeout(() => dlg?.close(), 0);
            }
            manageGateByChartRef.current.delete(chartId);
        };
    }, [chartId]);

    const toolsToShow = useMemo(
        () => (editingGateId ? ToolArray.filter((t) => EDIT_MODE_TOOL_NAMES.includes(t.name)) : ToolArray),
        [editingGateId],
    );

    // biome-ignore lint/correctness/useExhaustiveDependencies: need to check that gate id is correct here?
    const setSelectedTool = useCallback(
        (tool: Tool) => {
            // pending refactor
            const mode = Object.values(Tools).find((t) => t.name === tool)?.mode;
            if (!mode) {
                console.error("no mode found for tool", tool);
                return;
            }

            //same composite mode order doesn't work for all tools, so making `mode()` be more explicit for each
            setSelectionMode(new mode() as GeoJsonEditMode);
            setSelectedToolX(tool);
        },
        [setSelectionMode, editingGateId, setSelectionFeatureCollection, setSelectedToolX],
    );
    // add a row of buttons to the top of the chart
    // rectangle, circle, polygon, lasso, magic wand, etc.
    // (thanks copilot, that may be over-ambitious)
    // It would be good to have a poly-line tool with draggable points, though.
    // Also a brush tool for painting on masks with variable radius.
    const toolButtons = useMemo(() => {
        const tools = showAllTools ? toolsToShow : toolsToShow.filter((t) => t.name === selectedTool);
        return tools.map(({ name, ToolIcon }) => (
            <ToolButton
                key={name}
                name={name}
                ToolIcon={ToolIcon}
                selectedTool={selectedTool}
                setSelectedTool={setSelectedTool}
            />
        ));
    }, [selectedTool, setSelectedTool, toolsToShow, showAllTools]);

    const onManageGatesClick = useCallback(() => {
        if (manageGateByChartRef.current.get(chartId)) {
            // a dialog already exists
            return;
        }

        // Create a new manage gates dialog for a chart
        const dialog = new ManageGateDialogWrapper({
            chart,
            xField: cx?.field,
            yField: cy?.field,
            onDelete: onDeleteGate,
            onColorChange,
            onRenameGate,
            onExportClick,
            onEdit: onEditGate,
            onClose: () => {
                // remove the entry of the dialog
                manageGateByChartRef.current.delete(chartId);
            },
            activeRegion: config?.region,
        });

        // create a new entry for the dialog for a chart
        manageGateByChartRef.current.set(chartId, dialog);
    }, [
        cx?.field,
        cy?.field,
        onDeleteGate,
        onEditGate,
        onColorChange,
        onRenameGate,
        onExportClick,
        chart,
        chartId,
        config?.region,
    ]);

    const hasSelection = useMemo(
        () => selectionFeatureCollection.features.length > 0,
        [selectionFeatureCollection.features.length],
    );

    return (
        <>
            <ButtonGroup
                variant="contained"
                aria-label="choose tool for manipulating view or selection"
                //moving this to the top right corner and absolute to avoid interfering with axes
                className={`z-[2] absolute top-0 right-0 opacity-90 ${showAllTools ? "p-2" : "p-1"}`}
                sx={{
                    backgroundColor: theme === "dark" ? "rgba(37,37,37,0.92)" : "rgba(245,245,245,0.92)",
                    borderRadius: 0,
                    borderBottomLeftRadius: "5px",
                }}
            >
                {!editingGateId && (
                    <IconWithTooltip
                        tooltipText={expanded ? "Hide tools" : "Show tools"}
                        onClick={() => setExpanded((e) => !e)}
                        iconButtonProps={{
                            sx: { color: "var(--text_color)", borderRadius: 10, zIndex: 2, px: 0.25 },
                            "aria-label": expanded ? "Hide tools" : "Show tools",
                            "aria-expanded": expanded,
                        }}
                    >
                        {expanded ? <ChevronRightIcon /> : <ChevronLeftIcon />}
                    </IconWithTooltip>
                )}
                {toolButtons}
                {showAllTools && (
                    <>
                        {canShowDensityGrid && (
                            <IconWithTooltip
                                tooltipText={
                                    densityGridShown
                                        ? "Show density fields on one plot"
                                        : "Show density fields as a grid"
                                }
                                onClick={toggleDensityGrid}
                                iconButtonProps={{
                                    sx: { color: "var(--text_color)", borderRadius: 10, zIndex: 2 },
                                    "aria-label": densityGridShown ? "Density overlay" : "Density grid",
                                }}
                            >
                                {densityGridShown ? <LayersOutlinedIcon /> : <GridViewOutlined />}
                            </IconWithTooltip>
                        )}
                        {canUseGates && (
                            <Divider
                                orientation="vertical"
                                sx={{
                                    color: "inherit",
                                    width: "5px",
                                    height: "35px",
                                    padding: "2px",
                                }}
                            />
                        )}
                        {canUseGates && (
                            <IconWithTooltip
                                tooltipText={"Manage Gates"}
                                onClick={onManageGatesClick}
                                iconButtonProps={{
                                    sx: {
                                        color: "var(--text_color)",
                                        borderRadius: 10,
                                        zIndex: 2,
                                        ml: 1,
                                    },
                                    "aria-label": "Manage Gates",
                                }}
                            >
                                <TuneOutlined />
                            </IconWithTooltip>
                        )}
                    </>
                )}
                {canUseGates && hasSelection && !editingGateId && (
                    <IconWithTooltip
                        tooltipText={"Save selection as Gate"}
                        onClick={() => setGateDialogOpen(true)}
                        iconButtonProps={{
                            "aria-label": "Save Gate",
                        }}
                    >
                        <AddCircleOutlinedIcon />
                    </IconWithTooltip>
                )}
                {editingGateId && (
                    <>
                        <IconWithTooltip
                            tooltipText={"Confirm"}
                            onClick={() => {
                                onConfirmEditGate().catch(console.error);
                            }}
                            iconButtonProps={{
                                sx: {
                                    color: "green",
                                    borderRadius: 10,
                                    zIndex: 2,
                                },
                                "aria-label": "Confirm",
                            }}
                        >
                            <DoneOutlinedIcon />
                        </IconWithTooltip>
                        <IconWithTooltip
                            tooltipText={"Cancel"}
                            onClick={onCancelEditGate}
                            iconButtonProps={{
                                sx: {
                                    color: "red",
                                    borderRadius: 10,
                                    zIndex: 2,
                                },
                                "aria-label": "Cancel",
                            }}
                        >
                            <CloseOutlinedIcon />
                        </IconWithTooltip>
                    </>
                )}
            </ButtonGroup>

            <GateNameDialog open={gateDialogOpen} onClose={() => setGateDialogOpen(false)} onSaveGate={onSaveGate} />
        </>
    );
});
