import {
    Delete,
    DriveFileRenameOutline,
    Image as ImageIcon,
    Info,
    MoreVert,
    Share as ShareIcon,
    ViewColumn,
    WarningAmber,
} from "@mui/icons-material";
import {
    Checkbox,
    IconButton,
    ListItemIcon,
    ListItemText,
    Menu,
    MenuItem,
    Paper,
    Table,
    TableBody,
    TableCell,
    TableContainer,
    TableHead,
    TableRow,
    TableSortLabel,
    Tooltip,
    Typography,
    useTheme,
} from "@mui/material";
import { useCallback, useState } from "react";
import { buildProjectUrl } from "@/utils/mdvRouting";
import usePermissions from "./PermissionsContext";
import ProjectDeleteModal from "./ProjectDeleteModal";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectRenameModal from "./ProjectRenameModal";
import ProjectShareModal from "./ProjectShareModal";
import { LIST_COLUMNS, type ListColumn } from "./utils/catalogPreferences";
import type { Project, ProjectAccessType, SortBy, SortOrder } from "./utils/projectUtils";

export type Pvoid = Promise<void>;
export type ProjectListViewProps = {
    projects: Project[];
    columns: ListColumn[];
    onColumnsChange: (columns: ListColumn[]) => void;
    sortBy: SortBy;
    sortOrder: SortOrder;
    onSortChange: (sortBy: SortBy) => void;
    onDelete: (id: string) => Pvoid;
    onRename: (id: string, name: string) => Pvoid;
    onChangeType: (id: string, type: ProjectAccessType) => Pvoid;
    onExport: (id: string, name: string) => Pvoid;
    selectionMode: boolean;
    selectedProjectIds: Set<string>;
    onSelectionChange: (id: string, selected: boolean) => void;
};
export type MouseEv = React.MouseEvent<HTMLElement>;
const ProjectListView = ({
    projects,
    columns,
    onColumnsChange,
    sortBy,
    sortOrder,
    onSortChange,
    onDelete,
    onRename,
    onExport,
    onChangeType,
    selectionMode,
    selectedProjectIds,
    onSelectionChange,
}: ProjectListViewProps) => {
    const [columnsAnchorEl, setColumnsAnchorEl] = useState<HTMLElement | null>(null);
    const [anchorEl, setAnchorEl] = useState<HTMLElement | null>();
    const [selectedProject, setSelectedProject] = useState<Project>();

    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    const [isRenameModalOpen, setIsRenameModalOpen] = useState(false);
    const [isDeleteModalOpen, setIsDeleteModalOpen] = useState(false);
    const [isShareModalOpen, setIsShareModalOpen] = useState(false);
    const theme = useTheme();
    const { permissions: operationPermissions, isPublicPage } = usePermissions();

    const handleMenuClick = useCallback((event: MouseEv, project: Project) => {
        event.stopPropagation();
        setAnchorEl(event.currentTarget);
        setSelectedProject(project);
    }, []);

    const handleMenuClose = useCallback(() => {
        setAnchorEl(null);
    }, []);

    const handleRowClick = useCallback((projectId: Project["id"]) => {
        window.location.href = buildProjectUrl(projectId);
    }, []);

    const handleModalOpen = useCallback(
        (modalSetter: React.Dispatch<React.SetStateAction<boolean>>) => {
            handleMenuClose();
            modalSetter(true);
        },
        [handleMenuClose],
    );

    const show = useCallback(
        (col: ListColumn) =>
            columns.includes(col) && !(isPublicPage && (col === "lastModified" || col === "filesystem")),
        [columns, isPublicPage],
    );

    const ariaSort = (key: SortBy) => (sortBy === key ? (sortOrder === "asc" ? "ascending" : "descending") : "none");

    const hasReadme = useCallback((proj: Project) => !!proj.readme, []);

    const hasPermissions = useCallback(
        (proj: Project) =>
            proj.permissions.edit &&
            (operationPermissions.renameProject ||
                operationPermissions.deleteProject ||
                // operationPermissions.exportProject ||
                operationPermissions.shareProject),
        [operationPermissions],
    );

    // todo: Update this component for more cleaner code
    return (
        <>
            <TableContainer component={Paper}>
                <Table sx={{ "& td:not(.name-cell), & th:not(.name-cell)": { whiteSpace: "nowrap" } }}>
                    <TableHead>
                        <TableRow>
                            {selectionMode && <TableCell padding="checkbox" />}
                            <TableCell className="name-cell" sx={{ width: "100%" }} aria-sort={ariaSort("name")}>
                                <TableSortLabel
                                    active={sortBy === "name"}
                                    direction={sortBy === "name" ? sortOrder : "asc"}
                                    onClick={() => onSortChange("name")}
                                >
                                    Name
                                </TableSortLabel>
                            </TableCell>
                            {show("owner") && <TableCell>Owner</TableCell>}
                            {show("lastModified") && (
                                <TableCell aria-sort={ariaSort("lastModified")}>
                                    <TableSortLabel
                                        active={sortBy === "lastModified"}
                                        direction={sortBy === "lastModified" ? sortOrder : "desc"}
                                        onClick={() => onSortChange("lastModified")}
                                    >
                                        Last Modified
                                    </TableSortLabel>
                                </TableCell>
                            )}
                            {show("type") && <TableCell>Type</TableCell>}
                            {show("filesystem") && <TableCell>Filesystem</TableCell>}
                            <TableCell align="right" sx={{ whiteSpace: "nowrap" }}>
                                Actions
                                <Tooltip title="Choose columns">
                                    <IconButton
                                        size="small"
                                        onClick={(e) => setColumnsAnchorEl(e.currentTarget)}
                                        aria-label="choose columns"
                                        data-testid="choose_columns_button"
                                        sx={{ ml: 0.5 }}
                                    >
                                        <ViewColumn fontSize="small" />
                                    </IconButton>
                                </Tooltip>
                            </TableCell>
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {projects.map((project) => (
                            <TableRow
                                key={project.id}
                                hover
                                onClick={() => handleRowClick(project.id)}
                                sx={{ cursor: "pointer" }}
                                data-testid="project_list_row"
                                data-project-id={project.id}
                                aria-label={`open project ${project.name} (${project.id})`}
                            >
                                {selectionMode && (
                                    <TableCell padding="checkbox">
                                        {operationPermissions.deleteProject &&
                                            project.permissions.edit &&
                                            project.permissions.owner && (
                                                <Checkbox
                                                    checked={selectedProjectIds.has(project.id)}
                                                    onChange={(event) =>
                                                        onSelectionChange(project.id, event.target.checked)
                                                    }
                                                    onClick={(event) => event.stopPropagation()}
                                                    inputProps={{ "aria-label": `select project ${project.name}` }}
                                                />
                                            )}
                                    </TableCell>
                                )}
                                <TableCell className="name-cell">
                                    <div className="flex items-center gap-3">
                                        {project.thumbnail ? (
                                            <img
                                                src={project.thumbnail}
                                                alt=""
                                                style={{ width: 48, height: 32, objectFit: "contain", flexShrink: 0 }}
                                            />
                                        ) : (
                                            <ImageIcon sx={{ width: 48, color: "text.secondary" }} />
                                        )}
                                        <Typography>{project.name}</Typography>
                                    </div>
                                </TableCell>
                                {show("owner") && (
                                    <TableCell>{project.owner ? project.owner.join(", ") : ""}</TableCell>
                                )}
                                {show("lastModified") && <TableCell>{project.lastModified}</TableCell>}
                                {show("type") && <TableCell>{project.type}</TableCell>}
                                {show("filesystem") && (
                                    <TableCell>
                                        {project.writable ? (
                                            "Writable"
                                        ) : (
                                            <Tooltip title="The MDV server cannot write this project's files. It will open read-only.">
                                                <Typography
                                                    component="span"
                                                    color="warning.main"
                                                    sx={{ display: "inline-flex", alignItems: "center", gap: 0.5 }}
                                                >
                                                    <WarningAmber fontSize="small" />
                                                    Server read-only
                                                </Typography>
                                            </Tooltip>
                                        )}
                                    </TableCell>
                                )}
                                <TableCell align="right">
                                    {(hasReadme(project) || hasPermissions(project)) && (
                                        <IconButton
                                            onClick={(e) => handleMenuClick(e, project)}
                                            size="small"
                                            data-testid={`project_menu_${project.id}`}
                                        >
                                            <MoreVert />
                                        </IconButton>
                                    )}
                                </TableCell>
                            </TableRow>
                        ))}
                    </TableBody>
                </Table>
            </TableContainer>

            <Menu anchorEl={columnsAnchorEl} open={Boolean(columnsAnchorEl)} onClose={() => setColumnsAnchorEl(null)}>
                {LIST_COLUMNS.filter((c) => !(isPublicPage && c.hiddenOnPublicPage)).map((c) => (
                    <MenuItem
                        key={c.id}
                        dense
                        onClick={() =>
                            onColumnsChange(
                                columns.includes(c.id) ? columns.filter((x) => x !== c.id) : [...columns, c.id],
                            )
                        }
                    >
                        <Checkbox checked={columns.includes(c.id)} size="small" disableRipple />
                        {c.label}
                    </MenuItem>
                ))}
            </Menu>

            {selectedProject && (
                <>
                    <Menu anchorEl={anchorEl} open={Boolean(anchorEl)} onClose={handleMenuClose}>
                        {hasReadme(selectedProject) && (
                            <MenuItem onClick={() => handleModalOpen(setIsInfoModalOpen)}>
                                <ListItemIcon>
                                    <Info fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Project Info</ListItemText>
                            </MenuItem>
                        )}
                        {operationPermissions.renameProject && (
                            <MenuItem
                                onClick={() => handleModalOpen(setIsRenameModalOpen)}
                                disabled={!selectedProject?.permissions.edit}
                            >
                                <ListItemIcon>
                                    <DriveFileRenameOutline fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Rename Project</ListItemText>
                            </MenuItem>
                        )}
                        {/* {operationPermissions.exportProject && (
                            <MenuItem 
                                onClick={() => {
                                    if (selectedProject) {
                                        handleMenuClose();
                                        onExport(selectedProject.id, selectedProject.name);
                                    }
                                }}
                                disabled={!selectedProject?.permissions.edit}
                            >
                                <ListItemIcon>
                                    <UploadIcon fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Export Project (as *.mdv.zip)</ListItemText>
                            </MenuItem>
                        )} */}
                        {operationPermissions.shareProject && (
                            <MenuItem
                                onClick={() => handleModalOpen(setIsShareModalOpen)}
                                disabled={!selectedProject?.permissions.edit}
                            >
                                <ListItemIcon>
                                    <ShareIcon fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Share Project</ListItemText>
                            </MenuItem>
                        )}
                        {operationPermissions.deleteProject && (
                            <MenuItem
                                onClick={() => handleModalOpen(setIsDeleteModalOpen)}
                                sx={{ color: theme.palette.error.main }}
                                disabled={!selectedProject?.permissions.edit}
                            >
                                <ListItemIcon>
                                    <Delete color="error" fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Move to Recycle Bin</ListItemText>
                            </MenuItem>
                        )}
                    </Menu>

                    {hasReadme(selectedProject) && (
                        <ProjectInfoModal
                            open={isInfoModalOpen}
                            onClose={() => setIsInfoModalOpen(false)}
                            readme={selectedProject.readme}
                        />
                    )}

                    {selectedProject.permissions.edit && (
                        <>
                            {operationPermissions.renameProject && (
                                <ProjectRenameModal
                                    id={selectedProject.id}
                                    name={selectedProject.name}
                                    open={isRenameModalOpen}
                                    onRename={onRename}
                                    onClose={() => setIsRenameModalOpen(false)}
                                />
                            )}

                            {operationPermissions.shareProject && (
                                <ProjectShareModal
                                    open={isShareModalOpen}
                                    onClose={() => setIsShareModalOpen(false)}
                                    projectId={selectedProject.id}
                                />
                            )}

                            {operationPermissions.deleteProject && (
                                <ProjectDeleteModal
                                    id={selectedProject.id}
                                    open={isDeleteModalOpen}
                                    onDelete={onDelete}
                                    onClose={() => setIsDeleteModalOpen(false)}
                                />
                            )}
                        </>
                    )}
                </>
            )}
        </>
    );
};

export default ProjectListView;
