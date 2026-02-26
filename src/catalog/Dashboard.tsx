import { useColorMode } from "@/ThemeProvider";
import {
    Add,
    Download as DownloadIcon,
    ExpandMore,
    GridView,
    Reorder as ReorderIcon,
    Search,
} from "@mui/icons-material";
import Brightness4Icon from "@mui/icons-material/Brightness4";
import Brightness7Icon from "@mui/icons-material/Brightness7";
import {
    AppBar,
    Backdrop,
    Box,
    Button,
    ButtonBase,
    CircularProgress,
    Container,
    Divider,
    Grid2 as Grid,
    IconButton,
    InputBase,
    Menu,
    MenuItem,
    Paper,
    Toolbar,
    Tooltip,
    Typography,
    useTheme,
} from "@mui/material";
import React, { useCallback, useEffect, useMemo, useState } from "react";
import ProjectCard from "./ProjectCard";
import ProjectListView from "./ProjectListView";
import UserProfile from "./UserProfile";
import mdvLogo from "./assets/mdv_logo.png";
import useProjects from "./hooks/useProjects";
import {
    type SortBy,
    type SortOrder,
    sortProjects,
} from "./utils/projectUtils";
import AlertErrorComponent from "@/charts/dialogs/AlertErrorComponent";
import ImportProjectDialog from "@/react/components/ImportProjectDialog";
import usePermissions from "./PermissionsContext";
import useAuthEnabled from "./hooks/useAuthEnabled";
import { RefreshCwIcon } from "lucide-react";
import ReusableAlertDialog from "@/charts/dialogs/ReusableAlertDialog";

// todo: Refactor the code into different components and hooks for cleaner and readable code
// Maybe use a design pattern? As displaying certain components depend on some states
const Dashboard: React.FC = () => {
    const {
        projects,
        isLoading: projectsLoading,
        error,
        isErrorModalOpen,
        closeErrorModal,
        fetchProjects,
        createProject,
        deleteProject,
        renameProject,
        changeProjectType,
        setFilter,
        exportProject,
        rescanProjects,
    } = useProjects();
    const { permissions, isLoading: permissionsLoading, isPublicPage } = usePermissions();

    // Check if auth is enabled
    const authEnabled = useAuthEnabled();
    const { mode, toggleColorMode } = useColorMode();
    const [viewMode, setViewMode] = useState<"grid" | "list">("grid");
    const [sortBy, setSortBy] = useState<SortBy>("lastModified");
    const [sortOrder, setSortOrder] = useState<SortOrder>("desc");
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
    const [open, setOpen] = useState(false);
    const theme = useTheme();

    const isLoading = projectsLoading || permissionsLoading;

    React.useEffect(() => {
        fetchProjects();
    }, [fetchProjects]);

    useEffect(() => {
        if (!permissionsLoading && isPublicPage) {
            setSortBy("name");
            setSortOrder("asc");
        }
    }, [isPublicPage, permissionsLoading]);

    const handleCreateProject = async () => {
        try {
            const newProject = await createProject();
            const base = import.meta.env.DEV
                ? "http://localhost:5170?dir=/"
                : "";
            window.location.href = `${base}project/${newProject?.id}`;
        } catch (error) {
            console.error("Failed to create project:", error);
        }
    };

    const sortedProjects = useMemo(() => {
        return sortProjects(projects, sortBy, sortOrder);
    }, [projects, sortBy, sortOrder]);

    const handleSort = useCallback(
        (newSortBy: SortBy) => {
            if (newSortBy === sortBy) {
                // Toggle sort order if clicking the same sort option
                setSortOrder(sortOrder === "asc" ? "desc" : "asc");
            } else {
                // Set new sort type with default desc order for both cases
                setSortBy(newSortBy);
                setSortOrder("desc");
            }
            setAnchorEl(null);
        },
        [sortBy, sortOrder],
    );

    const toggleDropdown = (event: React.MouseEvent<HTMLButtonElement>) => {
        setAnchorEl(event.currentTarget);
    };

    const onRescanClick = useCallback(
        async () => {
            await rescanProjects();
        }, [rescanProjects]
    );

    return (
        <>
            <Box
                sx={{
                    flexGrow: 1,
                    bgcolor: "background.default",
                    minHeight: "100vh",
                }}
            >
                <AppBar position="static" color="default" elevation={0}>
                    <Toolbar>
                        <Box
                            component="img"
                            sx={{
                                height: 40,
                                mr: 2,
                            }}
                            alt="MDV Projects Logo"
                            src={mdvLogo}
                        />
                        <Box sx={{ flexGrow: 1 }} />
                        <Paper
                            component="form"
                            sx={{
                                p: "2px 4px",
                                display: "flex",
                                alignItems: "center",
                                width: 400,
                                mr: 2, // Add margin to the right
                            }}
                        >
                            <InputBase
                                sx={{ ml: 1, flex: 1 }}
                                placeholder="Search projects"
                                inputProps={{ "aria-label": "search projects" }}
                                onChange={(e) => setFilter(e.target.value)}
                            />
                            <IconButton
                                type="submit"
                                sx={{ p: "10px" }}
                                aria-label="search"
                            >
                                <Search />
                            </IconButton>
                        </Paper>
                        {authEnabled && <UserProfile />}
                        <IconButton
                            sx={{ ml: 1 }}
                            onClick={toggleColorMode}
                            color="inherit"
                            data-testid="theme_toggle_catalog"
                        >
                            {mode === "dark" ? (
                                <Brightness4Icon />
                            ) : (
                                <Brightness7Icon />
                            )}
                        </IconButton>
                    </Toolbar>
                </AppBar>

                <Container maxWidth="lg" sx={{ mt: 4, mb: 4 }}>
                    <Grid container spacing={3} sx={{ mb: 4 }}>
                        {permissions.createProject && (
                            <Grid size={{xs: 12, sm: 6, md: 4, lg: 3}}>
                                <ButtonBase
                                    sx={{
                                        width: "100%",
                                        display: "block",
                                        textAlign: "center",
                                    }}
                                >
                                    <Paper
                                        sx={{
                                            p: 2,
                                            display: "flex",
                                            flexDirection: "column",
                                            alignItems: "center",
                                        }}
                                        onClick={() => handleCreateProject()}
                                    >
                                        <Add
                                            sx={{
                                                fontSize: 40,
                                                color: "primary.main",
                                                mb: 1,
                                            }}
                                        />
                                        <Typography variant="subtitle1" align="center">
                                            Create new project
                                        </Typography>
                                    </Paper>
                                </ButtonBase>
                            </Grid>
                        )}
                        {permissions.importProject && (
                            <Grid size={{xs: 12, sm: 6, md: 4, lg: 3}}>
                                <ButtonBase
                                    sx={{
                                        width: "100%",
                                        display: "block",
                                        textAlign: "center",
                                    }}
                                >
                                    <Paper
                                        sx={{
                                            p: 2,
                                            display: "flex",
                                            flexDirection: "column",
                                            alignItems: "center",
                                        }}
                                        onClick={() => setOpen(true)}
                                    >
                                        <DownloadIcon
                                            sx={{
                                                fontSize: 40,
                                                color: "primary.main",
                                                mb: 1,
                                            }}
                                        />
                                        <Typography variant="subtitle1" align="center" sx={{}}>
                                            Import an existing project
                                        </Typography>
                                    </Paper>
                                </ButtonBase>
                            </Grid>
                        )}
                    </Grid>

                    <Box
                        sx={{
                            display: "flex",
                            justifyContent: "space-between",
                            alignItems: "center",
                            mb: 1,
                        }}
                    >
                        <Typography variant="h5">{!isPublicPage ? "Recent Projects" : "Published Projects"}</Typography>
                        <Box sx={{ display: "flex", alignItems: "center" }}>
                            {/* Hide rescan projects on public page */}
                            {!isPublicPage && (
                                <Paper
                                    elevation={1}
                                    sx={{
                                        display: "flex",
                                        alignItems: "center",
                                        justifyContent: "center",
                                        borderRadius: "4px",
                                        bgcolor: "background.paper",
                                        marginRight: 2,
                                        height: "50px",
                                    }}
                                >
                                    <Button
                                        sx={{
                                            padding: 1,
                                            display: "flex",
                                            justifyContent: "space-around",
                                            height: "100%",
                                            width: "100%"
                                        }}
                                        onClick={onRescanClick}
                                    >
                                        <RefreshCwIcon />
                                        <Typography 
                                            sx={{
                                                marginLeft: 1
                                            }}
                                        >
                                            Rescan Projects
                                        </Typography>
                                    </Button>
                                </Paper>
                            )}
                            <Paper
                                elevation={1}
                                sx={{
                                    display: "flex",
                                    alignItems: "center",
                                    justifyContent: "center",
                                    borderRadius: "4px",
                                    width: "205px",
                                    height: "50px",
                                    bgcolor: "background.paper",
                                }}
                            >
                                <Button
                                    endIcon={<ExpandMore />}
                                    onClick={toggleDropdown}
                                    sx={{
                                        padding: 1,
                                        textTransform: "none",
                                        width: "100%",
                                        height: "100%",
                                        display: "flex",
                                        justifyContent: "space-between",
                                        alignItems: "center",
                                    }}
                                >
                                    <Box
                                        sx={{
                                            display: "flex",
                                            alignItems: "center",
                                        }}
                                    >
                                        {/* Hide last modified on public page */}
                                        Sort by:{" "}
                                        {!isPublicPage ? 
                                            sortBy === "lastModified"
                                                ? `Last modified (${sortOrder === "desc" ? "Newest first" : "Oldest first"})`
                                                : `Name (${sortOrder === "desc" ? "Z to A" : "A to Z"})`
                                            :
                                            `Name (${sortOrder === "desc" ? "Z to A" : "A to Z"})`
                                        }
                                    </Box>
                                </Button>
                            </Paper>
                            <Tooltip
                                title={
                                    viewMode === "grid" ? "List View" : "Grid View"
                                }
                            >
                                <IconButton
                                    onClick={() =>
                                        setViewMode(
                                            viewMode === "grid" ? "list" : "grid",
                                        )
                                    }
                                    sx={{ ml: 2 }}
                                >
                                    {viewMode === "grid" ? (
                                        <ReorderIcon sx={{ fontSize: 32 }} />
                                    ) : (
                                        <GridView sx={{ fontSize: 32 }} />
                                    )}
                                </IconButton>
                            </Tooltip>
                        </Box>

                        <Menu
                            anchorEl={anchorEl}
                            open={Boolean(anchorEl)}
                            onClose={() => setAnchorEl(null)}
                        >
                            {!isPublicPage && 
                                (<MenuItem
                                    onClick={() => handleSort("lastModified")}
                                    sx={{
                                        display: "flex",
                                        justifyContent: "space-between",
                                        width: "200px",
                                        gap: 1,
                                    }}
                                >
                                    <span>Last modified</span>
                                    {sortBy === "lastModified" && (
                                        <span>
                                            {sortOrder === "desc"
                                                ? "(Oldest first)"
                                                : "(Newest first)"
                                            }
                                        </span>
                                    )}
                                </MenuItem>)
                            }
                            <MenuItem
                                onClick={() => handleSort("name")}
                                sx={{
                                    display: "flex",
                                    justifyContent: "space-between",
                                    width: "200px",
                                    gap: 1,
                                }}
                            >
                                <span>Name</span>
                                {sortBy === "name" && (
                                    <span>
                                        {sortOrder === "desc" ? "A to Z" : "Z to A"}
                                    </span>
                                )}
                            </MenuItem>
                        </Menu>
                    </Box>

                    <Divider sx={{ mb: 2 }} />

                    {viewMode === "grid" ? (
                        <Grid container spacing={4}>
                            {sortedProjects.map((project) => (
                                <Grid
                                    key={project.id}
                                    size={{
                                        xs: 12,
                                        sm: 6,
                                        md: 4,
                                        lg: 3
                                    }}
                                >
                                    <ProjectCard
                                        {...project}
                                        onDelete={deleteProject}
                                        onRename={renameProject}
                                        onChangeType={changeProjectType}
                                        onAddCollaborator={(email) => {}}
                                        onExport={exportProject}
                                    />
                                </Grid>
                            ))}
                        </Grid>
                    ) : (
                        <ProjectListView
                            projects={sortedProjects}
                            onDelete={deleteProject}
                            onRename={renameProject}
                            onChangeType={changeProjectType}
                            onExport={exportProject}
                        />
                    )}
                </Container>
                {error && (
                    <ReusableAlertDialog
                        open={isErrorModalOpen}
                        handleClose={closeErrorModal}
                        isAlertErrorComponent
                        component={
                            <AlertErrorComponent
                                message={error as string}
                            />
                        }
                    />
                )}
                {open && permissions.importProject && (
                    <ImportProjectDialog open={open} setOpen={setOpen} />
                )}
            </Box>
            {isLoading && (
                <Backdrop
                    open={isLoading}
                    sx={{ zIndex: theme.zIndex.modal + 1 }}
                >
                    <CircularProgress color="inherit" />
                </Backdrop>
            )}
        </>
    );
};

export default Dashboard;
