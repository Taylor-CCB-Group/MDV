import { useColorMode } from "@/ThemeProvider";
import {
    Add,
    Cloud,
    CloudUpload,
    Download as DownloadIcon,
    ExpandMore,
    Folder,
    GridView,
    Reorder as ReorderIcon,
    Search,
    Upload,
} from "@mui/icons-material";
import Brightness4Icon from "@mui/icons-material/Brightness4";
import Brightness7Icon from "@mui/icons-material/Brightness7";
import {
    AppBar,
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
} from "@mui/material";
import React, { useCallback, useMemo, useState } from "react";
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
import ReusableDialog from "@/charts/dialogs/ReusableDialog";
import AlertErrorComponent from "@/charts/dialogs/AlertErrorComponent";
import ImportProjectDialog from "@/react/components/ImportProjectDialog";

const Dashboard: React.FC = () => {
    const {
        projects,
        isLoading,
        error,
        isErrorModalOpen,
        closeErrorModal,
        fetchProjects,
        createProject,
        deleteProject,
        renameProject,
        changeProjectType,
        setFilter,
    } = useProjects();

    const { mode, toggleColorMode } = useColorMode();
    const [viewMode, setViewMode] = useState<"grid" | "list">("grid");
    const [sortBy, setSortBy] = useState<SortBy>("lastModified");
    const [sortOrder, setSortOrder] = useState<SortOrder>("desc");
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
    const [open, setOpen] = useState(false);

    React.useEffect(() => {
        fetchProjects();
    }, [fetchProjects]);

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

    return (
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
                    <UserProfile />
                    <IconButton
                        sx={{ ml: 1 }}
                        onClick={toggleColorMode}
                        color="inherit"
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
                    <Grid size={{xs: 12, sm: 6, md: 3}}>
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
                    <Grid size={{xs: 12, sm: 6, md: 3}}>
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
                                <Typography variant="subtitle1" align="center">
                                    Import an existing project
                                </Typography>
                            </Paper>
                        </ButtonBase>
                    </Grid>
                </Grid>

                <Box
                    sx={{
                        display: "flex",
                        justifyContent: "space-between",
                        alignItems: "center",
                        mb: 1,
                    }}
                >
                    <Typography variant="h5">Recent Projects</Typography>
                    <Box sx={{ display: "flex", alignItems: "center" }}>
                        <Paper
                            elevation={1}
                            sx={{
                                padding: "8px",
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
                                    Sort by:{" "}
                                    {sortBy === "lastModified"
                                        ? `Last modified (${sortOrder === "desc" ? "Newest first" : "Oldest first"})`
                                        : `Name (${sortOrder === "desc" ? "Z to A" : "A to Z"})`}
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
                        <MenuItem
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
                                        ? "Newest first"
                                        : "Oldest first"}
                                </span>
                            )}
                        </MenuItem>
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
                                    {sortOrder === "desc" ? "Z to A" : "A to Z"}
                                </span>
                            )}
                        </MenuItem>
                    </Menu>
                </Box>

                <Divider sx={{ mb: 2 }} />

                {isLoading ? (
                    <CircularProgress />
                ) : viewMode === "grid" ? (
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
                    />
                )}
            </Container>
            <ReusableDialog
                open={isErrorModalOpen}
                handleClose={closeErrorModal}
                isAlertErrorComponent
                component={
                    <AlertErrorComponent
                        message={error as string}
                    />
                }
            />
            <ImportProjectDialog open={open} setOpen={setOpen} />
        </Box>
    );
};

export default Dashboard;
