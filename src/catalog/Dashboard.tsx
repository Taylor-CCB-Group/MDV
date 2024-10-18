import React, { useCallback, useState } from "react";
import {
    Add,
    Search,
    Menu as MenuIcon,
    GridView,
    ViewList,
    ExpandMore,
    Folder,
} from "@mui/icons-material";
import {
    AppBar,
    Toolbar,
    Typography,
    IconButton,
    InputBase,
    Paper,
    Button,
    Grid,
    Menu,
    MenuItem,
    Tooltip,
    Box,
    Container,
    Divider,
    useTheme,
    ThemeProvider,
    createTheme,
    ButtonBase,
    CircularProgress,
} from "@mui/material";
import Brightness4Icon from "@mui/icons-material/Brightness4";
import Brightness7Icon from "@mui/icons-material/Brightness7";
import { alpha } from "@mui/material/styles";
import ProjectCard from "./ProjectCard";
import useProjects from "./hooks/useProjects";
import UserProfile from "./UserProfile";
import mdvLogo from "./assets/mdv_logo.png";
import { useColorMode } from "@/ThemeProvider";


const Dashboard: React.FC = () => {
    const {
        projects,
        isLoading,
        error,
        fetchProjects,
        createProject,
        deleteProject,
        renameProject,
        changeProjectType,
        setFilter,
        setSortBy,
        setSortOrder,
        sortBy,
    } = useProjects();

    const { mode, toggleColorMode } = useColorMode();
    const [viewMode, setViewMode] = useState<"grid" | "list">("grid");
    const [projectType, setProjectType] = useState<"Editable" | "Read-Only">(
        "Editable",
    );
    const [isDropdownOpen, setIsDropdownOpen] = useState(false);
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);

    React.useEffect(() => {
        fetchProjects();
    }, [fetchProjects]);

    const handleCreateProject = async () => {
        // in future we should change the signature of createProject to not need the project name...
        // mdv_server_app will ignore it anyway... but for now I'm making it so in dev server we can enter the project name with a prompt...
        const newProjectName = import.meta.env.DEV ? prompt("Enter project name") : "<<UNUSED>>";
        if (newProjectName.trim()) {
            try {
                const newProject = await createProject(newProjectName.trim());
                // Redirect to the new project page
                // note that we may not be hosted in the root of the domain, so we should use relative paths
                window.location.href = `project/${newProject.id}`;
            } catch (error) {
                console.error("Failed to create project:", error);
                alert("Failed to create project. Please try again.");
            }
        }
    };

    const handleSort = useCallback((newSortBy: "lastModified" | "name") => {
        // toggle sort order if changing sort option
        if (newSortBy === sortBy) setSortOrder((prev) => (prev === "asc" ? "desc" : "asc"));
        // default to ascending order if changing sort option to name, descending order otherwise (most recent first)
        else setSortOrder(newSortBy === "name" ? "asc" : "desc");
        setSortBy(newSortBy);
        setIsDropdownOpen(false);
    }, [sortBy]);

    const toggleDropdown = (event: React.MouseEvent<HTMLButtonElement>) => {
        setAnchorEl(event.currentTarget);
        setIsDropdownOpen(!isDropdownOpen);
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
                                <Brightness7Icon />
                            ) : (
                                <Brightness4Icon />
                            )}
                        </IconButton>
                    </Toolbar>
                </AppBar>

                <Container maxWidth="lg" sx={{ mt: 4 }}>
                    <Grid container spacing={3} sx={{ mb: 4 }}>
                        <Grid item xs={12} sm={6} md={3}>
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
                                    <Typography
                                        variant="subtitle1"
                                        align="center"
                                    >
                                        Create new project
                                    </Typography>
                                </Paper>
                            </ButtonBase>
                        </Grid>

                        <Grid item xs={12} sm={6} md={3}>
                            <ButtonBase
                                onClick={() => {
                                    // Handle the create from template action
                                }}
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
                                >
                                    <Folder
                                        sx={{
                                            fontSize: 40,
                                            color: "secondary.main",
                                            mb: 1,
                                        }}
                                    />
                                    <Typography
                                        variant="subtitle1"
                                        align="center"
                                    >
                                        Create from template 1
                                    </Typography>
                                </Paper>
                            </ButtonBase>
                        </Grid>

                        <Grid item xs={12} sm={6} md={3}>
                            <ButtonBase
                                onClick={() => {}}
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
                                >
                                    <Folder
                                        sx={{
                                            fontSize: 40,
                                            color: "secondary.main",
                                            mb: 1,
                                        }}
                                    />
                                    <Typography
                                        variant="subtitle1"
                                        align="center"
                                    >
                                        Create from template 2
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
                                    justifyContent: "center",
                                }}
                            >
                                Sort by:{" "}
                                {sortBy === "lastModified"
                                    ? "Last modified"
                                    : "Name"}
                            </Button>
                        </Paper>
                        <Menu
                            anchorEl={anchorEl}
                            open={isDropdownOpen}
                            onClose={() => setIsDropdownOpen(false)}
                        >
                            <MenuItem
                                onClick={() => handleSort("lastModified")}
                            >
                                Last modified
                            </MenuItem>
                            <MenuItem onClick={() => handleSort("name")}>
                                Name
                            </MenuItem>
                        </Menu>
                    </Box>

                    <Divider sx={{ mb: 2 }} />

                    {isLoading ? (
                        <CircularProgress />
                    ) : error ? (
                        <Typography color="error">{error}</Typography>
                    ) : (
                        <Grid container spacing={4}>
                            {projects.map((project) => (
                                <Grid
                                    item
                                    key={project.id}
                                    xs={12}
                                    sm={6}
                                    md={4}
                                    lg={3}
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
                    )}
                </Container>
            </Box>
    );
};

export default Dashboard;
