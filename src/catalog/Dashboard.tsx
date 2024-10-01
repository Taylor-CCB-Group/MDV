import React, { useState } from "react";
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
import { alpha } from "@mui/material/styles";
import ProjectCard from "./ProjectCard";
import useProjects from "./useProjects";

const theme = createTheme({
    palette: {
        primary: {
            main: "#2c3e50",
        },
        secondary: {
            main: "#34495e",
        },
        background: {
            default: "#ecf0f1",
            paper: "#ffffff",
        },
        text: {
            primary: "#2c3e50",
            secondary: "#7f8c8d",
        },
    },
    typography: {
        fontFamily: '"Roboto", "Helvetica", "Arial", sans-serif',
        h5: {
            fontWeight: 500,
        },
        body1: {
            fontSize: "0.9rem",
        },
    },
    components: {
        MuiButton: {
            styleOverrides: {
                root: {
                    textTransform: "none",
                },
            },
        },
        MuiPaper: {
            styleOverrides: {
                root: {
                    boxShadow:
                        "0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24)",
                },
            },
        },
    },
});

const Dashboard: React.FC = () => {
    const {
        projects,
        isLoading,
        error,
        fetchProjects,
        createProject,
        deleteProject,
        setFilter,
        setSortBy,
        setSortOrder,
        sortBy,
    } = useProjects();

    const [viewMode, setViewMode] = useState<"grid" | "list">("grid");
    const [isModalOpen, setIsModalOpen] = useState(false);
    const [newProjectName, setNewProjectName] = useState("");
    const [projectType, setProjectType] = useState<"Editable" | "Read-Only">(
        "Editable",
    );
    const [isDropdownOpen, setIsDropdownOpen] = useState(false);
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);

    React.useEffect(() => {
        fetchProjects();
    }, [fetchProjects]);

    const handleCreateProject = async () => {
        if (newProjectName.trim()) {
            try {
                const newProject = await createProject(newProjectName.trim());
                setNewProjectName("");
                setIsModalOpen(false);
                // Redirect to the new project page
                window.location.href = `/project/${newProject.id}`;
            } catch (error) {
                console.error("Failed to create project:", error);
                alert("Failed to create project. Please try again.");
            }
        }
    };

    const handleSort = (option: "lastModified" | "name") => {
        setSortBy(option);
        setSortOrder((prev) => (prev === "asc" ? "desc" : "asc"));
        setIsDropdownOpen(false);
    };

    const toggleDropdown = (event: React.MouseEvent<HTMLButtonElement>) => {
        setAnchorEl(event.currentTarget);
        setIsDropdownOpen(!isDropdownOpen);
    };

    return (
        <ThemeProvider theme={theme}>
            <Box
                sx={{
                    flexGrow: 1,
                    bgcolor: "background.default",
                    minHeight: "100vh",
                }}
            >
                <AppBar position="static" color="default" elevation={0}>
                    <Toolbar>
                        <IconButton
                            size="large"
                            edge="start"
                            color="inherit"
                            aria-label="menu"
                            sx={{ mr: 2 }}
                        >
                            <MenuIcon />
                        </IconButton>
                        <Typography
                            variant="h6"
                            component="div"
                            sx={{ flexGrow: 1 }}
                        >
                            MDV Projects
                        </Typography>
                        <Paper
                            component="form"
                            sx={{
                                p: "2px 4px",
                                display: "flex",
                                alignItems: "center",
                                width: 400,
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
                    </Toolbar>
                </AppBar>

                <Container maxWidth="lg" sx={{ mt: 4 }}>
                    <Grid container spacing={3} sx={{ mb: 4 }}>
                        <Grid item xs={12} sm={6} md={3}>
                            <ButtonBase
                                onClick={() => setIsModalOpen(true)}
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
                                        onRename={(id, newName) => {
                                        }}
                                        onChangeType={(id, newType) => {
                                        }}
                                        onAddCollaborator={(email) => {
                                        }}
                                    />
                                </Grid>
                            ))}
                        </Grid>
                    )}
                </Container>
            </Box>
        </ThemeProvider>
    );
};

export default Dashboard;