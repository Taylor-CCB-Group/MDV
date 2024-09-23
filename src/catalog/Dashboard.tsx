import type React from "react";
import { useState, useEffect, useCallback, useMemo } from "react";
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
} from "@mui/material";
import { alpha } from "@mui/material/styles";
import ProjectCard from "./ProjectCard";
interface Project {
    id: string;
    name: string;
    type: "Editable" | "Read-Only";
    lastModified: string;
}

const theme = createTheme({
    palette: {
        primary: {
            main: "#2c3e50", // A deep, muted blue
        },
        secondary: {
            main: "#34495e", // A slightly lighter shade for contrast
        },
        background: {
            default: "#ecf0f1", // A light gray background
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

type SortOption = "lastModified" | "name";

const Dashboard: React.FC = () => {
    const [projects, setProjects] = useState<Project[]>([]);
    const [filter, setFilter] = useState("");
    const [viewMode, setViewMode] = useState<"grid" | "list">("grid");
    const [isModalOpen, setIsModalOpen] = useState(false);
    const [newProjectName, setNewProjectName] = useState("");
    const [projectType, setProjectType] = useState<"Editable" | "Read-Only">(
        "Editable",
    );
    const [sortBy, setSortBy] = useState<SortOption>("lastModified");
    const [sortOrder, setSortOrder] = useState<"asc" | "desc">("desc");
    const [isDropdownOpen, setIsDropdownOpen] = useState(false);
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);

    useEffect(() => {
        // Mocked project data
        setProjects([
            {
                id: "1",
                name: "Project A",
                type: "Editable",
                lastModified: "Apr 1, 2024",
            },
            {
                id: "2",
                name: "Project B",
                type: "Read-Only",
                lastModified: "Apr 2, 2024",
            },
            {
                id: "3",
                name: "Project C",
                type: "Editable",
                lastModified: "Apr 3, 2024",
            },
        ]);
    }, []);

    const filteredAndSortedProjects = useMemo(() => {
        const result = projects.filter((p) =>
            p.name.toLowerCase().includes(filter.toLowerCase()),
        );

        result.sort((a, b) => {
            if (sortBy === "name") {
                return sortOrder === "asc"
                    ? a.name.localeCompare(b.name)
                    : b.name.localeCompare(a.name);
            }
            return sortOrder === "asc"
                ? new Date(a.lastModified).getTime() -
                      new Date(b.lastModified).getTime()
                : new Date(b.lastModified).getTime() -
                      new Date(a.lastModified).getTime();
        });

        return result;
    }, [projects, filter, sortBy, sortOrder]);

    const handleCreateProject = useCallback(() => {
        if (newProjectName.trim()) {
            const newProject = {
                id: Date.now().toString(),
                name: newProjectName.trim(),
                type: projectType,
                lastModified: new Date().toLocaleDateString("en-US", {
                    month: "short",
                    day: "numeric",
                    year: "numeric",
                }),
            };
            setProjects((prevProjects) => [...prevProjects, newProject]);
            setNewProjectName("");
            setProjectType("Editable");
            setIsModalOpen(false);
        }
    }, [newProjectName, projectType]);

    const handleDeleteProject = useCallback((id: string) => {
        setProjects((prevProjects) => prevProjects.filter((p) => p.id !== id));
    }, []);

    const handleSort = (option: SortOption) => {
        if (sortBy === option) {
            setSortOrder(sortOrder === "asc" ? "desc" : "asc");
        } else {
            setSortBy(option);
            setSortOrder("desc");
        }
        setIsDropdownOpen(false); // Close the dropdown after selection
    };

    const toggleDropdown = () => {
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
                                value={filter}
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
                        {/* Create new project card */}
                        <Grid item xs={12} sm={6} md={3}>
                            <ButtonBase
                                onClick={() => {
                                    setIsModalOpen(true); // Handle the create new project action
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

                        {/* Create from template card */}
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

                        {/* Create from template card 2*/}
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
                            elevation={1} // Reduced elevation for less conspicuousness
                            sx={{
                                padding: "8px", // Small padding for internal spacing
                                display: "flex",
                                alignItems: "center",
                                justifyContent: "center",
                                borderRadius: "4px", // Minimal border radius for a cleaner look
                                width: "205px", // Fixed width to ensure consistency
                                height: "50px", // Fixed height for uniformity
                                bgcolor: "background.paper", // Use background color to blend with the UI
                            }}
                        >
                            <Button
                                endIcon={<ExpandMore />}
                                onClick={(event) => {
                                    setAnchorEl(event.currentTarget);
                                    setIsDropdownOpen(true);
                                }}
                                sx={{
                                    textTransform: "none",
                                    width: "100%", // Ensures the button takes full width
                                    height: "100%", // Ensures the button takes full height
                                    display: "flex",
                                    justifyContent: "center", // Center the text
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
                                Last modified{" "}
                                {sortBy === "lastModified" &&
                                    (sortOrder === "asc" ? "↑" : "↓")}
                            </MenuItem>
                            <MenuItem onClick={() => handleSort("name")}>
                                Name{" "}
                                {sortBy === "name" &&
                                    (sortOrder === "asc" ? "↑" : "↓")}
                            </MenuItem>
                        </Menu>
                    </Box>

                    <Divider sx={{ mb: 2 }} />

                    <Grid container spacing={4}>
                        {filteredAndSortedProjects.map((project) => (
                            <Grid
                                item
                                key={project.id}
                                xs={12}
                                sm={6}
                                md={4}
                                lg={3}
                            >
                                <ProjectCard
                                    createdAt={"Jan 1, 2024"}
                                    owner={"Ben"}
                                    collaborators={[
                                        "lorem@gmail.com",
                                        "ipsum@gmail.com",
                                        "dolor@gmail.com",
                                    ]}
                                    numberOfStructures={"32"}
                                    numberOfImages={"64"}
                                    {...project}
                                    onDelete={handleDeleteProject}
                                />
                            </Grid>
                        ))}
                    </Grid>
                </Container>
            </Box>
        </ThemeProvider>
    );
};

export default Dashboard;
