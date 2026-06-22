import { useEffect, useState } from "react";
import {
    Alert,
    Box,
    Chip,
    CircularProgress,
    Container,
    Paper,
    Stack,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow,
    Typography,
} from "@mui/material";
import {
    adminApi,
    type AdminProject,
    type AdminSession,
    type AdminUser,
} from "./api";

function formatUserName(user: AdminUser) {
    const name = [user.firstName, user.lastName].filter(Boolean).join(" ");
    return name || user.email;
}

export default function AdminApp() {
    const [session, setSession] = useState<AdminSession | null>(null);
    const [users, setUsers] = useState<AdminUser[]>([]);
    const [projects, setProjects] = useState<AdminProject[]>([]);
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        let cancelled = false;

        async function loadInitialState() {
            setLoading(true);
            setError(null);
            try {
                const [sessionResult, usersResult, projectsResult] = await Promise.all([
                    adminApi.session(),
                    adminApi.users(),
                    adminApi.projects(),
                ]);
                if (cancelled) return;
                setSession(sessionResult);
                setUsers(usersResult.users);
                setProjects(projectsResult.projects);
            } catch (err) {
                if (!cancelled) {
                    setError(err instanceof Error ? err.message : "Failed to load admin data");
                }
            } finally {
                if (!cancelled) {
                    setLoading(false);
                }
            }
        }

        void loadInitialState();
        return () => {
            cancelled = true;
        };
    }, []);

    return (
        <Box sx={{ minHeight: "100vh", bgcolor: "background.default" }}>
            <Box
                component="header"
                sx={{
                    bgcolor: "background.paper",
                    borderBottom: 1,
                    borderColor: "divider",
                    py: 2,
                }}
            >
                <Container maxWidth="xl">
                    <Stack
                        direction={{ xs: "column", md: "row" }}
                        spacing={2}
                        justifyContent="space-between"
                        alignItems={{ xs: "flex-start", md: "center" }}
                    >
                        <Box>
                            <Typography variant="h5">MDV Admin</Typography>
                            <Typography variant="body2" color="text.secondary">
                                Read-only deployment overview
                            </Typography>
                        </Box>
                        {session ? (
                            <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
                                <Chip
                                    color={session.isAdmin ? "success" : "default"}
                                    label={session.isAdmin ? "Admin session" : "No admin access"}
                                    size="small"
                                />
                                <Chip
                                    label={session.authEnabled ? "Auth enabled" : "Local dev auth"}
                                    size="small"
                                />
                                <Chip label={session.user.email} size="small" />
                            </Stack>
                        ) : null}
                    </Stack>
                </Container>
            </Box>

            <Container maxWidth="xl" sx={{ py: 3 }}>
                <Stack spacing={2}>
                    {error ? <Alert severity="error">{error}</Alert> : null}
                    {loading ? (
                        <Stack alignItems="center" sx={{ py: 8 }}>
                            <CircularProgress />
                        </Stack>
                    ) : (
                        <Stack
                            direction={{ xs: "column", lg: "row" }}
                            spacing={2}
                            alignItems="flex-start"
                        >
                            <Paper sx={{ p: 2, flex: 1, width: "100%" }}>
                                <Stack spacing={2}>
                                    <Stack
                                        direction="row"
                                        justifyContent="space-between"
                                        alignItems="center"
                                    >
                                        <Typography variant="h6">Projects</Typography>
                                        <Chip label={projects.length} size="small" />
                                    </Stack>
                                    {projects.length === 0 ? (
                                        <Alert severity="info">
                                            No projects found in the MDV database.
                                        </Alert>
                                    ) : (
                                        <Table size="small">
                                            <TableHead>
                                                <TableRow>
                                                    <TableCell>Name</TableCell>
                                                    <TableCell>ID</TableCell>
                                                    <TableCell>Access</TableCell>
                                                    <TableCell>Status</TableCell>
                                                </TableRow>
                                            </TableHead>
                                            <TableBody>
                                                {projects.map((project) => (
                                                    <TableRow key={project.id}>
                                                        <TableCell>{project.name}</TableCell>
                                                        <TableCell>{project.id}</TableCell>
                                                        <TableCell>
                                                            {project.accessLevel || "unknown"}
                                                        </TableCell>
                                                        <TableCell>
                                                            {project.isDeleted
                                                                ? "Deleted"
                                                                : "Active"}
                                                        </TableCell>
                                                    </TableRow>
                                                ))}
                                            </TableBody>
                                        </Table>
                                    )}
                                </Stack>
                            </Paper>

                            <Paper sx={{ p: 2, flex: 1, width: "100%" }}>
                                <Stack spacing={2}>
                                    <Stack
                                        direction="row"
                                        justifyContent="space-between"
                                        alignItems="center"
                                    >
                                        <Typography variant="h6">Users</Typography>
                                        <Chip label={users.length} size="small" />
                                    </Stack>
                                    {users.length === 0 ? (
                                        <Alert severity="info">
                                            No users found in the MDV database yet.
                                        </Alert>
                                    ) : (
                                        <Table size="small">
                                            <TableHead>
                                                <TableRow>
                                                    <TableCell>User</TableCell>
                                                    <TableCell>Email</TableCell>
                                                    <TableCell>Status</TableCell>
                                                    <TableCell>Admin</TableCell>
                                                </TableRow>
                                            </TableHead>
                                            <TableBody>
                                                {users.map((user) => (
                                                    <TableRow key={user.id}>
                                                        <TableCell>{formatUserName(user)}</TableCell>
                                                        <TableCell>{user.email}</TableCell>
                                                        <TableCell>
                                                            {user.isActive ? "Active" : "Inactive"}
                                                        </TableCell>
                                                        <TableCell>
                                                            {user.isAdmin ? "Yes" : "No"}
                                                        </TableCell>
                                                    </TableRow>
                                                ))}
                                            </TableBody>
                                        </Table>
                                    )}
                                </Stack>
                            </Paper>
                        </Stack>
                    )}
                </Stack>
            </Container>
        </Box>
    );
}
