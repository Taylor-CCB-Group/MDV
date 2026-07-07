import { Alert, Box, Paper, Stack, Typography } from "@mui/material";
import type { AdminSession } from "../api";

type ProfileViewProps = {
    session: AdminSession | null;
};

export function ProfileView({ session }: ProfileViewProps) {
    return (
        <Paper className="admin-showcase-panel admin-profile-panel" variant="outlined">
            <Stack spacing={2.5}>
                <Box>
                    <Typography variant="h5">Admin Profile</Typography>
                    {/* <Typography variant="body2">
                        Mock profile details for the current admin session.
                    </Typography> */}
                </Box>
                {session ? (
                    <Stack spacing={1.5}>
                        <Box className="admin-profile-row">
                            <Typography variant="body2">Email</Typography>
                            <Typography>{session.user.email}</Typography>
                        </Box>
                        <Box className="admin-profile-row">
                            <Typography variant="body2">Role</Typography>
                            <Typography>{session.isAdmin ? "Admin" : "No admin access"}</Typography>
                        </Box>
                        <Box className="admin-profile-row">
                            <Typography variant="body2">Authentication</Typography>
                            <Typography>{session.authEnabled ? "Auth0 enabled" : "Local dev mode"}</Typography>
                        </Box>
                        <Box className="admin-profile-row">
                            <Typography variant="body2">Assigned workspace</Typography>
                            <Typography>MDV deployment administrator</Typography>
                        </Box>
                        <Alert severity="info">
                            Profile data is mocked for now and can be connected to real account metadata later.
                        </Alert>
                    </Stack>
                ) : null}
            </Stack>
        </Paper>
    );
}
