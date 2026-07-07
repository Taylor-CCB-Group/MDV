import { Alert, Button, Paper, Stack, Typography } from "@mui/material";
import { LogOut } from "lucide-react";
import type { AdminAccessState } from "../adminUtils";

type AccessStatePanelProps = {
    accessState: Exclude<AdminAccessState, null>;
    onSignIn: () => void;
    onSignOut: () => void;
};

export function AccessStatePanel({ accessState, onSignIn, onSignOut }: AccessStatePanelProps) {
    return (
        <Paper className="admin-showcase-panel" variant="outlined">
            <Stack spacing={1.5}>
                <Typography variant="h6">
                    {accessState === "login_required"
                        ? "Sign in required"
                        : "Admin access required"}
                </Typography>
                <Typography color="text.secondary">
                    {accessState === "login_required"
                        ? "You need to sign in before using the MDV Admin Portal."
                        : "Your account is signed in, but it is not marked as an MDV admin for this deployment."}
                </Typography>
                <Alert severity={accessState === "login_required" ? "info" : "warning"}>
                    {accessState === "login_required"
                        ? "After signing in, reload this page to continue."
                        : "Ask an existing admin or deployment owner to grant admin access before continuing."}
                </Alert>
                <Stack direction="row" spacing={1}>
                    {accessState === "login_required" ? (
                        <Button
                            variant="contained"
                            size="small"
                            onClick={onSignIn}
                        >
                            Sign in
                        </Button>
                    ) : (
                        <Button
                            variant="outlined"
                            size="small"
                            startIcon={<LogOut size={16} strokeWidth={2.1} />}
                            onClick={onSignOut}
                        >
                            Sign out
                        </Button>
                    )}
                </Stack>
            </Stack>
        </Paper>
    );
}
