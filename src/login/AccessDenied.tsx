import { getAppMountPath } from "@/utils/mdvRouting";
import { Button, Stack, Typography } from "@mui/material";
import type React from "react";
import AuthCard from "./AuthCard";
import type { AccessDenial } from "./accessDenial";

/**
 * Shown on the login page when the server turned a signed-in person away.
 *
 * Deliberately has no "Sign in" button. Auth0 still holds a session, so signing
 * in again would complete at once and land back here - the loop this replaces.
 * "Try again" does the same thing, but is offered as a retry after the person
 * has fixed something; "Sign out" ends the Auth0 session so another account can
 * be used.
 */
const AccessDenied: React.FC<{ denial: AccessDenial }> = ({ denial }) => {
    const mountPath = getAppMountPath();

    return (
        <AuthCard>
            <Typography
                variant="h5"
                sx={{
                    mb: 2,
                    color: "#202124",
                    fontWeight: 400,
                    fontSize: "24px",
                    textAlign: "center",
                }}
            >
                {denial.title}
            </Typography>

            <Typography
                sx={{
                    mb: 4,
                    color: "#202124",
                    fontSize: "16px",
                    textAlign: "center",
                }}
            >
                {denial.detail}
            </Typography>

            <Stack direction={{ xs: "column", sm: "row" }} spacing={1.5}>
                <Button
                    variant="contained"
                    href={`${mountPath}login`}
                    sx={{
                        py: 1.5,
                        px: 4,
                        fontSize: "0.875rem",
                        textTransform: "none",
                        borderRadius: 5,
                        bgcolor: "#1a73e8",
                        "&:hover": {
                            bgcolor: "#1557b0",
                        },
                        minWidth: 160,
                    }}
                >
                    Try again
                </Button>
                <Button
                    variant="outlined"
                    href={`${mountPath}logout`}
                    sx={{
                        py: 1.5,
                        px: 4,
                        fontSize: "0.875rem",
                        textTransform: "none",
                        borderRadius: 5,
                        minWidth: 160,
                    }}
                >
                    Sign out
                </Button>
            </Stack>
        </AuthCard>
    );
};

export default AccessDenied;
