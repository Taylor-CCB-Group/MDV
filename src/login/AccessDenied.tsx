import { Button, Stack, Typography } from "@mui/material";
import type React from "react";
import { getAppMountPath } from "@/utils/mdvRouting";
import AuthCard from "./AuthCard";
import type { AccessDenial, DenialAction } from "./accessDenial";

const buttonBase = {
    py: 1.5,
    px: 4,
    fontSize: "0.875rem",
    textTransform: "none",
    borderRadius: 5,
} as const;

// The card is light in both themes, so every colour is set here: the theme's
// dark-mode colours would make the text unreadable on this card.
const mainButton = {
    ...buttonBase,
    color: "#fff",
    bgcolor: "#1a73e8",
    "&:hover": { bgcolor: "#1557b0" },
};
const otherButton = {
    ...buttonBase,
    color: "#1a73e8",
    borderColor: "#dadce0",
    "&:hover": { borderColor: "#1a73e8", bgcolor: "rgba(26, 115, 232, 0.04)" },
};

/**
 * Shown on the login page when the server turned a signed-in person away.
 *
 * Deliberately has no plain "Sign in" button. Auth0 still holds a session, so
 * signing in again would complete at once and land back here - the loop this
 * replaces. A retry is offered only as the step after fixing something, and
 * switching account goes through /logout, which ends the Auth0 session first.
 */
const AccessDenied: React.FC<{ denial: AccessDenial; account: string | null }> = ({ denial, account }) => {
    const mountPath = getAppMountPath();
    const hrefFor = (action: DenialAction) => `${mountPath}${action.kind === "retry" ? "login" : "logout"}`;

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

            {account ? (
                <Typography sx={{ mb: 1.5, color: "#202124", fontSize: "16px", textAlign: "center" }}>
                    You signed in as <strong>{account}</strong>.
                </Typography>
            ) : null}

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

            <Stack spacing={1.5} sx={{ width: "100%", maxWidth: 280 }}>
                {denial.actions.map((action, index) => (
                    <Button
                        key={action.kind}
                        fullWidth
                        variant={index === 0 ? "contained" : "outlined"}
                        href={hrefFor(action)}
                        sx={index === 0 ? mainButton : otherButton}
                    >
                        {action.label}
                    </Button>
                ))}
            </Stack>
        </AuthCard>
    );
};

export default AccessDenied;
