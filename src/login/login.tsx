import type React from "react";
import { useEffect } from "react";
import { Button, Typography } from "@mui/material";
import { buildDashboardUrl } from "@/utils/mdvRouting";
import AuthCard from "./AuthCard";

const Login: React.FC = () => {
    useEffect(() => {
        // Check if this is a callback from Auth0
        const params = new URLSearchParams(window.location.search);
        if (params.get('callback')) {
            window.location.href = buildDashboardUrl();
        }
    }, []);

    const handleLogin = async () => {
        try {
            // Redirects to Auth0's login interface
            // The correct backend endpoint must be set here
            window.location.href = "login";
        } catch (err) {
            console.error("Login error:", err);
        }
    };

    // Needs to be changed later to standard routing
    const handleLoginWithSSO = async () => {
        try {
            window.location.href = "login_sso";
        } catch (err) {
            console.error("Login error:", err);
        }
    };

    return (
        <AuthCard>
            <Typography
                variant="h5"
                sx={{
                    mb: 1,
                    color: "#202124",
                    fontWeight: 400,
                    fontSize: "24px",
                }}
            >
                Sign in
            </Typography>

            <Typography
                sx={{
                    mb: 4,
                    color: "#202124",
                    fontSize: "16px",
                }}
            >
                to continue to MDV
            </Typography>

            <Button
                variant="contained"
                onClick={handleLogin}
                sx={{
                    py: 1.5,
                    px: 6,
                    fontSize: "0.875rem",
                    textTransform: "none",
                    borderRadius: 5,
                    bgcolor: "#1a73e8",
                    "&:hover": {
                        bgcolor: "#1557b0",
                    },
                    minWidth: 200,
                }}
            >
                Sign in
            </Button>
            {/* <Typography
                sx={{
                    mt: 6,
                    mb: 1,
                    color: "#202124",
                    fontSize: "16px",
                }}
            >
                If you are an Oxford University User:
            </Typography>

                <Button
                    variant="contained"
                    sx={{
                        mt: 1.5,
                        py: 1.5,
                        px: 6,
                        fontSize: "0.875rem",
                        textTransform: "none",
                        borderRadius: 5,
                        minWidth: 200,
                        backgroundColor: "white",
                        color: "black",
                        "&:hover": {
                            bgcolor: "grey",
                        },
                    }}
                    onClick={handleLoginWithSSO}
                >
                    Sign in with SSO
                </Button> */}
        </AuthCard>
    );
};

export default Login;
