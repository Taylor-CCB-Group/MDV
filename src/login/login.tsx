import type React from "react";
import { useEffect } from "react";
import { Box, Button, Paper, Typography } from "@mui/material";
import mdvLogo from "./assets/mdv_logo.png";

const Login: React.FC = () => {
    useEffect(() => {
        // Check if this is a callback from Auth0
        const params = new URLSearchParams(window.location.search);
        if (params.get('callback')) {
            const base = import.meta.env.DEV ? "http://localhost:5170?dir=/" : "";
            window.location.href = `${base}/catalog_dev`;
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
        <Box
            sx={{
                minHeight: "100vh",
                display: "flex",
                alignItems: "center",
                justifyContent: "center",
                bgcolor: "#f1f3f4",
            }}
        >
            <Paper
                elevation={1}
                sx={{
                    p: 6,
                    width: "100%",
                    maxWidth: 450,
                    display: "flex",
                    flexDirection: "column",
                    alignItems: "center",
                    borderRadius: 2,
                    border: "1px solid #dadce0",
                    backgroundColor: "unset",
                }}
            >
                <Box
                    component="img"
                    src={mdvLogo}
                    alt="Logo"
                    sx={{
                        width: 300,
                        height: "auto",
                        mb: 4,
                    }}
                />

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
                <Typography
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
                        }}
                        onClick={handleLoginWithSSO}
                    >
                        Sign in with SSO
                    </Button>
            </Paper>
        </Box>
    );
};

export default Login;
