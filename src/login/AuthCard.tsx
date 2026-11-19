import { Box, Paper } from "@mui/material";
import type React from "react";
import mdvLogo from "./assets/mdv_logo.png";

/**
 * The page around signing in: background, card and logo. Shared so that the
 * pages someone can land on while signing in all look like the same place.
 */
const AuthCard: React.FC<{ children: React.ReactNode }> = ({ children }) => (
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
            {children}
        </Paper>
    </Box>
);

export default AuthCard;
