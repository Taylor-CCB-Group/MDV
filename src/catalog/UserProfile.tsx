import { Person } from "@mui/icons-material";
import {
    Avatar,
    Box,
    CircularProgress,
    Divider,
    IconButton,
    Menu,
    MenuItem,
    Typography,
} from "@mui/material";
import type React from "react";
import { useState } from "react";
import useUser from "./hooks/useUser";

const UserProfile: React.FC = () => {
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
    const { user, isLoading, error } = useUser();
    const open = Boolean(anchorEl);

    const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
        setAnchorEl(event.currentTarget);
    };

    const handleClose = () => {
        setAnchorEl(null);
    };

    const handleSignOut = () => {
        try {
            window.location.href = "logout";
        } catch (err) {
            console.error("Logout error:", err);
        }
    };

    return (
        <>
            <IconButton
                onClick={handleClick}
                size="small"
                sx={{ ml: 1, mr: 1.5 }}
                aria-controls={open ? "account-menu" : undefined}
                aria-haspopup="true"
                aria-expanded={open ? "true" : undefined}
                aria-label="Account menu"
            >
                {isLoading ? (
                    <CircularProgress size={32} />
                ) : (
                    <Avatar
                        sx={{ width: 32, height: 32, color: "inherit" }}
                        src={user?.avatarUrl}
                        alt={user?.name ?? "Account"}
                    >
                        <Person />
                    </Avatar>
                )}
            </IconButton>
            <Menu
                anchorEl={anchorEl}
                id="account-menu"
                open={open}
                onClose={handleClose}
                sx={{
                    overflow: "visible",
                    filter: "drop-shadow(0px 2px 8px rgba(0,0,0,0.32))",
                    mt: 1.5,
                    "& .MuiAvatar-root": {
                        width: 60,
                        height: 60,
                        ml: -1,
                        mr: 1,
                    },
                }}
                transformOrigin={{ horizontal: "right", vertical: "top" }}
                anchorOrigin={{ horizontal: "right", vertical: "bottom" }}
            >
                {isLoading && (
                    <Box
                        sx={{
                            p: 2,
                            display: "flex",
                            justifyContent: "center",
                            minWidth: 200,
                        }}
                    >
                        <CircularProgress size={24} />
                    </Box>
                )}
                {error && (
                    <>
                    <Box sx={{ p: 2, maxWidth: 280 }}>
                        <Typography color="error" variant="body2">
                            {error}
                        </Typography>
                    </Box>
                        <Divider />
                        </>
                )}
                {user && (
                    <Box sx={{ p: 2 }}>
                        <Box
                            sx={{
                                display: "flex",
                                alignItems: "flex-start",
                                mb: 2,
                                gap: 1,
                                px: 1,
                            }}
                        >
                            <Avatar
                                sx={{ width: 100, height: 100, color: "inherit" }}
                                src={user.avatarUrl}
                                alt={user.name}
                            >
                                {!user.avatarUrl && (
                                    <Person sx={{ fontSize: 50 }} />
                                )}
                            </Avatar>
                            <Box>
                                <Typography
                                    variant="subtitle1"
                                    sx={{ fontWeight: 600 }}
                                >
                                    {user.name}
                                </Typography>
                                <Typography
                                    variant="body2"
                                    color="text.secondary"
                                    sx={{ mt: 0.5 }}
                                >
                                    {user.email}
                                </Typography>
                                <Typography
                                    variant="body2"
                                    color="text.secondary"
                                    sx={{ mt: 1 }}
                                >
                                    {user.association}
                                </Typography>
                            </Box>
                        </Box>
                        <Divider />
                    </Box>
                )}
                <MenuItem onClick={handleSignOut}>Sign out</MenuItem>
            </Menu>
        </>
    );
};

export default UserProfile;
