import { useState } from "react";
import { useColorMode } from "@/ThemeProvider";
import {
    AppBar,
    Box,
    Button,
    Chip,
    IconButton,
    Stack,
    Toolbar,
    Tooltip,
} from "@mui/material";
import Brightness4Icon from "@mui/icons-material/Brightness4";
import Brightness7Icon from "@mui/icons-material/Brightness7";
import { ArrowLeft, LogOut, RefreshCw } from "lucide-react";
import { buildApiUrl } from "@/utils/mdvRouting";
import type { AdminSession } from "../api";
import { adminRoutes, type AdminAccessState, type AdminRoute } from "../adminUtils";

const adminToolbarButtonSx = {
    borderWidth: 1.5,
    borderRadius: 2,
    textTransform: "none",
    px: 1.5,
    py: 0.75,
} as const;

type AdminHeaderProps = {
    activeRoute: AdminRoute;
    accessState: AdminAccessState;
    authEnabled: boolean;
    session: AdminSession | null;
    syncingUsers: boolean;
    onNavigate: (route: AdminRoute) => void;
    onGoToProjects: () => void;
    onSignIn: () => void;
    onSignOut: () => void;
    onSyncUsers: () => void;
};

export function AdminHeader({
    activeRoute,
    accessState,
    authEnabled,
    session,
    syncingUsers,
    onNavigate,
    onGoToProjects,
    onSignIn,
    onSignOut,
    onSyncUsers,
}: AdminHeaderProps) {
    const [customLogoVisible, setCustomLogoVisible] = useState(true);
    const { mode, toggleColorMode } = useColorMode();

    return (
        <AppBar position="static" color="default" elevation={0} className="admin-app-bar">
            <Toolbar className="admin-toolbar">
                <Box className="admin-header-brand">
                    <Button
                        variant="text"
                        onClick={onGoToProjects}
                        className="admin-brand-button"
                    >
                        MDV
                    </Button>
                    <Box
                        component="img"
                        className="admin-secondary-logo"
                        alt="Custom deployment logo"
                        src={buildApiUrl("secondary_logo")}
                        onError={() => setCustomLogoVisible(false)}
                        onClick={onGoToProjects}
                        sx={{ display: customLogoVisible ? "block" : "none" }}
                    />
                </Box>
                <Stack
                    direction="row"
                    spacing={1.5}
                    flexWrap="wrap"
                    useFlexGap
                    className="admin-header-nav"
                >
                    {adminRoutes.map((route) => {
                        const Icon = route.icon;
                        const selected = activeRoute === route.id;
                        return (
                            <Button
                                key={route.id}
                                size="small"
                                variant={selected ? "contained" : "outlined"}
                                color={selected ? "primary" : "inherit"}
                                className="admin-nav-button"
                                startIcon={<Icon size={16} strokeWidth={2.1} />}
                                onClick={() => onNavigate(route.id)}
                                sx={selected ? { textTransform: "none" } : adminToolbarButtonSx}
                            >
                                {route.label}
                            </Button>
                        );
                    })}
                </Stack>
                <Box className="admin-header-actions">
                    <Tooltip
                        title={
                            authEnabled
                                ? "Sync Auth0 users into the MDV database"
                                : "User sync is only available when Auth0 is enabled"
                        }
                    >
                        <span>
                            <Button
                                variant="outlined"
                                startIcon={<RefreshCw size={18} />}
                                disabled={!authEnabled || syncingUsers}
                                onClick={onSyncUsers}
                                sx={adminToolbarButtonSx}
                            >
                                {syncingUsers ? "Syncing..." : "Sync users"}
                            </Button>
                        </span>
                    </Tooltip>
                    <Button
                        variant="outlined"
                        startIcon={<ArrowLeft size={18} />}
                        onClick={onGoToProjects}
                        sx={adminToolbarButtonSx}
                    >
                        Back to catalog
                    </Button>
                    {accessState === "login_required" ? (
                        <Button
                            variant="outlined"
                            onClick={onSignIn}
                            sx={adminToolbarButtonSx}
                        >
                            Sign in
                        </Button>
                    ) : null}
                    {authEnabled && session ? (
                        <Button
                            variant="outlined"
                            startIcon={<LogOut size={18} />}
                            onClick={onSignOut}
                            sx={adminToolbarButtonSx}
                        >
                            Sign Out
                        </Button>
                    ) : null}
                    <Tooltip title="Toggle theme">
                        <IconButton
                            onClick={toggleColorMode}
                            color="inherit"
                            aria-label="Toggle theme"
                        >
                            {mode === "dark" ? <Brightness4Icon /> : <Brightness7Icon />}
                        </IconButton>
                    </Tooltip>
                </Box>
            </Toolbar>
        </AppBar>
    );
}
