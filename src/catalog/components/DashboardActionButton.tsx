import { Button, Paper, Typography } from "@mui/material";
import type { ReactNode } from "react";

interface DashboardActionButtonProps {
    icon: ReactNode;
    label: string;
    onClick: () => void;
    active?: boolean;
}

const DashboardActionButton = ({
    icon,
    label,
    onClick,
    active = false,
}: DashboardActionButtonProps) => (
    <Paper
        elevation={1}
        sx={{
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            borderRadius: "4px",
            bgcolor: active ? "primary.main" : "background.paper",
            color: active ? "primary.contrastText" : "text.primary",
            marginRight: 2,
            height: "50px",
        }}
    >
        <Button
            sx={{
                color: "inherit",
                padding: 1,
                display: "flex",
                justifyContent: "space-around",
                height: "100%",
                width: "100%",
                textTransform: "none",
            }}
            onClick={onClick}
        >
            {icon}
            <Typography sx={{ marginLeft: 1 }}>
                {label}
            </Typography>
        </Button>
    </Paper>
);

export default DashboardActionButton;
