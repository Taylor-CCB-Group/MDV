import { Alert, AlertTitle, Paper, Typography } from "@mui/material";

export type AlertType = "error" | "info" | "success" | "warning";

export type AlertErrorComponentProps = {
    message: string;
    title?: string;
    alertType?: AlertType;
};

const AlertErrorComponent = ({ message, title, alertType = "error" }: AlertErrorComponentProps) => {
    const backgroundColor = `var(--background_color_${alertType})`;
    const textColor = `var(--text_color_${alertType})`;

    return (
        <div
            style={{
                margin: "20px auto",
                maxWidth: "1000px",
                minHeight: "15vh",
                display: "flex",
                justifyContent: "center",
                alignItems: "center",
            }}
        >
            <Paper
                elevation={3}
                sx={{
                    bgcolor: backgroundColor,
                    color: "var(--text_color)",
                    height: "100%",
                    width: "100%",
                }}
            >
                <>
                    <Alert
                        severity={alertType}
                        sx={{
                            "& .MuiAlert-message": {
                                width: "100%",
                            },
                            bgcolor: "transparent",
                            color: textColor,
                            padding: 2,
                            display: "flex",
                            alignItems: title ? "flex-start" : "center",
                        }}
                    >
                        {title && <AlertTitle variant="h6">{title}</AlertTitle>}
                        <Typography variant="subtitle1">{message}</Typography>
                    </Alert>
                </>
            </Paper>
        </div>
    );
};

export default AlertErrorComponent;
