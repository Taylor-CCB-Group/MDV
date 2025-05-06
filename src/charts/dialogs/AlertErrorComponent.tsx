import { Alert, AlertTitle, Paper, Typography } from "@mui/material";

export type AlertErrorComponentProps = {
    message: string;
    title?: string;
};

const AlertErrorComponent = ({ message, title }: AlertErrorComponentProps) => {
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
                    bgcolor: "var(--background_color_error)",
                    color: "var(--text_color)",
                    height: "100%",
                    width: "100%",
                }}
            >
                <>
                    <Alert
                        severity="error"
                        sx={{
                            "& .MuiAlert-message": {
                                width: "100%",
                            },
                            bgcolor: "transparent",
                            "& .MuiAlert-icon": {
                                color: "var(--icon_color_error)",
                            },
                            color: "var(--text_color_error)",
                            padding: 2,
                            display: "flex",
                            alignItems: title ? "flex-start" : "center",
                        }}
                    >
                        {title && (
                            <AlertTitle variant="h6">
                                {title}
                            </AlertTitle>
                        )}
                        <Typography 
                            variant="subtitle1"
                        >
                            {message}
                        </Typography>
                    </Alert>
                </>
            </Paper>
        </div>
    );
};

export default AlertErrorComponent;
