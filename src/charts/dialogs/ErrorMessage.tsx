import { Alert, AlertTitle, Paper, Typography } from "@mui/material";
import "../css/charts.css";

export type ErrorMessageProps = {
    message: string;
    title?: string;
};

const ErrorMessage = ({ message, title = "Error Occurred" }: ErrorMessageProps) => {
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
                            padding: 2,
                            display: "flex",
                            alignItems: "center",
                        }}
                    >
                        <Typography 
                            variant="subtitle1"
                            color="error"
                        >
                            {message}
                        </Typography>
                    </Alert>
                </>
            </Paper>
        </div>
    );
};

export default ErrorMessage;
