import type React from "react";
import { createContext, useContext, useEffect, useMemo, useState } from "react";
import {
    ThemeProvider as MuiThemeProvider,
    createTheme,
    useMediaQuery,
    type PaletteMode,
    CssBaseline,
} from "@mui/material";

type ColorMode = {
    toggleColorMode: () => void;
    mode: PaletteMode;
};

const ColorModeContext = createContext<ColorMode>({
    toggleColorMode: () => {},
    mode: "light",
});

const getDesignTokens = (mode: PaletteMode) =>
    createTheme({
        palette: {
            mode,
            ...(mode === "light"
                ? {
                      // Light mode
                      primary: { main: "#2c3e50" },
                      secondary: { main: "#34495e" },
                      background: { default: "#ecf0f1", paper: "#ffffff" },
                      text: { primary: "#2c3e50", secondary: "#7f8c8d" },
                  }
                : {
                      // Dark mode
                      primary: { main: "#ffffff" },
                      secondary: { main: "#222" },
                      background: { default: "#353536", paper: "#404040" },
                      text: { primary: "#ffffff", secondary: "#ffffff" },
                  }),
        },
        typography: {
            fontFamily: '"Roboto", "Helvetica", "Arial", sans-serif',
            h5: { fontWeight: 500 },
            body1: { fontSize: "0.9rem" },
        },
        components: {
            MuiButton: {
                styleOverrides: {
                    root: { textTransform: "none" },
                },
            },
            MuiPaper: {
                styleOverrides: {
                    root: {
                        boxShadow:
                            "0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24)",
                    },
                },
            },
            MuiDialog: {
                styleOverrides: {
                    paper: {
                        backgroundImage: 'none',
                    },
                },
            },
        },
    });

export const CustomThemeProvider: React.FC<{ children: React.ReactNode }> = ({
    children,
}) => {
    const prefersDarkMode = useMediaQuery("(prefers-color-scheme: dark)");
    const [mode, setMode] = useState<PaletteMode>(
        prefersDarkMode ? "dark" : "light",
    );

    useEffect(() => {
        // Automatically update theme based on system preference
        setMode(prefersDarkMode ? "dark" : "light");
    }, [prefersDarkMode]);

    useEffect(() => {
        const htmlElement = document.querySelector("html");
        if (htmlElement) {
            htmlElement.classList.remove("light", "dark");
            htmlElement.classList.add(mode);
        }
    }, [mode]);

    const colorMode = useMemo(
        () => ({
            toggleColorMode: () => {
                setMode((prevMode) =>
                    prevMode === "light" ? "dark" : "light",
                );
            },
            mode,
        }),
        [mode],
    );

    const theme = useMemo(() => getDesignTokens(mode), [mode]);

    return (
        <ColorModeContext.Provider value={colorMode}>
            <MuiThemeProvider theme={theme}>
                <CssBaseline />
                {children}
            </MuiThemeProvider>
        </ColorModeContext.Provider>
    );
};

export const useColorMode = () => useContext(ColorModeContext);
