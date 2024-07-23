import useMediaQuery from '@mui/material/useMediaQuery';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import { useMemo } from 'react';

export function MaterialWrapper({children}) {
    // todo make this work with user override
    const prefersDarkMode = useMediaQuery('(prefers-color-scheme: dark)');

    const theme = useMemo(() => createTheme({
        palette: {
            mode: prefersDarkMode ? 'dark' : 'light',
        },
    }), [prefersDarkMode]);

    return (
        <ThemeProvider theme={theme}>
            <CssBaseline enableColorScheme>
                {children}
            </CssBaseline>
        </ThemeProvider>
    )
}