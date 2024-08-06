import { createRoot } from "react-dom/client";
import useMediaQuery from '@mui/material/useMediaQuery';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import { useMemo } from 'react';

/** This makes sure any material-ui components have an appropriate theme applied.
 * Not clear that we want to use keep using material-ui, but for now it seems not to mix too badly.
 */
function MaterialWrapper({ children }) {
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

/**
 * todo - this is a placeholder for refactoring so that there is a single react root,
 * then charts/dialogs etc will be rendered into it with portals.
 */

const createMdvPortal = (component: JSX.Element, container: HTMLElement) => {
    const root = createRoot(container);
    root.render((
        <MaterialWrapper>
            {component}
        </MaterialWrapper>
    ));
    return root;
}


export { createMdvPortal }