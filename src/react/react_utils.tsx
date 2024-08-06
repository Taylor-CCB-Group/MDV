import { createRoot } from "react-dom/client";
import useMediaQuery from '@mui/material/useMediaQuery';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import { type PropsWithChildren, StrictMode, useMemo } from 'react';
import { ProjectProvider } from "@/modules/ProjectContext";
import { observer } from "mobx-react-lite";

/** This makes sure any material-ui components have an appropriate theme applied.
 * Not clear that we want to use keep using material-ui, but for now it seems not to mix too badly.
 */
const MaterialWrapper = observer(function MaterialWrapper({ children }: PropsWithChildren) {
    // todo make this work with user override
    // const prefersDarkMode = useMediaQuery('(prefers-color-scheme: dark)');
    const prefersDarkMode = window.mdv.chartManager.theme === 'dark';

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
});

/**
 * It also makes sure that any common global context, styles etc are applied.
 * todo - this is a placeholder for refactoring so that there is a single react root,
 * then charts/dialogs etc will be rendered into it with portals.
 */
const createMdvPortal = (component: JSX.Element, container: HTMLElement) => {
    const root = createRoot(container);
    root.render((
        <StrictMode>
            <MaterialWrapper>
                <ProjectProvider>
                    {component}
                </ProjectProvider>
            </MaterialWrapper>
        </StrictMode>
    ));
    return root;
}


export { createMdvPortal }