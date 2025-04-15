import { createRoot } from "react-dom/client";
import { createTheme, ThemeProvider } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";
import { type PropsWithChildren, StrictMode, useMemo } from "react";
import { ProjectProvider } from "@/modules/ProjectContext";
import { observer } from "mobx-react-lite";
import type BaseChart from "@/charts/BaseChart";
import type { BaseDialog } from "@/utilities/Dialog";
import { OuterContainerProvider, useOuterContainer } from "./screen_state";
import { createFilterOptions } from "@mui/material";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
// import { ReactQueryDevtools } from '@tanstack/react-query-devtools';

// todo - think about whether this might lead to unexpected future issues
// consider virtualization etc
const filterOptions = createFilterOptions<any>({ limit: 256 });
/** This makes sure any material-ui components have an appropriate theme applied, including
 * setting `container` in `defaultProps` of any `<Popper>` and `<Popover>` components
 * which should ensure that tooltips / dropdowns etc work correctly without needing to manually
 * `useOuterContainer`.
 */
const MaterialWrapper = observer(function MaterialWrapper({
    children,
}: PropsWithChildren) {
    const prefersDarkMode = window.mdv?.chartManager?.theme === "dark";
    const container = useOuterContainer();
    const defaultProps = useMemo(() => ({ container }), [container]);

    const theme = useMemo(
        () =>
            createTheme({
                palette: {
                    mode: prefersDarkMode ? "dark" : "light",
                },
                components: {
                    MuiPopper: { defaultProps },
                    MuiPopover: { defaultProps },
                    MuiDialog: { defaultProps },
                    MuiAutocomplete: { defaultProps: { filterOptions } },
                    // There is a bit of a general issue with material-ui favouring less dense layouts
                    // than we would tend to want - it is aiming for a more mobile-friendly UX...
                    MuiTextField: { defaultProps: { size: "small" }},
                    MuiRadio: {
                        styleOverrides: {
                            root: {
                                padding: '5px',
                            },
                        }
                    }
                },
                typography: {
                    fontFamily: [
                        "Roboto",
                        '"Helvetica Neue"',
                        "Arial",
                        "sans-serif",
                    ].join(","),
                    fontSize: 12,
                }
            }),
        [prefersDarkMode, defaultProps],
    );

    return (
        <ThemeProvider theme={theme}>
            <CssBaseline enableColorScheme>{children}</CssBaseline>
        </ThemeProvider>
    );
});

// should we be keeping this local to each react root? - not sure if it matters either way.
// using coderabbitai suggestion for options - may review if we're using it for other things than histogram.
const queryClient = new QueryClient({
    defaultOptions: {
        queries: {
            refetchOnWindowFocus: false,
        }
    }
});

/**
 * It also makes sure that any common global context, styles etc are applied.
 * todo - this is a placeholder for refactoring so that there is a single react root,
 * then charts/dialogs etc will be rendered into it with portals.
 *
 * @param component - the component to render
 * @param container - the container to render into
 * @param parent - 'parent' chart or dialog that can be used to determine the container to be used for anything that needs to be rendered outside of the main container.
 * Currently this will be something that has `{ observable: { container: HTMLElement } }` as a property.
 * If not provided, the default is document.body.
 * @returns the root element that was created
 */
const createMdvPortal = (component: JSX.Element, container: HTMLElement, parent?: BaseChart<any> | BaseDialog) => {
    const root = createRoot(container);
    root.render(
        <StrictMode>
            <QueryClientProvider client={queryClient}>
                {/* <ReactQueryDevtools initialIsOpen={false} /> */}
                <OuterContainerProvider parent={parent}>
                    <MaterialWrapper>
                        <ProjectProvider>{component}</ProjectProvider>
                    </MaterialWrapper>
                </OuterContainerProvider>
            </QueryClientProvider>
        </StrictMode>,
    );
    return root;
};

export { createMdvPortal };
