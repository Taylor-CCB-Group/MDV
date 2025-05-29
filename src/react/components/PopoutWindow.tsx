import { useEffect, useRef } from "react";
import { createPortal } from "react-dom";

export type PopoutWindowProps = {
    children: React.ReactNode;
    onClose: () => void;
    features?: string;
    title?: string;
}

export function PopoutWindow({ children, features = "width=600,height=400", title, onClose }: PopoutWindowProps) {
    const containerEl = useRef<HTMLElement>(document.createElement("div"));
    const externalWindow = useRef<Window | null>(null);
    // had to use these refs as the style changes were not reflecting properly
    const themeObserver = useRef<MutationObserver>();
    const observer = useRef<MutationObserver>();
    const onCloseRef = useRef(onClose);

    // Keep onClose ref updated
    useEffect(() => {
        onCloseRef.current = onClose;
    }, [onClose]);

    useEffect(() => {
        // open a new window
        externalWindow.current = window.open("", "_blank", features);
        if (!externalWindow.current) return;

        // assign a title to it
        externalWindow.current.document.title = title || "Chat Window";

        // Function to synchronize the theme by syncing the class names
        const syncTheme = () => {
            if (!externalWindow.current) return;
            externalWindow.current.document.documentElement.className = document.documentElement.className;
        };

        // Initial sync
        syncTheme();

        // Set up a MutationObserver on the main window
        themeObserver.current = new MutationObserver(syncTheme);
        themeObserver.current.observe(document.documentElement, {
            attributes: true,
            attributeFilter: ["class"],
        });

        // Function to add a stylesheet or style element to the popout window
        const addStyleElement = (styleElement: HTMLStyleElement | HTMLLinkElement) => {
            if (!externalWindow.current) {
                return;
            }
            if (styleElement instanceof HTMLLinkElement) {
                const newLink = externalWindow.current.document.createElement("link");
                newLink.rel = "stylesheet";
                newLink.href = styleElement.href;
                externalWindow.current.document.head.appendChild(newLink);
            } else if (styleElement instanceof HTMLStyleElement) {
                const newStyle = externalWindow.current.document.createElement("style");
                newStyle.innerHTML = styleElement.innerHTML;
                externalWindow.current.document.head.appendChild(newStyle);
            }
        };

        // Initial copy of all existing styles
        Array.from(document.querySelectorAll('link[rel="stylesheet"], style')).forEach((styleElement) => {
            addStyleElement(styleElement as HTMLStyleElement | HTMLLinkElement);
        });

        // Observe the main window's head for new stylesheets and style elements
        observer.current = new MutationObserver((mutations) => {
            mutations.forEach((mutation) => {
                if (mutation.type === "childList" && mutation.addedNodes.length > 0) {
                    mutation.addedNodes.forEach((node) => {
                        if (node instanceof HTMLLinkElement || node instanceof HTMLStyleElement) {
                            addStyleElement(node);
                        }
                    });
                }
            });
        });

        // Start observing the main window's head for changes
        observer.current.observe(document.head, {
            childList: true,
            subtree: true,
        });

        // append the container to the window
        externalWindow.current.document.body.appendChild(containerEl.current);

        // call onClose before unload
        const handleUnload = () => {
            if (onCloseRef.current) onCloseRef.current();
        };

        externalWindow.current.addEventListener("beforeunload", handleUnload);

        // Cleanup
        return () => {
            if (observer.current) observer.current.disconnect();
            if (themeObserver.current) themeObserver.current.disconnect();
            if (externalWindow.current) {
                externalWindow.current.removeEventListener("beforeunload", handleUnload);
                externalWindow.current.close();
            }
        };
    }, [features, title]);

    if (!containerEl.current) {
        return null;
    }

    // create a new portal
    return createPortal(children, containerEl.current);
}

export default PopoutWindow;
