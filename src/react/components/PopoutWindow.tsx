import { useCallback, useEffect, useRef, useState } from "react";
import { createPortal } from "react-dom";
import createCache, { type EmotionCache } from "@emotion/cache";
import { CacheProvider } from "@emotion/react";

export type PopoutWindowProps = {
    children: React.ReactNode;
    onClose: () => void;
    features?: string;
    title?: string;
};

export function PopoutWindow({ children, features = "width=600,height=400", title, onClose }: PopoutWindowProps) {
    const [isLoading, setIsLoading] = useState(true);
    // Container div
    const containerEl = useRef<HTMLElement>(document.createElement("div"));
    // Popout window
    const externalWindow = useRef<Window | null>(null);
    // Theme and style observers
    const themeObserver = useRef<MutationObserver | null>(null);
    const observer = useRef<MutationObserver | null>(null);
    const onCloseRef = useRef(onClose);
    const emotionCacheRef = useRef<EmotionCache | null>(null);

    useEffect(() => {
        // Update onClose ref when onClose is changed
        onCloseRef.current = onClose;
    }, [onClose]);

    const initializePopout = useCallback(async () => {
        // Open a new window
        externalWindow.current = window.open("", "_blank", features);
        if (!externalWindow.current) return;

        externalWindow.current.document.title = title || "Chat MDV";
        const bodyStyles = window.getComputedStyle(document.body);

        // Copying font styles
        externalWindow.current.document.body.style.fontFamily = bodyStyles.fontFamily;
        externalWindow.current.document.body.style.fontSize = bodyStyles.fontSize;
        externalWindow.current.document.body.style.lineHeight = bodyStyles.lineHeight;
        externalWindow.current.document.body.style.color = bodyStyles.color;
        externalWindow.current.document.body.style.margin = bodyStyles.margin;
        externalWindow.current.document.body.style.padding = bodyStyles.padding;

        // Create loading screen
        const loadingDiv = externalWindow.current.document.createElement("div");
        loadingDiv.style.cssText = `
            position: fixed; top: 0; left: 0; width: 100%; height: 100%;
            background: white; display: flex; align-items: center; justify-content: center;
            font-family: ${bodyStyles.fontFamily};
            font-size: 16px; color: #666; z-index: 9999;
        `;
        loadingDiv.innerHTML = `
            <div style="text-align: center;">
                <div style="width: 40px; height: 40px; border: 4px solid #f3f3f3; border-top: 4px solid #3498db; border-radius: 50%; animation: spin 1s linear infinite; margin: 0 auto 16px;"></div>
                <div>Loading...</div>
            </div>
            <style>@keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }</style>
        `;
        externalWindow.current.document.body.appendChild(loadingDiv);

        // Create Emotion cache for mui
        emotionCacheRef.current = createCache({
            key: "popout-mui",
            container: externalWindow.current.document.head,
        });

        // Sync theme
        const syncTheme = () => {
            if (!externalWindow.current) return;
            externalWindow.current.document.documentElement.className = document.documentElement.className;
        };
        syncTheme();

        themeObserver.current = new MutationObserver(syncTheme);
        themeObserver.current.observe(document.documentElement, {
            attributes: true,
            attributeFilter: ["class"],
        });

        // Copy all styles
        const styleElements = Array.from(document.querySelectorAll('link[rel="stylesheet"], style'));

        await Promise.all(
            styleElements.map((element) => {
                if (!externalWindow.current) return;
                if (element instanceof HTMLLinkElement) {
                    const newLink = externalWindow.current.document.createElement("link");
                    newLink.rel = "stylesheet";
                    newLink.href = element.href;
                    externalWindow.current.document.head.appendChild(newLink);
                } else if (element instanceof HTMLStyleElement) {
                    const newStyle = externalWindow.current.document.createElement("style");
                    newStyle.innerHTML = element.innerHTML;
                    externalWindow.current.document.head.appendChild(newStyle);
                }
            }),
        );

        // Style Observer
        observer.current = new MutationObserver((mutations) => {
            mutations.forEach((mutation) => {
                mutation.addedNodes.forEach(async (node) => {
                    if (!externalWindow.current) return;
                    if (node instanceof HTMLLinkElement) {
                        const newLink = externalWindow.current.document.createElement("link");
                        newLink.rel = "stylesheet";
                        newLink.href = node.href;
                        externalWindow.current.document.head.appendChild(newLink);
                    } else if (node instanceof HTMLStyleElement) {
                        const newStyle = externalWindow.current.document.createElement("style");
                        newStyle.innerHTML = node.innerHTML;
                        externalWindow.current.document.head.appendChild(newStyle);
                    }
                });
            });
        });

        observer.current.observe(document.head, { childList: true });

        // Append the container to window
        externalWindow.current.document.body.appendChild(containerEl.current);

        // Remove loading and show content
        loadingDiv.remove();
        setIsLoading(false);

        const handleUnload = () => {
            if (onCloseRef.current) onCloseRef.current();
        };
        externalWindow.current.addEventListener("beforeunload", handleUnload);
    }, [features, title]);

    useEffect(() => {
        initializePopout();

        // Cleanup
        return () => {
            if (observer.current) observer.current.disconnect();
            if (themeObserver.current) themeObserver.current.disconnect();
            if (externalWindow.current) {
                externalWindow.current.close();
            }
        };
    }, [initializePopout]);

    if (!containerEl.current || !emotionCacheRef.current || isLoading) {
        return null;
    }

    return createPortal(<CacheProvider value={emotionCacheRef.current}>{children}</CacheProvider>, containerEl.current);
}

export default PopoutWindow;
