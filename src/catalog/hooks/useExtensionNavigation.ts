import { apiFetch } from "@/utils/mdvRouting";
import { useEffect, useState } from "react";

export interface ExtensionNavigationItem {
    id: string;
    label: string;
    url: string;
}

function isRecord(value: unknown): value is Record<string, unknown> {
    return typeof value === "object" && value !== null;
}

function isExtensionNavigationItem(value: unknown): value is ExtensionNavigationItem {
    return (
        isRecord(value) &&
        typeof value.id === "string" &&
        typeof value.label === "string" &&
        typeof value.url === "string"
    );
}

export function parseExtensionNavigation(value: unknown): ExtensionNavigationItem[] {
    if (!isRecord(value) || !Array.isArray(value.extensions)) {
        return [];
    }
    return value.extensions.filter(isExtensionNavigationItem);
}

export default function useExtensionNavigation(): ExtensionNavigationItem[] {
    const [items, setItems] = useState<ExtensionNavigationItem[]>([]);

    useEffect(() => {
        let cancelled = false;

        async function fetchNavigation(): Promise<void> {
            try {
                const response = await apiFetch("extension_navigation");
                if (!response.ok) {
                    return;
                }
                const result: unknown = await response.json();
                if (!cancelled) {
                    setItems(parseExtensionNavigation(result));
                }
            } catch {
                // Navigation is optional; the catalog remains usable without it.
            }
        }

        void fetchNavigation();
        return () => {
            cancelled = true;
        };
    }, []);

    return items;
}
