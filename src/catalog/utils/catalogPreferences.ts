import type { SortBy, SortOrder } from "./projectUtils";

export type ViewMode = "grid" | "list";
export type ListColumn = "owner" | "lastModified" | "type" | "filesystem";

/** Optional list columns, in display order. Name and Actions are always shown. */
export const LIST_COLUMNS: { id: ListColumn; label: string; hiddenOnPublicPage?: boolean }[] = [
    { id: "owner", label: "Owner" },
    { id: "lastModified", label: "Last Modified", hiddenOnPublicPage: true },
    { id: "type", label: "Type" },
    { id: "filesystem", label: "Filesystem", hiddenOnPublicPage: true },
];

export type CatalogPreferences = {
    viewMode: ViewMode;
    sortBy: SortBy;
    sortOrder: SortOrder;
    columns: ListColumn[];
};

export const DEFAULT_CATALOG_PREFERENCES: CatalogPreferences = {
    viewMode: "grid",
    sortBy: "lastModified",
    sortOrder: "desc",
    columns: ["lastModified"],
};

export const CATALOG_PREFERENCES_KEY = "mdv.catalog.preferences";

const isViewMode = (v: unknown): v is ViewMode => v === "grid" || v === "list";
const isSortBy = (v: unknown): v is SortBy => v === "name" || v === "lastModified";
const isSortOrder = (v: unknown): v is SortOrder => v === "asc" || v === "desc";
const isColumn = (v: unknown): v is ListColumn => LIST_COLUMNS.some((c) => c.id === v);

/** Reads stored preferences, falling back to defaults per field if missing or invalid. */
export function loadCatalogPreferences(): CatalogPreferences {
    const defaults = DEFAULT_CATALOG_PREFERENCES;
    try {
        const raw = localStorage.getItem(CATALOG_PREFERENCES_KEY);
        if (!raw) return defaults;
        const parsed: unknown = JSON.parse(raw);
        if (typeof parsed !== "object" || parsed === null) return defaults;
        const p: Record<string, unknown> = { ...parsed };
        return {
            viewMode: isViewMode(p.viewMode) ? p.viewMode : defaults.viewMode,
            sortBy: isSortBy(p.sortBy) ? p.sortBy : defaults.sortBy,
            sortOrder: isSortOrder(p.sortOrder) ? p.sortOrder : defaults.sortOrder,
            columns: Array.isArray(p.columns) ? p.columns.filter(isColumn) : defaults.columns,
        };
    } catch {
        return defaults;
    }
}

export function saveCatalogPreferences(prefs: CatalogPreferences): void {
    try {
        localStorage.setItem(CATALOG_PREFERENCES_KEY, JSON.stringify(prefs));
    } catch {
        // storage unavailable or full: preferences just won't persist
    }
}

/** Whether the user has ever stored a preference (used so public-page defaults don't override it). */
export function hasStoredCatalogPreferences(): boolean {
    try {
        return localStorage.getItem(CATALOG_PREFERENCES_KEY) !== null;
    } catch {
        return false;
    }
}
