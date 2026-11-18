import { buildApiUrl } from "@/utils/mdvRouting";

export function formatFileSizeMb(size: number): string {
    return (size / (1024 * 1024)).toFixed(2);
}

export function getProjectIdFromRoot(root: string): string | null {
    const projectId = root.split("/").filter(Boolean).pop();
    return projectId ?? null;
}

/** Socket.IO `path` option for the API root. buildApiUrl resolves an api root of "/"
 * against the page's mount path, which a prefix-stripping proxy hides from the server. */
export function getSocketPath(mainApiRoute: string): string {
    return new URL(buildApiUrl("socket.io", mainApiRoute), window.location.href).pathname;
}

export function buildViewRedirectUrl(root: string, viewName: string): string {
    const projectPath = root.startsWith("http") ? new URL(root).pathname : root;
    const { origin } = window.location;

    if (import.meta.env.DEV) {
        const params = new URLSearchParams();
        params.set("dir", projectPath);
        params.set("view", viewName);
        return `${origin}/?${params.toString()}`;
    }

    return `${origin}${projectPath}?view=${encodeURIComponent(viewName)}`;
}
