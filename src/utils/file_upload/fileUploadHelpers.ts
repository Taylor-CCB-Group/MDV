export function formatFileSizeMb(size: number): string {
    return (size / (1024 * 1024)).toFixed(2);
}

export function getProjectIdFromRoot(root: string): string | null {
    const projectId = root.split("/").filter(Boolean).pop();
    return projectId ?? null;
}

export function getSocketPath(mainApiRoute: string): string {
    return `${mainApiRoute.replace(/\/$/, "")}/socket.io`;
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
