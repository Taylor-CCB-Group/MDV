function isAbsoluteUrl(value: string) {
    return /^[a-zA-Z][a-zA-Z\d+\-.]*:/.test(value);
}

function trimTrailingSlash(value: string) {
    if (value === "/") return value;
    return value.replace(/\/+$/, "");
}

function ensureTrailingSlash(value: string) {
    return value.endsWith("/") ? value : `${value}/`;
}

function normalizeSameOriginPath(value: string) {
    const trimmed = trimTrailingSlash(value);
    return trimmed || "/";
}

function asUrl(value: string) {
    return new URL(value, window.location.href);
}

function buildAppShellUrl(search = "") {
    const mountPath = getAppMountPath();
    return `${mountPath}${search}`;
}

function isLocalHost(hostname = window.location.hostname) {
    return hostname === "localhost" || hostname === "127.0.0.1";
}

export function isHostedPreviewHost(hostname = window.location.hostname) {
    return hostname.endsWith(".netlify.app");
}

export function getDefaultPreviewApiRoot() {
    return "http://localhost:5055";
}

export function isDefaultPreviewApiRoot(apiRoot: string) {
    return ensureTrailingSlash(apiRoot) === ensureTrailingSlash(getDefaultPreviewApiRoot());
}

export function getDirParam() {
    return new URLSearchParams(window.location.search).get("dir");
}

export function isProjectDir(dir: string) {
    return /\/project\/[^/]+\/?$/.test(asUrl(dir).pathname);
}

export function isProjectPath(pathname = window.location.pathname) {
    return /\/project\/[^/]+\/?$/.test(pathname);
}

export function getAppMountPath(pathname = window.location.pathname) {
    const normalizedPath = pathname || "/";

    if (isProjectPath(normalizedPath)) {
        const basePath = normalizedPath.replace(/\/project\/[^/]+\/?$/, "");
        return ensureTrailingSlash(basePath || "/");
    }

    for (const suffix of ["/login_dev", "/login_dev.html", "/catalog_dev", "/catalog_dev.html", "/index.html"]) {
        if (normalizedPath.endsWith(suffix)) {
            const basePath = normalizedPath.slice(0, -suffix.length);
            return ensureTrailingSlash(basePath || "/");
        }
    }

    if (normalizedPath.endsWith("/")) {
        return normalizedPath;
    }

    return ensureTrailingSlash(normalizedPath);
}

export function getApiRootFromDir(dir: string) {
    const url = asUrl(dir);
    const projectMarker = url.pathname.lastIndexOf("/project/");
    const apiPath =
        projectMarker >= 0 ? url.pathname.slice(0, projectMarker) : url.pathname;

    if (url.origin === window.location.origin) {
        return normalizeSameOriginPath(apiPath);
    }

    return trimTrailingSlash(`${url.origin}${apiPath}`) || url.origin;
}

export function normalizeApiRoot(apiRoot: string) {
    if (!apiRoot || apiRoot === "/") return "/";

    if (isAbsoluteUrl(apiRoot)) {
        return ensureTrailingSlash(trimTrailingSlash(apiRoot));
    }

    const normalized = apiRoot.startsWith("/") ? apiRoot : `/${apiRoot}`;
    return ensureTrailingSlash(trimTrailingSlash(normalized));
}

export function getDashboardApiRoot() {
    const dir = getDirParam();
    if (!dir) {
        if (isHostedPreviewHost() && !isLocalHost()) {
            return getDefaultPreviewApiRoot();
        }
        return "/";
    }
    return getApiRootFromDir(dir);
}

export function getProjectDirFromLocation() {
    const dir = getDirParam();
    if (dir) return dir;

    if (isHostedPreviewHost() && isProjectPath() && !isLocalHost()) {
        return `${getDefaultPreviewApiRoot()}${window.location.pathname}`;
    }

    return `${window.location.origin}${window.location.pathname}`;
}

export function buildApiUrl(path: string, apiRoot = getDashboardApiRoot()) {
    const normalizedPath = path.replace(/^\/+/, "");

    if (apiRoot === "/") {
        return `${getAppMountPath()}${normalizedPath}`;
    }

    if (isAbsoluteUrl(apiRoot)) {
        return new URL(normalizedPath, ensureTrailingSlash(apiRoot)).toString();
    }

    return `${trimTrailingSlash(apiRoot)}/${normalizedPath}`;
}

export function apiFetch(input: string, init?: RequestInit) {
    return fetch(buildApiUrl(input), init);
}

export function buildDashboardUrl(apiRoot = getDashboardApiRoot()) {
    if (!getDirParam() && isHostedPreviewHost() && isDefaultPreviewApiRoot(apiRoot)) {
        return "/";
    }
    if (apiRoot === "/") return getAppMountPath();
    return buildAppShellUrl(`?dir=${encodeURIComponent(ensureTrailingSlash(apiRoot))}`);
}

export function buildProjectUrl(projectId: string | number, apiRoot = getDashboardApiRoot()) {
    if (!getDirParam() && isHostedPreviewHost() && isDefaultPreviewApiRoot(apiRoot)) {
        return `/project/${projectId}`;
    }
    if (apiRoot === "/") {
        return `${getAppMountPath()}project/${projectId}`;
    }

    return buildAppShellUrl(`?dir=${encodeURIComponent(buildApiUrl(`project/${projectId}`, apiRoot))}`);
}

export function shouldShowLocalBackendNotice() {
    return isHostedPreviewHost() && !getDirParam();
}

export function shouldRenderDashboard() {
    const params = new URLSearchParams(window.location.search);
    if (params.get("popout") === "true") return false;
    if (isProjectPath()) return false;

    const dir = params.get("dir");
    return !(dir && isProjectDir(dir));
}
