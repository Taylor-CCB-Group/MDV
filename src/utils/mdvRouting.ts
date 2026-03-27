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

export function getDirParam() {
    return new URLSearchParams(window.location.search).get("dir");
}

export function isProjectDir(dir: string) {
    return /\/project\/[^/]+\/?$/.test(asUrl(dir).pathname);
}

export function isProjectPath(pathname = window.location.pathname) {
    return /^\/project\/[^/]+\/?$/.test(pathname);
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

export function getDashboardApiRoot() {
    const dir = getDirParam();
    if (!dir) return "/";
    return getApiRootFromDir(dir);
}

export function buildApiUrl(path: string, apiRoot = getDashboardApiRoot()) {
    const normalizedPath = path.replace(/^\/+/, "");

    if (apiRoot === "/") {
        return `/${normalizedPath}`;
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
    if (apiRoot === "/") return "/";
    return `/?dir=${encodeURIComponent(ensureTrailingSlash(apiRoot))}`;
}

export function buildProjectUrl(projectId: string | number, apiRoot = getDashboardApiRoot()) {
    if (apiRoot === "/") {
        return `/project/${projectId}`;
    }

    return `/?dir=${encodeURIComponent(buildApiUrl(`project/${projectId}`, apiRoot))}`;
}

export function shouldRenderDashboard() {
    const params = new URLSearchParams(window.location.search);
    if (params.get("popout") === "true") return false;
    if (isProjectPath()) return false;

    const dir = params.get("dir");
    return !(dir && isProjectDir(dir));
}
