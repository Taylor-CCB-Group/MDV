/**
 * Explains why a piece of the app's code failed to load.
 *
 * The browser's own message ("Failed to fetch dynamically imported module: <url>") says
 * neither what the server answered nor whether the file exists, and the browser cannot be
 * asked. So the named file is fetched again and the answer described: a proxy's error page,
 * a missing file, or a file that is there after all (the first request failed on the way).
 */

/** The file a failed import or CSS preload names, in Chrome's, Firefox's and Vite's wording. */
export function failedModuleUrl(error: unknown): string | undefined {
    const message = error instanceof Error ? error.message : typeof error === "string" ? error : "";
    const match =
        /(?:Failed to fetch dynamically imported module|error loading dynamically imported module|Unable to preload CSS for):?\s+(\S+)/i.exec(
            message,
        );
    return match?.[1];
}

/** What the server answers for the file now. */
export type ImportProbe = {
    url: string;
    status?: number;
    statusText?: string;
    contentType?: string;
    /** `<title>` of an HTML answer, such as a proxy's "502 Proxy Error" page. */
    title?: string;
    /** Why the request itself failed, when it got no answer at all. */
    fetchError?: string;
};

/** Just the part of `fetch` this needs, so tests can pass a plain function. */
export type FetchFn = (url: string, init?: RequestInit) => Promise<Response>;

export async function probeUrl(url: string, fetchFn: FetchFn = fetch): Promise<ImportProbe> {
    try {
        const response = await fetchFn(url, { cache: "no-store" });
        const contentType = response.headers.get("content-type")?.split(";")[0].trim() || undefined;
        let title: string | undefined;
        if (contentType === "text/html") {
            title = htmlTitle(await response.text());
        }
        return { url, status: response.status, statusText: response.statusText || undefined, contentType, title };
    } catch (error) {
        return { url, fetchError: error instanceof Error ? error.message : String(error) };
    }
}

/** `<title>` of an HTML page, such as a proxy's "502 Proxy Error". */
export function htmlTitle(text: string): string | undefined {
    return /<title>([^<]*)<\/title>/i.exec(text)?.[1].trim() || undefined;
}

/** "502 Proxy Error" rather than "502 (502 Proxy Error)" when the page title already has the status. */
export function describeAnswer(status: number, title?: string): string {
    if (!title) return `${status}`;
    return title.startsWith(`${status}`) ? title : `${status} (${title})`;
}

/** A sentence or two for someone who has just seen the failure. */
export function explainProbe(probe: ImportProbe): string {
    const { status, title, contentType, fetchError } = probe;
    if (fetchError !== undefined) {
        return `Fetching the file again failed too (${fetchError}). The server may have stopped, or the connection to it is down.`;
    }
    if (status === undefined) return "";
    const answer = describeAnswer(status, title);
    // Code and styles never come as HTML: a page here is a sign-in page, an app shell served
    // for unknown paths, or some other stand-in for the file.
    if (status >= 200 && status < 300 && contentType === "text/html") {
        return (
            `Fetching the file again got a web page${title ? ` ("${title}")` : ""} instead of the file. ` +
            "The server, or a proxy in front of it, answers this path with a page: the file may " +
            "not exist, or the session may need signing in again."
        );
    }
    if (status >= 200 && status < 300) {
        return (
            "The file is on the server: fetching it again worked. The first request failed on the " +
            "way, which happens when a proxy between the browser and MDV drops requests while " +
            "many files load at once. Reloading usually works."
        );
    }
    if (status === 502 || status === 503 || status === 504) {
        return (
            `Fetching the file again got ${answer} from a proxy between the browser and MDV, ` +
            "so the proxy could not reach MDV or MDV did not answer it."
        );
    }
    if (status === 404) {
        return (
            "The server does not have this file. The page may come from a different build of MDV " +
            "than the server is running now; reloading the page should fix it."
        );
    }
    return `Fetching the file again got ${answer}.`;
}

/** Fetches the file a load error names and explains the answer; undefined when it names none. */
export async function diagnoseImportFailure(
    error: unknown,
    fetchFn: FetchFn = fetch,
): Promise<{ probe: ImportProbe; explanation: string } | undefined> {
    const url = failedModuleUrl(error);
    if (!url) return undefined;
    const probe = await probeUrl(url, fetchFn);
    return { probe, explanation: explainProbe(probe) };
}
