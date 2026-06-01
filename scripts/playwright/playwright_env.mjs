const CURSOR_SANDBOX_BROWSERS_MARKER = "cursor-sandbox-cache";

/**
 * Env for spawning Playwright. Cursor agent shells set PLAYWRIGHT_BROWSERS_PATH to a
 * sandbox cache under /var/folders/... that usually has no browser binaries; unset it
 * so Playwright uses the default install (e.g. ~/Library/Caches/ms-playwright).
 */
export function playwrightSpawnEnv(baseEnv = process.env) {
    const env = { ...baseEnv };
    const browsersPath = env.PLAYWRIGHT_BROWSERS_PATH;
    if (browsersPath?.includes(CURSOR_SANDBOX_BROWSERS_MARKER)) {
        delete env.PLAYWRIGHT_BROWSERS_PATH;
    }
    return env;
}
