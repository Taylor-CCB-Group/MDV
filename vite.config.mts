import { execSync } from "node:child_process";
import * as fs from "node:fs";
import * as path from "node:path";
import { fileURLToPath } from "node:url";
import rolldownBabel from "@rolldown/plugin-babel";
import { babel as rollupBabel } from "@rollup/plugin-babel";
// import vitePluginSocketIO from 'vite-plugin-socket.io';
import react, { reactCompilerPreset } from "@vitejs/plugin-react";
import { visualizer } from "rollup-plugin-visualizer";
import { type ProxyOptions, type UserConfig, defineConfig } from "vite";
import glsl from "vite-plugin-glsl";

const configDir = path.dirname(fileURLToPath(import.meta.url));

/**
 * Pick the file an `exports` entry points at, preferring the ESM condition. Only the
 * shapes these packages actually publish: a string, or an object of conditions.
 */
function exportsTarget(entry: unknown): string | undefined {
    if (typeof entry === "string") return entry;
    if (typeof entry !== "object" || entry === null) return undefined;
    const conditions = entry as Record<string, unknown>;
    for (const condition of ["import", "default"]) {
        const value = conditions[condition];
        if (typeof value === "string") return value;
    }
    return undefined;
}

const escapeRegExp = (value: string) => value.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");

/**
 * Optional local SpatialData.js checkout (see `pnpm link:spatialdata`).
 *
 * Returns the packages that were found, each with the alias entries that point it at the
 * checkout: one per `exports` subpath, then the bare package name.
 *
 * Subpaths need their own entries and need to come first. Vite matches a *string* alias
 * only on an exact hit or a `/` boundary, so `@spatialdata/core/parquet-worker?worker&url`
 * — the query is part of the id — misses a `@spatialdata/core/parquet-worker` string alias
 * and falls through to the bare `@spatialdata/core → <pkgRoot>` one, which rewrites it to
 * `<pkgRoot>/parquet-worker`. Nothing is there: the real file is
 * `<pkgRoot>/dist/parquet-worker.js`. Hence a regex anchored to end-or-`?`, which both
 * matches and carries the query through to the replacement.
 *
 * Without this, every subpath import breaks under `pnpm link:spatialdata` —
 * `zarrextra/workers`, `@spatialdata/core/parquet-worker`, `@spatialdata/core/parquet-wasm`
 * — while a registry install is fine, so it only ever bites while testing an upstream
 * change. Reading `exports` rather than listing subpaths keeps it right when one is added.
 */
function spatialdataLinkAliases(): Array<{ name: string; alias: Array<{ find: RegExp; replacement: string }> }> {
    const root = process.env.SPATIALDATA_ROOT?.trim();
    if (!root) return [];
    const abs = path.resolve(root.startsWith("~/") ? path.join(process.env.HOME ?? "", root.slice(2)) : root);
    const packages: Array<[string, string]> = [
        ["@spatialdata/avivatorish", "packages/avivatorish"],
        ["@spatialdata/core", "packages/core"],
        ["@spatialdata/layers", "packages/layers"],
        ["@spatialdata/react", "packages/react"],
        ["@spatialdata/vis", "packages/vis"],
        ["zarrextra", "packages/zarrextra"],
    ];
    const linked: Array<{ name: string; alias: Array<{ find: RegExp; replacement: string }> }> = [];
    for (const [name, rel] of packages) {
        const pkgRoot = path.join(abs, rel);
        const manifestPath = path.join(pkgRoot, "package.json");
        if (!fs.existsSync(manifestPath)) continue;
        const alias: Array<{ find: RegExp; replacement: string }> = [];
        const exportsMap: unknown = JSON.parse(fs.readFileSync(manifestPath, "utf8")).exports;
        if (typeof exportsMap === "object" && exportsMap !== null) {
            for (const [key, entry] of Object.entries(exportsMap as Record<string, unknown>)) {
                // Wildcards would need a pattern alias, and none of these packages use one.
                if (!key.startsWith("./") || key.includes("*")) continue;
                const target = exportsTarget(entry);
                if (target) {
                    alias.push({
                        find: new RegExp(`^${escapeRegExp(`${name}/${key.slice(2)}`)}(?=$|\\?)`),
                        // `$` is a backreference in the replacement string; paths rarely
                        // contain one, but a checkout under `~/code/$work` would silently
                        // resolve to nonsense.
                        replacement: path.join(pkgRoot, target).replace(/\$/g, "$$$$"),
                    });
                }
            }
        }
        alias.push({ find: new RegExp(`^${escapeRegExp(name)}(?=$|/|\\?)`), replacement: pkgRoot });
        linked.push({ name, alias });
    }
    return linked;
}

const linkedSpatialdata = spatialdataLinkAliases();
/** Every alias entry, subpaths ahead of their bare package name — order is load-bearing. */
const spatialdataAliasEntries = linkedSpatialdata.flatMap((pkg) => pkg.alias);
const linkedSpatialdataNames = linkedSpatialdata.map((pkg) => pkg.name);

/**
 * Checkout roots Vite must be allowed to READ from when @spatialdata/* is a local
 * `link:`. The packages resolve through a symlink out of this project, so their
 * non-JS assets — core's vendored `parquet_wasm_bg.wasm`, zarrextra's
 * `codec-worker.js` — are outside the default `server.fs.allow` root and come back
 * 404, which leaves a spatial chart stuck on "Loading" with no error of its own.
 *
 * Derived from the resolved symlink rather than SPATIALDATA_ROOT, so it follows
 * whatever `pnpm link:spatialdata` actually wired up. Empty for a registry install.
 */
function linkedSpatialdataRoots(): string[] {
    const roots = new Set<string>();
    for (const name of ["@spatialdata/vis", "@spatialdata/core", "zarrextra"]) {
        try {
            const real = fs.realpathSync(path.join(configDir, "node_modules", name));
            // Path-aware, not `startsWith`: a sibling checkout at `<configDir>-other`
            // has the config dir as a string prefix, and would be read as "inside" —
            // silently dropping the one root that needed allowing.
            const rel = path.relative(configDir, real);
            const inside = rel === "" || (!rel.startsWith("..") && !path.isAbsolute(rel));
            // <checkout>/packages/<pkg> → <checkout>
            if (!inside) roots.add(path.resolve(real, "../.."));
        } catch {
            // not installed / not linked — nothing to allow
        }
    }
    return [...roots];
}

/**
 * `worker_format` comes from the shell (see the `build-flask-vite-jbrowse` script).
 * A typo used to be cast straight through, so Vite emitted workers in whatever it
 * made of the value and the failure surfaced far from the cause.
 */
function workerFormat(): "es" | "iife" {
    const requested = process.env.worker_format;
    if (requested === undefined || requested === "") return "iife";
    if (requested !== "es" && requested !== "iife") {
        throw new Error(`worker_format must be "es" or "iife", got "${requested}"`);
    }
    return requested;
}

const spatialdataFsAllow = linkedSpatialdataRoots();

// zarrita needed a polyfill for Buffer - seems like a bug
// seems ok without as long we don't use ZipFileStore (marked experimental anyway)
// having the polyfill means the build works, but devserver fails with 'cannot import outside a module'
// (not in the code using zarrita, but in unrelated worker modules)
// import { nodePolyfills } from 'vite-plugin-node-polyfills'

const flaskURL = "http://127.0.0.1:5055";
const port = Number(process.env.PORT || process.env.VITE_PORT || 5170);
const build = (process.env.build || "desktop_pt") as "production" | "dev_pt" | "desktop" | "desktop_pt";
// setting output path: use --outDir
// todo review --assetsDir / nofont / cleanup & consolidate entrypoints
// maybe also the various build configurations at some point.

const reactCompiler = reactCompilerPreset();
reactCompiler.rolldown.filter ??= {};
// this doesn't like decorators it seems, for now if we only point it at tsx files then it doesn't have a problem
// (as far as we've noticed)
reactCompiler.rolldown.filter.id = /\.tsx(?:$|\?)/;
const useReactCompiler = process.env.VITE_USE_REACT_COMPILER !== "0" && process.env.VITE_USE_REACT_COMPILER !== "false";
const enableBundleAnalysis = process.env.VITE_BUNDLE_ANALYZE === "1" || process.env.VITE_BUNDLE_ANALYZE === "true";

/** Same rules as main's per-build assetFileNames, plus fonts under assets/ (Rolldown emits url(./font) next to assets/mdv.css). */
function flaskAssetFileNames(assetInfo: { name?: string }): string {
    const name = assetInfo.name ?? "";
    // project_bootstrap / desktop_index import ./all_css → emitted as all_css.css
    if (name.includes("index.css") || name === "all_css.css" || name === "mdv.css" || name === "desktop_index.css") {
        return "assets/mdv.css";
    }
    if (name === "catalog.css") return "assets/catalog.css";
    if (process.env.VITE_ENTRYPOINT) {
        const { name: entryBase } = path.parse(process.env.VITE_ENTRYPOINT);
        if (name === `${entryBase}.css`) return "assets/mdv.css";
    }
    const ext = path.extname(name).slice(1).toLowerCase();
    // Wasm is emitted by the main graph AND by a worker graph — @spatialdata/core's
    // 6.6MB parquet-wasm goes in both. Worker builds do not take this function, so an
    // unhashed name here yields `img/parquet_wasm_bg.wasm` for one graph and the worker
    // default `assets/parquet_wasm_bg-<hash>.wasm` for the other: two copies of identical
    // bytes. Matching the worker default collapses them onto one file.
    if (ext === "wasm") return "assets/[name]-[hash][extname]";
    if (["woff", "woff2", "ttf", "eot"].includes(ext)) return "assets/[name][extname]";
    if (ext === "svg" && /^fa-(brands|regular|solid)-/.test(name)) return "assets/[name][extname]";
    return "img/[name][extname]";
}

/**
 * shim for various different build configurations.
 *
 * nb `dev_pt` which returns `{}` suffices for many things - vite devserver, netlify build.
 * `desktop_pt` just has a little more config setting a .ts input & making sure other things will work with flask template.
 * other methods are supposed to be for replacing other webpack configs.
 */
function getRollupOptions() {
    if (build === "production") {
        const version = process.env.mdv_version ? "-" + process.env.mdv_version : "";

        // somewhat equivalent to original webpack production build - not the current 'production' with new features.
        return {
            input: process.env.nofont ? "src/modules/basic_index_nf.js" : "src/modules/basic_index.js",
            output: {
                entryFileNames: `mdv${version}.js`,
                assetFileNames: flaskAssetFileNames,
            },
        };
    }
    if (build === "desktop_pt") {
        // currently there are different versions of entrypoint, this is the one with react etc.
        // used for Flask...
        return {
            input: {
                mdv: "src/modules/project_bootstrap.tsx",
                catalog: "src/catalog/catalog_index.tsx",
                login: "src/login/login_index.tsx",
            },
            output: {
                entryFileNames: "js/[name].js",
                assetFileNames: flaskAssetFileNames,
            },
        };
    }
    if (build === "dev_pt") {
        // version of vite build used for netlify deploy preview & devserver, using default 'index.html' entrypoint
        // which now picks dashboard vs project viewer at runtime.
        return {};
    }
    if (build === "desktop") {
        return {
            input: "src/modules/desktop_index.js",
            output: {
                entryFileNames: "js/mdv.js",
                assetFileNames: flaskAssetFileNames,
            },
        };
    }
    if (process.env.VITE_ENTRYPOINT) {
        // If you want a custom entrypoint - in particular, in order to have a custom DataLoader for interfacing
        // with another backend, you can specify it with VITE_ENTRYPOINT environment variable, e.g.
        // `VITE_ENTRYPOINT=path/to/my_index.js pnpm exec vite build --outDir path/to/output`
        // (nb, we may change the logic in this config...)
        return {
            input: process.env.VITE_ENTRYPOINT,
            output: {
                entryFileNames: "js/mdv.js",
                assetFileNames: flaskAssetFileNames,
            },
        };
    }
    throw new Error(`Unknown build type '${build}' and no VITE_ENTRYPOINT specified.`);
}

// avoiding some repition by defining a proxyOptions object used for all proxied routes.
// Flask already adds CORS + Range expose headers via add_safe_headers; keep changeOrigin so
// cross-origin Range requests (parquet / OME-TIFF / zarr via spatialdata.js) work through the proxy.
const proxyOptions = { target: flaskURL, changeOrigin: true };
// ... and then this is a bit more concise than
const proxy = [
    "^/(get_|images|tracks|save|chat|spatial).*", // single-project Flask routes (incl. /spatial zarr etc.)
    "^/project/[^/]+/.+", // proxy nested project routes, but keep /project/:id for the Vite app shell
    "^/.*\\.(json|b|gz)$",
    "/projects",
    "/create_project",
    "/import_project",
    "/export_project",
    "/delete_project",
    "/extension_config",
    "/ucsc_proxy",
    "/enable_auth",
    "/api_root",
    "/rescan_projects",
    "/login_dev",
    "/secondary_logo",
    // biome-ignore lint/performance/noAccumulatingSpread: don't care about performance in vite config
].reduce((acc, route) => ({ ...acc, [route]: proxyOptions }), {}) as Record<string, ProxyOptions>;
// (failed) attempt to let this proxy without cors_allowed_origins wildcard on server
// using more specific socketio in python for now
proxy["/socket.io"] = {
    target: flaskURL.replace("http:", "ws:"),
    changeOrigin: true,
    ws: true,
    // rewriteWsOrigin: true, // copilot says this is not needed in the current version of socket.io
    // ^^ not helping either way...
};

export default defineConfig(async (): Promise<UserConfig> => {
    // For local development, try to get Git info. This is guarded by a check for the .git directory
    // to prevent errors in environments where git is not available (like during Docker build, where
    // even in dev environment, the build happens before .git is copied in, for cache purposes).
    if (fs.existsSync(".git")) {
        try {
            const commitDate = execSync("git log -1 --format=%cI").toString().trimEnd();
            const branchName = execSync("git rev-parse --abbrev-ref HEAD").toString().trimEnd();
            const commitHash = execSync("git rev-parse HEAD").toString().trimEnd();
            const lastCommitMessage = execSync("git show -s --format=%s").toString().trimEnd();

            process.env.VITE_GIT_COMMIT_DATE = commitDate;
            process.env.VITE_GIT_BRANCH_NAME = branchName;
            process.env.VITE_GIT_COMMIT_HASH = commitHash;
            process.env.VITE_GIT_LAST_COMMIT_MESSAGE = lastCommitMessage;
            process.env.VITE_GIT_DIRTY = execSync('git diff --quiet || echo "dirty"').toString().trimEnd();
        } catch (e) {
            console.error("Failed to get git info:", e);
        }
    }
    process.env.VITE_BUILD_DATE = new Date().toISOString();

    return {
        base: process.env.asset_base || (build === "dev_pt" ? "/" : "./"),
        server: {
            headers: {
                // Match python/mdvtools/server_utils.add_safe_headers (SharedArrayBuffer + Range CORS).
                "Cross-Origin-Embedder-Policy": "require-corp",
                "Cross-Origin-Opener-Policy": "same-origin",
                "Cross-Origin-Resource-Policy": "cross-origin",
                "Access-Control-Allow-Origin": "*",
                "Access-Control-Allow-Methods": "GET",
                "Access-Control-Allow-Headers": "Content-Type, Range, X-Requested-With, Authorization",
                "Access-Control-Expose-Headers": "Content-Range, Content-Length, Accept-Ranges",
                "Access-Control-Max-Age": "86400",
            },
            port,
            strictPort: true,
            proxy,
            ...(spatialdataFsAllow.length ? { fs: { allow: [configDir, ...spatialdataFsAllow] } } : {}),
        },
        publicDir: process.env.exclude_dir ? false : "examples", //used for netlify.toml??... the rest is noise.
        build: {
            sourcemap: !process.env.nomap,
            rollupOptions: {
                ...getRollupOptions(),
                external: ["./python/**"],
            },
        },
        plugins: [
            glsl(),
            rollupBabel({
                babelHelpers: "bundled",
                extensions: [".js", ".jsx", ".ts", ".tsx", ".mjs", ".cjs"],
                include: ["src/**/*"],
                plugins: [
                    ["@babel/plugin-transform-typescript", { allowDeclareFields: true }],
                    "@babel/plugin-transform-class-static-block",
                    ["@babel/plugin-proposal-decorators", { version: "2023-05" }],
                    ["@babel/plugin-transform-class-properties", { loose: true }],
                ],
            }),
            react({
                include: [/\.tsx?$/, /\.jsx?$/],
            }),
            ...(useReactCompiler
                ? [
                      rolldownBabel({
                          presets: [reactCompiler],
                      }),
                  ]
                : []),
            ...(enableBundleAnalysis
                ? [
                      visualizer({
                          filename: process.env.VITE_BUNDLE_ANALYZE_OUTPUT || "dist/bundle-analysis.html",
                          template: "treemap",
                          gzipSize: true,
                          brotliSize: true,
                      }),
                  ]
                : []),
        ],
        worker: {
            format: workerFormat(),
        },
        resolve: {
            // Array form, not an object: the linked-checkout entries are regexes (see
            // `spatialdataLinkAliases`), and an object alias map can only key on strings.
            alias: [{ find: /^@(?=\/)/, replacement: path.resolve(configDir, "./src") }, ...spatialdataAliasEntries],
            // A linked checkout resolves its own bare imports from ITS node_modules, so the
            // renderer ends up with two of everything even at identical versions: two
            // @deck.gl/core (the shader hooks a layer declares are not the ones the assembler
            // knows — "DECKGL_FILTER_COLOR: no matching overloaded function"), two Reacts
            // (invalid hook call), two Matrix4. Collapse the packages that cross the
            // boundary onto MDV's copy. Only when linked, so a registry install is untouched.
            ...(spatialdataFsAllow.length
                ? {
                      dedupe: [
                          "react",
                          "react-dom",
                          "deck.gl",
                          "@deck.gl/core",
                          "@deck.gl/layers",
                          "@deck.gl/extensions",
                          "@deck.gl/geo-layers",
                          "@deck.gl/mesh-layers",
                          "@luma.gl/core",
                          "@luma.gl/engine",
                          "@luma.gl/shadertools",
                          "@luma.gl/webgl",
                          "@luma.gl/constants",
                          "@math.gl/core",
                          "@hms-dbmi/viv",
                          "@vivjs/views",
                          "@vivjs/constants",
                          "zarrita",
                      ],
                  }
                : {}),
        },
        // Vite 7 crawls all **/*.html for dep pre-bundling. Python/Flask templates and
        // public/index.html reference static/js/mdv.js (built output), which is not a
        // resolvable package — limit scanning to real Vite entry HTML files.
        optimizeDeps: {
            entries: [
                path.resolve(configDir, "index.html"),
                path.resolve(configDir, "src/static.html"),
                path.resolve(configDir, "src/obvios.html"),
                path.resolve(configDir, "login_dev.html"),
                path.resolve(configDir, "catalog_dev.html"),
            ],
            // @spatialdata/core requires Zod 4 while MDV uses Zod 3; excluding both preserves
            // each package's own dependency resolution instead of sharing the optimized
            // Zod 3 cache.
            // When SPATIALDATA_ROOT is set, exclude the whole linked set so Vite does
            // not freeze a stale prebundle of local dist builds.
            //
            // Do NOT extend this to every @spatialdata package for debuggability's
            // sake, tempting as it is. They pull CJS transitives — arrow's `long`, via
            // core's vendored parquet — that only work in the browser because
            // prebundling converts them to ESM. Excluding the packages leaves those
            // undiscovered and the app dies at load with "does not provide an export
            // named 'default'". Debuggability is paid for upstream instead: the
            // packages ship `dist/index.js.map` from 0.6.0 on, and Rolldown emits a
            // prebundle map without being asked, so a crash inside one names the
            // function rather than reading `Le (…/.vite/deps/@spatialdata_layers.js)`.
            // (An `optimizeDeps.esbuildOptions.sourcemap` lived here for that; under
            // Vite 8 it is deprecated, it warns, and removing it changes no output.)
            exclude: ["@spatialdata/core", "zod", ...linkedSpatialdataNames],
        },
    } as UserConfig;
});
