/**
 * Guards the class of bug where a bundled chunk keeps a literal runtime path to a file
 * the build never emitted, so the app 404s in production while dev works fine (the dev
 * server serves the dependency's own tree from node_modules).
 *
 * The instance that prompted it was @spatialdata/core's parquet-wasm loader, which
 * reached the glue by a `@vite-ignore`d relative path no consumer bundler could resolve
 * (MDV#539 copied core's `vendor/` tree in to satisfy it). Core 0.9.0 loads it through
 * the `@spatialdata/core/parquet-wasm` export instead, so the bundler emits the wasm
 * itself and that reference is gone.
 *
 * Deliberately not an assertion about parquet-wasm specifically — any chunk pointing at
 * anything unemitted fails. It does mean the test passes vacuously if nothing escapes
 * its own directory, which is the healthy state.
 */
import * as fs from "node:fs";
import * as path from "node:path";
import { describe, expect, test } from "vitest";

const repoRoot = path.resolve(__dirname, "../../..");

/** Every outDir the build scripts in package.json write to. */
const BUILD_OUT_DIRS = [
    "vite-dist",
    "dist",
    "dist/flask",
    "dist/mdv",
    "python/mdvtools/static",
    "python/dist_hotfix/flask",
];

/**
 * Rolldown chunk output only — not the build root, `examples/` (publicDir) or the copied
 * vendor tree. Excluding the root is not just about noise: a `../` reference from there
 * escapes the build tree by definition, so `isServedFromBuild` rejects it and scanning the
 * root cannot flag anything. Root entries (`production` emits `mdv-<version>.js`) would need
 * a different rule to be worth reading.
 */
const CHUNK_DIRS = ["assets", "js"];

/**
 * Paths escaping the chunk's own directory. Same-directory references are emitted by the
 * bundler and are hashed, so they cannot dangle; `../` is how an unresolved package-relative
 * path survives into the output.
 */
const ESCAPING_PATH = /["'`](\.\.\/[^"'`\n]*?\.(?:js|mjs|wasm))["'`]/g;

/**
 * References that survive into a browser chunk but are unreachable there, so they never
 * become a request. Each needs the runtime guard that makes it dead quoted alongside it —
 * an entry here is a claim about control flow, not a way to silence a finding.
 */
const KNOWN_INERT = new Set([
    // @spatialdata/core's parquet-wasm loader initialises from disk under Vitest/Node,
    // where undici cannot fetch a file:// URL. That branch is behind
    // `typeof process !== 'undefined' && process.versions?.node != null &&
    // typeof window === 'undefined'`, which no browser or worker satisfies — Vite defines
    // no global `process`, and the `node:fs` imports inside it resolved to
    // `__vite-browser-external` stubs. Only the string reaches the bundle.
    "../vendor/parquet-wasm/parquet_wasm_bg.wasm",
]);

/**
 * Only references landing inside the build tree are things the server would serve. Anything
 * resolving above it is a string that never becomes a request — Rolldown's CommonJS interop
 * registers modules under keys like `"../../node_modules/…/main.js"`, which are inert.
 */
function isServedFromBuild(target: string, buildDir: string): boolean {
    const relative = path.relative(buildDir, target);
    return relative !== "" && !relative.startsWith("..") && !path.isAbsolute(relative);
}

function jsFilesIn(dir: string): string[] {
    if (!fs.existsSync(dir)) return [];
    return fs
        .readdirSync(dir, { recursive: true, encoding: "utf8" })
        .filter((entry) => entry.endsWith(".js") || entry.endsWith(".mjs"))
        .map((entry) => path.join(dir, entry))
        .filter((file) => fs.statSync(file).isFile());
}

function danglingReferences(buildDir: string): string[] {
    const dangling: string[] = [];
    for (const chunkDir of CHUNK_DIRS) {
        for (const file of jsFilesIn(path.join(buildDir, chunkDir))) {
            const source = fs.readFileSync(file, "utf8");
            for (const [, reference] of source.matchAll(ESCAPING_PATH)) {
                if (KNOWN_INERT.has(reference)) continue;
                const target = path.resolve(path.dirname(file), reference);
                if (isServedFromBuild(target, buildDir) && !fs.existsSync(target)) {
                    dangling.push(`${path.relative(buildDir, file)} → ${reference}`);
                }
            }
        }
    }
    return dangling;
}

const builtDirs = BUILD_OUT_DIRS.map((dir) => path.join(repoRoot, dir)).filter((dir) =>
    fs.existsSync(path.join(dir, "assets")),
);

describe("build output asset references", () => {
    // CI builds immediately before running vitest, so finding nothing there means the
    // check silently stopped gating — worse than a failure. Locally a build is optional.
    test.runIf(process.env.CI)("a production build is present to check", () => {
        expect(builtDirs, `no build output found in ${BUILD_OUT_DIRS.join(", ")}`).not.toHaveLength(0);
    });

    test.skipIf(builtDirs.length === 0)("no chunk references a file the build did not emit", () => {
        const dangling = builtDirs.flatMap((dir) =>
            danglingReferences(dir).map((entry) => `${path.relative(repoRoot, dir)}/${entry}`),
        );
        const detail = `built chunks point at files missing from the output (rebuild if a listed\nbuild is just stale):\n${dangling.join("\n")}`;
        expect(dangling, detail).toEqual([]);
    });
});
