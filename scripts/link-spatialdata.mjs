#!/usr/bin/env node
/**
 * Point MDV at a local SpatialData.js checkout (dist builds) via pnpm `link:`.
 *
 * Usage:
 *   pnpm link:spatialdata
 *   SPATIALDATA_ROOT=~/code/www/SpatialData.ts pnpm link:spatialdata
 *   pnpm unlink:spatialdata
 *
 * After linking, rebuild upstream when you change it:
 *   (cd "$SPATIALDATA_ROOT" && pnpm --filter @spatialdata/vis build)
 * or keep a watch running in that repo.
 */

import { execSync } from "node:child_process";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";
import { fileURLToPath } from "node:url";

const mdvRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const packageJsonPath = path.join(mdvRoot, "package.json");

const PACKAGE_DIRS = {
    "@spatialdata/avivatorish": "packages/avivatorish",
    "@spatialdata/core": "packages/core",
    "@spatialdata/layers": "packages/layers",
    "@spatialdata/react": "packages/react",
    "@spatialdata/vis": "packages/vis",
    zarrextra: "packages/zarrextra",
};

const PUBLISHED_RANGES = {
    "@spatialdata/avivatorish": "^0.5.0",
    "@spatialdata/core": "^0.5.0",
    "@spatialdata/layers": "^0.5.0",
    "@spatialdata/react": "^0.5.0",
    "@spatialdata/vis": "^0.5.0",
    zarrextra: "^0.4.0",
};

function expandHome(p) {
    if (p.startsWith("~/")) return path.join(os.homedir(), p.slice(2));
    return p;
}

function candidateRoots() {
    const fromEnv = process.env.SPATIALDATA_ROOT?.trim();
    const candidates = [];
    if (fromEnv) candidates.push(expandHome(fromEnv));
    candidates.push(
        path.resolve(mdvRoot, "../SpatialData.ts"),
        path.resolve(mdvRoot, "../SpatialData.js"),
        path.resolve(os.homedir(), "code/www/SpatialData.ts"),
        path.resolve(os.homedir(), "code/www/SpatialData.js"),
    );
    return candidates;
}

function resolveSpatialRoot() {
    for (const root of candidateRoots()) {
        const visPkg = path.join(root, "packages/vis/package.json");
        if (fs.existsSync(visPkg)) return root;
    }
    throw new Error(
        `Could not find SpatialData.js checkout. Set SPATIALDATA_ROOT (tried:\n  ${candidateRoots().join("\n  ")})`,
    );
}

function readPackageJson() {
    return JSON.parse(fs.readFileSync(packageJsonPath, "utf8"));
}

function writePackageJson(pkg) {
    fs.writeFileSync(packageJsonPath, `${JSON.stringify(pkg, null, 2)}\n`);
}

function assertDistBuilt(root) {
    const missing = [];
    for (const [name, rel] of Object.entries(PACKAGE_DIRS)) {
        const distIndex = path.join(root, rel, "dist/index.js");
        if (!fs.existsSync(distIndex)) missing.push(`${name} (${rel}/dist)`);
    }
    if (missing.length) {
        throw new Error(
            `Local SpatialData packages need a build before linking. Missing:\n  ${missing.join("\n  ")}\n` +
                `From ${root}: pnpm build`,
        );
    }
}

function link() {
    const root = resolveSpatialRoot();
    assertDistBuilt(root);
    const pkg = readPackageJson();
    pkg.dependencies ??= {};
    pkg.pnpm ??= {};
    pkg.pnpm.overrides ??= {};

    for (const [name, rel] of Object.entries(PACKAGE_DIRS)) {
        const target = path.join(root, rel);
        const spec = `link:${target}`;
        if (name in pkg.dependencies) pkg.dependencies[name] = spec;
        pkg.pnpm.overrides[name] = spec;
    }

    writePackageJson(pkg);
    console.log(`Linked @spatialdata/* + zarrextra → ${root}`);
    console.log("Running pnpm install…");
    execSync("pnpm install", { cwd: mdvRoot, stdio: "inherit" });
    console.log("\nLinked. Restart Vite (clear .vite cache if needed: rm -rf node_modules/.vite).");
}

function unlinkPackages() {
    const pkg = readPackageJson();
    pkg.dependencies ??= {};
    pkg.pnpm ??= {};
    pkg.pnpm.overrides ??= {};

    for (const [name, range] of Object.entries(PUBLISHED_RANGES)) {
        if (name in pkg.dependencies) pkg.dependencies[name] = range;
        if (pkg.pnpm.overrides[name]?.startsWith("link:")) {
            delete pkg.pnpm.overrides[name];
        }
    }

    writePackageJson(pkg);
    console.log("Restored published @spatialdata/* / zarrextra ranges");
    console.log("Running pnpm install…");
    execSync("pnpm install", { cwd: mdvRoot, stdio: "inherit" });
}

const unlink = process.argv.includes("--unlink") || process.argv.includes("unlink");
if (unlink) unlinkPackages();
else link();
