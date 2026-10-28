#!/usr/bin/env node
/**
 * Point MDV at a local IntraSpatial checkout (dist builds) via pnpm `link:`.
 *
 * Usage:
 *   INTRASPATIAL_ROOT=~/code/intraspatial pnpm link:intraspatial
 *   pnpm unlink:intraspatial
 *
 * Unlink restores the published `@intraspatial/deck` range. Rebuild upstream when you
 * change it:
 *   (cd "$INTRASPATIAL_ROOT" && pnpm build:lib && cd packages/deck && pnpm build)
 *
 * `typegpu` is a peer of both packages and stays a real MDV dependency either way.
 * See `scripts/link-spatialdata.mjs` for the sibling SpatialData.js link.
 */

import { execSync } from "node:child_process";
import * as fs from "node:fs";
import * as os from "node:os";
import * as path from "node:path";
import { fileURLToPath } from "node:url";

const mdvRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const packageJsonPath = path.join(mdvRoot, "package.json");

/** package → directory within the checkout */
const PACKAGE_DIRS = {
    intraspatial: ".",
    "@intraspatial/deck": "packages/deck",
};
/**
 * Direct dependencies, and the range `unlink` restores for each; the others only need an
 * override. Keep in step with package.json: this is a separate copy of the range.
 */
const PUBLISHED_RANGES = {
    "@intraspatial/deck": "^0.1.0",
};

function expandHome(p) {
    if (p.startsWith("~/")) return path.join(os.homedir(), p.slice(2));
    return p;
}

function resolveRoot() {
    const fromEnv = process.env.INTRASPATIAL_ROOT?.trim();
    const candidates = fromEnv ? [expandHome(fromEnv)] : [path.resolve(os.homedir(), "code/intraspatial")];
    for (const root of candidates) {
        if (fs.existsSync(path.join(root, "packages/deck/package.json"))) return path.resolve(root);
    }
    throw new Error(
        `Could not find an IntraSpatial checkout with packages/deck. Set INTRASPATIAL_ROOT (tried:\n  ${candidates.join("\n  ")})`,
    );
}

function readPackageJson() {
    return JSON.parse(fs.readFileSync(packageJsonPath, "utf8"));
}

function writePackageJson(pkg) {
    fs.writeFileSync(packageJsonPath, `${JSON.stringify(pkg, null, 2)}\n`);
}

function assertDistBuilt(root) {
    const missing = Object.entries(PACKAGE_DIRS)
        .filter(([, rel]) => !fs.existsSync(path.join(root, rel, "dist/index.js")))
        .map(([name, rel]) => `${name} (${path.join(rel, "dist")})`);
    if (missing.length) {
        throw new Error(
            `Local IntraSpatial packages need a build before linking. Missing:\n  ${missing.join("\n  ")}\n` +
                `From ${root}: pnpm build:lib, then pnpm build in packages/deck`,
        );
    }
}

function link() {
    const root = resolveRoot();
    assertDistBuilt(root);
    const pkg = readPackageJson();
    pkg.dependencies ??= {};
    pkg.pnpm ??= {};
    pkg.pnpm.overrides ??= {};

    for (const [name, rel] of Object.entries(PACKAGE_DIRS)) {
        const spec = `link:${path.join(root, rel)}`;
        pkg.pnpm.overrides[name] = spec;
        if (name in PUBLISHED_RANGES) pkg.dependencies[name] = spec;
    }

    writePackageJson(pkg);
    console.log(`Linked intraspatial + @intraspatial/deck → ${root}`);
    console.log("Running pnpm install…");
    execSync("pnpm install", { cwd: mdvRoot, stdio: "inherit" });
    console.log("\nLinked. Restart Vite (clear .vite cache if needed: rm -rf node_modules/.vite).");
}

function unlinkPackages() {
    const pkg = readPackageJson();
    pkg.dependencies ??= {};
    pkg.pnpm ??= {};
    pkg.pnpm.overrides ??= {};

    for (const name of Object.keys(PACKAGE_DIRS)) {
        if (pkg.pnpm.overrides[name]?.startsWith("link:")) delete pkg.pnpm.overrides[name];
    }
    for (const [name, range] of Object.entries(PUBLISHED_RANGES)) pkg.dependencies[name] = range;

    writePackageJson(pkg);
    console.log("Restored the published @intraspatial/deck range");
    console.log("Running pnpm install…");
    execSync("pnpm install", { cwd: mdvRoot, stdio: "inherit" });
}

const unlink = process.argv.includes("--unlink") || process.argv.includes("unlink");
if (unlink) unlinkPackages();
else link();
