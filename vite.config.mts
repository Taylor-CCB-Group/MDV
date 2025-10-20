import { defineConfig, type ProxyOptions } from 'vite';
// import vitePluginSocketIO from 'vite-plugin-socket.io';
import react from '@vitejs/plugin-react';
import type { RollupOptions } from 'rollup'; // Import RollupOptions from rollup
import * as path from 'node:path';
import { execSync } from 'node:child_process';
import * as fs from 'node:fs';
// zarrita needed a polyfill for Buffer - seems like a bug
// seems ok without as long we don't use ZipFileStore (marked experimental anyway)
// having the polyfill means the build works, but devserver fails with 'cannot import outside a module'
// (not in the code using zarrita, but in unrelated worker modules)
// import { nodePolyfills } from 'vite-plugin-node-polyfills'


const flaskURL = "http://127.0.0.1:5055";
const port = 5170;
// setting output path: use --outDir
// todo review --assetsDir / nofont / cleanup & consolidate entrypoints
// maybe also the various build configurations at some point.

/**
 * shim for various different build configurations.
 * 
 * nb `dev_pt` which returns `{}` suffices for many things - vite devserver, netlify build.
 * `desktop_pt` just has a little more config setting a .ts input & making sure other things will work with flask template.
 * other methods are supposed to be for replacing other webpack configs.
 */
function getRollupOptions(): RollupOptions {
    const build = process.env.build || "desktop_pt" as 'production' | 'dev_pt' | 'desktop' | 'desktop_pt';
    if (build === 'production') {
        // somewhat equivalent to original webpack production build - not the current 'production' with new features.
        return {
            input: process.env.nofont ? 'src/modules/basic_index_nf.js' : 'src/modules/basic_index.js',
            output: {
                entryFileNames: 'js/mdv.js',
                assetFileNames: (assetInfo) => {
                    //todo: match webpack behaviour with assetsDir / css-loader.
                    if (assetInfo?.name?.includes('index.css')) return 'assets/mdv.css';
                    return 'img/[name][extname]';
                },
            }
        }
    } if (build === 'desktop_pt') {
        // currently there are different versions of entrypoint, this is the one with react etc.
        // used for Flask...
        return {
            input: {
                'mdv': 'src/modules/static_index.ts',
                'catalog': 'src/catalog/catalog_index.tsx',
                'login': 'src/login/login_index.tsx',
            },
            output: {
                entryFileNames: 'js/[name].js',
                assetFileNames: (assetInfo) => {
                    if (assetInfo.name === 'mdv.css') return 'assets/mdv.css';
                    if (assetInfo.name === 'catalog.css') return 'assets/catalog.css';
                    
                    //not including hash, may impact caching, but more similar to previous webpack behavior
                    return 'img/[name][extname]';
                },
            }
        }
    } if (build === 'dev_pt') {
        // version of vite build used for netlify deploy preview & devserver, using default 'index.html' entrypoint
        // (which as of writing refers to same static_index.ts as desktop_pt)
        // nb, localhost:5170/catalog_dev.html works on devserver without needing to change this.
        return {}
    } if (build === 'desktop') {
        return {
            input: 'src/modules/desktop_index.js',
            output: {
                entryFileNames: 'js/mdv.js',
                assetFileNames: (assetInfo) => {
                    if (assetInfo.name === 'desktop_index.css') return 'assets/mdv.css';
                    //not including hash, may impact caching, but more similar to previous webpack behavior
                    return 'img/[name][extname]';
                },
            }
        }
    }
    if (process.env.VITE_ENTRYPOINT) {
        // If you want a custom entrypoint - in particular, in order to have a custom DataLoader for interfacing
        // with another backend, you can specify it with VITE_ENTRYPOINT environment variable, e.g.
        // `VITE_ENTRYPOINT=path/to/my_index.js npx vite build --outDir path/to/output`
        // (nb, we may change the logic in this config...)
        const {name} = path.parse(process.env.VITE_ENTRYPOINT);
        return {
            input: process.env.VITE_ENTRYPOINT,
            output: {
                entryFileNames: 'js/mdv.js',
                assetFileNames: (assetInfo) => {
                    if (assetInfo.name === `${name}.css`) return 'assets/mdv.css';
                    //not including hash, may impact caching, but more similar to previous webpack behavior
                    return 'img/[name][extname]';
                },
            }
        }
    }
    throw new Error(`Unknown build type '${build}' and no VITE_ENTRYPOINT specified.`);
}

// avoiding some repition by defining a proxyOptions object used for all proxied routes.
const proxyOptions = { target: flaskURL, changeOrigin: true };
// ... and then this is a bit more concise than 
const proxy = [
    '^/(get_|images|tracks|save|chat).*', // these routes are proxied to flask server in 'single project' mode
    '^/project/.*/\?', //this works with ?dir=/project/id/foo...
    //will fail if url has search params <-- ? (what will fail?)
    //will cause problems if we have json files that don't want to be proxied
    '^/.*\\.(json|b|gz)$',
    '/projects',
    '/create_project',
    '/import_project',
    '/export_project',
    '/delete_project',
    '/extension_config',
    '/enable_auth',
    '/api_root',
    '/rescan_projects',
    '/login_dev',
// biome-ignore lint/performance/noAccumulatingSpread: don't care about performance in vite config
].reduce((acc, route) => ({...acc, [route]: proxyOptions}), {}) as Record<string, ProxyOptions>;
// (failed) attempt to let this proxy without cors_allowed_origins wildcard on server
// using more specific socketio in python for now
proxy['/socket.io'] = {
    target: flaskURL.replace("http:", "ws:"),
    changeOrigin: true,
    ws: true,
    // rewriteWsOrigin: true, // copilot says this is not needed in the current version of socket.io
    // ^^ not helping either way...
}

// not sure how we should make the root route work with vite devserver...
// proxy['/'] = {
//     target: "/catalog_dev.html"
// };

export default defineConfig(env => {
    // For local development, try to get Git info. This is guarded by a check for the .git directory
    // to prevent errors in environments where git is not available (like during Docker build, where
    // even in dev environment, the build happens before .git is copied in, for cache purposes).
    if (fs.existsSync('.git')) {
        try {
            const commitDate = execSync('git log -1 --format=%cI').toString().trimEnd();
            const branchName = execSync('git rev-parse --abbrev-ref HEAD').toString().trimEnd();
            const commitHash = execSync('git rev-parse HEAD').toString().trimEnd();
            const lastCommitMessage = execSync('git show -s --format=%s').toString().trimEnd();
        
            process.env.VITE_GIT_COMMIT_DATE = commitDate;
            process.env.VITE_GIT_BRANCH_NAME = branchName;
            process.env.VITE_GIT_COMMIT_HASH = commitHash;
            process.env.VITE_GIT_LAST_COMMIT_MESSAGE = lastCommitMessage;
            process.env.VITE_GIT_DIRTY = execSync('git diff --quiet || echo "dirty"').toString().trimEnd();
        } catch (e) {
            console.error('Failed to get git info:', e);
        }
    }
    process.env.VITE_BUILD_DATE = new Date().toISOString();

    return ({
    base: "./",
    server: {
        headers: {
            "Cross-Origin-Embedder-Policy": "require-corp",
            "Cross-Origin-Opener-Policy": "same-origin",
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "GET",
            "Access-Control-Allow-Headers": "X-Requested-With, content-type, Authorization",
        },
        port,
        strictPort: true,
        proxy,
    },
    publicDir: 'examples', //used for netlify.toml??... the rest is noise.
    build: {
        sourcemap: true,
        rollupOptions: { 
            ...getRollupOptions(),
            external: ['./python/**'],
        },
    },
    plugins: [
        react({
            include: [/\.tsx?$/, /\.jsx?$/],
            babel: {
                plugins: [
                    [
                        "@babel/plugin-proposal-decorators",
                        {
                            version: "2023-05"
                        }
                    ]
                ]
            }
        })
    ],
    worker: {
        format: 'iife',
    },
    resolve: {
        alias: {
            "@": path.resolve(__dirname, "./src"),
        }
    }
})})