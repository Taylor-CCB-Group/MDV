import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react';
import type { RollupOptions } from 'rollup'; // Import RollupOptions from rollup
import * as path from 'path';

const flaskURL = "http://127.0.0.1:5051";

// setting output path: use --outDir
// todo review --assetsDir / nofont / cleanup & consolidate entrypoints

let hasWarned1 = false;
let hasWarned2 = false;

/**
 * shim for various different build configurations.
 * 
 * nb `dev_pt` which returns `{}` suffices for many things - vite devserver, netlify build.
 * `desktop_pt` just has a little more config setting a .ts input & making sure other things will work with flask template.
 * other methods are supposed to be for replacing other webpack configs.
 */
function getRollupOptions(): RollupOptions {
    const build = process.env.build as 'production' | 'dev_pt' | 'desktop' | 'desktop_pt';
    if (build === 'production') {
        return {
            input: process.env.nofont ? 'src/modules/basic_index_nf.js' : 'src/modules/basic_index.js',
            output: {
                entryFileNames: 'js/mdv.js',
                assetFileNames: (assetInfo) => {
                    //todo: match webpack behaviour with assetsDir / css-loader.
                    if (assetInfo.name.includes('index.css')) return 'assets/mdv.css';
                    return 'img/[name][extname]';
                },
            }
        }
    } else if (build === 'desktop_pt') {
        // currently there are different versions of entrypoint, this is the one with react etc.
        return {
            input: 'src/modules/static_index.ts',
            output: {
                entryFileNames: 'js/mdv.js',
                assetFileNames: (assetInfo) => {
                    if (assetInfo.name === 'static_index.css') return 'assets/mdv.css';
                    //not including hash, may impact caching, but more similar to previous webpack behavior
                    return 'img/[name][extname]';
                },
            }
        }
    } else if (build === 'dev_pt') {
        // version of vite build used for netlify deploy preview & devserver, using default 'index.html' entrypoint
        // (which as of writing refers to same static_index.ts as desktop_pt)
        return {}
    } else if (build === 'desktop') {
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
    } else {
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
                        if (assetInfo.name === name + '.css') return 'assets/mdv.css';
                        //not including hash, may impact caching, but more similar to previous webpack behavior
                        return 'img/[name][extname]';
                    },
                }
            }
        }
        throw new Error(`Unknown build type '${build}' and no VITE_ENTRYPOINT specified.`);
    }
}

export default defineConfig(env => { return {
    base: "./",
    server: {
        headers: {
            "Cross-Origin-Embedder-Policy": "require-corp",
            "Cross-Origin-Opener-Policy": "same-origin",
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "GET",
            "Access-Control-Allow-Headers": "X-Requested-With, content-type, Authorization",
        },
        port: 5170,
        strictPort: true,
        proxy: {
            // these routes are proxied to flask server in 'single project' mode
            '^/(get_|images|tracks|save).*': {
                target: flaskURL,
                changeOrigin: true,
            },
            '^/project/.*/\?': { //this works with ?dir=/project/project_name
                target: flaskURL,
                changeOrigin: true,
            },
            //will fail if url has search params <-- ?
            //will cause problems if we have json files that don't want to be proxied
            '^/.*\\.(json|b|gz)$': {
                target: flaskURL,
                changeOrigin: true,
            },
        }
    },
    publicDir: 'examples', //used for netlify.toml??... the rest is noise.
    build: {
        sourcemap: true,
        rollupOptions: { 
            ...getRollupOptions(),
            external: ['./python/**'],
            // this is annoying... lots of warnings in console output otherwise...
            // https://github.com/vitejs/vite/issues/15012#issuecomment-1815854072
            onLog(level, log, handler) {
                if (log.cause && (log.cause as any).message === `Can't resolve original location of error.`) {
                    if (hasWarned1) return;
                    hasWarned1 = true;
                    console.warn('Ignoring "Can\'t resolve original location of error." warnings... see comments in vite.config.mts');
                    return;
                }
                // there are still lots of other warnings, particularly from '@loaders.gl'
                // `"requireFromFile" is not exported by "__vite-browser-external"` etc...
                // I think to do with parts of that codebase that are have node code that won't be hit at runtime.
                // Tried to update deck.gl & luma.gl, which are responsible for @loaders.gl being included,
                // but that leads to other errors...
                // (I think because of older @deck.gl/core=8.8.27 being a peerDependency of @vivjs).
                // `RollupError: "_deepEqual" is not exported by "node_modules/@deck.gl/core/dist/esm/index.js", 
                //  imported by "node_modules/@deck.gl/extensions/dist/esm/collision-filter/collision-filter-effect.js"`
                // (in the case of that particular issue, the `_deepEqual` seems to have been added to core at the same 
                // time as collision-filter-effect.js which uses it, but we have an indirect dependency on older version...)
                if (log.message.includes('@loaders.gl')) {
                    if (hasWarned2) return;
                    hasWarned2 = true;
                    console.warn('Ignoring "@loaders.gl" warnings... see comments in vite.config.mts');
                    return;
                }
                handler(level, log)
            }
        },
    },
    plugins: [
        react({
            include: [/\.tsx?$/, /\.jsx?$/],
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
}})