import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react';
import type { RollupOptions } from 'rollup'; // Import RollupOptions from rollup

const flaskURL = "http://127.0.0.1:5051";

// setting output path: use --outDir
// todo review --assetsDir / nofont / cleanup & consolidate entrypoints


function getRollupOptions(): RollupOptions {
    const build = process.env.build as 'production' | 'netlify' | 'desktop' | 'desktop_pt';
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
    } else if (build === 'netlify') {
        // version of vite build used for netlify deploy preview, using default 'index.html' entrypoint
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
        throw new Error(`Unknown build type: ${build}`);
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
        rollupOptions: getRollupOptions(),
    },
    plugins: [
        react({
            include: [/\.tsx?$/, /\.jsx?$/],
        })
    ],
}})