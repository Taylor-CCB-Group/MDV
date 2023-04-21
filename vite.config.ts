import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react';
import * as dotenv from 'dotenv';

dotenv.config();

const MDV_STATIC_DIR = process.env.MDV_STATIC_DIR || "../mdv_static";

export default defineConfig({
    server: {
        headers: {
            "Cross-Origin-Embedder-Policy": "require-corp",
            "Cross-Origin-Opener-Policy": "same-origin",
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "GET",
            "Access-Control-Allow-Headers": "X-Requested-With, content-type, Authorization",
        },
        fs: {
            allow: ["examples", ".", MDV_STATIC_DIR],
        },
    },
    resolve: {
        alias: {
            "/data": "/examples/data",
            "/static": MDV_STATIC_DIR,
        }
    },
    publicDir: "examples",
    build: {
        outDir: "vite-dist",
        rollupOptions: {
            input: {
                index: "./index.html",
                static: "src/static.html",
            }
        }
    },
    plugins: [
        react({
            include: [/\.tsx?$/, /\.jsx?$/],
        })
    ],
})