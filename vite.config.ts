import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react';
import * as dotenv from 'dotenv';
import { resolve } from 'path';

dotenv.config();

const MDV_STATIC_DIR = process.env.MDV_STATIC_DIR || resolve("../mdv_static");
console.log("MDV_STATIC_DIR", MDV_STATIC_DIR);
const STATIC_PROXY = process.env.MDV_STATIC_PROXY || "http://localhost:9000";

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
            allow: ["examples", "."],
        },
        proxy: {
            "/static": {
                target: STATIC_PROXY,
                secure: false,
                changeOrigin: true,
                rewrite(path) {
                    return path.replace(/^\/static/, "");
                },
            },
        }
    },
    resolve: {
        alias: {
            "/data": "/examples/data",
        }
    },
    publicDir: "examples",
    build: {
        outDir: "vite-dist",
        rollupOptions: {
            input: {
                index: "./index.html",
                test: "src/static.html",
            }
        }
    },
    plugins: [
        react({
            include: [/\.tsx?$/, /\.jsx?$/],
        })
    ],
})