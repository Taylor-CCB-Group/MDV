import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react';
import { resolve } from 'path';


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
        fs: {
            allow: ["examples", "."],
        },
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
                obvios: resolve(__dirname, "src/obvios.html"),
            }
        }
    },
    plugins: [
        react({
            include: [/\.tsx?$/, /\.jsx?$/],
        })
    ],
}})