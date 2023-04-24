import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react';
import { resolve } from 'path';


export default defineConfig(env => { return {
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
                obvios: "./index.html",
                index: resolve(__dirname, "src/static.html"),
            }
        }
    },
    plugins: [
        react({
            include: [/\.tsx?$/, /\.jsx?$/],
        })
    ],
}})