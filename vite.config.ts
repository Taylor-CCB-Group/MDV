import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react';
import { resolve } from 'path';

const flaskURL = "http://127.0.0.1:5050";

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
        port: 5170,
        strictPort: true,
        proxy: {
            '^/(get_|images|tracks|save).*': {
                target: flaskURL,
                changeOrigin: true,
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