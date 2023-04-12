import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react';

export default defineConfig({
    server: {
        headers: {
            "Cross-Origin-Embedder-Policy": "require-corp",
            "Cross-Origin-Opener-Policy": "same-origin",
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
    },
    plugins: [
        react({
            include: [/\.tsx?$/, /\.jsx?$/],
        })
    ],
})