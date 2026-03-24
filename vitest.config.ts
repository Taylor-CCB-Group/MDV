/**
 * Vitest is pinned to 4.0.18 in package.json: 4.1+ can load app `.ts` through paths that use
 * Node `stripTypeScriptTypes`, which breaks MobX decorators / `accessor` / parameter properties
 * (e.g. link_utils). See vitest#9873 / #9876. Peer range lists Vite 6–7; we use Vite 8 with --legacy-peer-deps.
 */
import { defineConfig } from 'vitest/config';
import react from '@vitejs/plugin-react';
import glsl from 'vite-plugin-glsl';
import * as path from 'node:path';
// This is the script that runs in Node.js and starts the browser
// needs to be revisited - haven't had time to get it working yet.
// import { BrowserTestDriver } from '@probe.gl/test-utils';
// new BrowserTestDriver().run({
//   server: {
//     // Bundles and serves the browser script
//     command: 'webpack-dev-server', //<<< vite?
//     arguments: ['--env.render-test']
//   },
//   headless: true
// });

export default defineConfig({
  plugins: [
    glsl(), // Handle GLSL shader imports in tests
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
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
    }
  },
  test: {
    environment: 'jsdom',
    globals: true,
    include: ['src/**/*.{test,spec}.{js,ts,jsx,tsx}'],
    deps: {
      optimizer: {
        client: {
          enabled: true,
        },
        ssr: {
          enabled: true,
        },
      },
    },
  },
}); 