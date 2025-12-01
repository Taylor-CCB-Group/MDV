import { defineConfig } from 'vitest/config';
import react from '@vitejs/plugin-react';
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
        web: {
          enabled: true,
        },
        ssr: {
          enabled: true,
        },
      },
    },
  },
}); 