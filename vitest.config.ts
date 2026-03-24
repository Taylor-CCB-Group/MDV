/**
 * Vitest is pinned to 4.0.18 in package.json: 4.1+ can load app `.ts` through paths that use
 * Node `stripTypeScriptTypes`, which breaks MobX decorators / `accessor` / parameter properties
 * (e.g. link_utils). See vitest#9873 / #9876. Vitest 4.0.18 lists Vite 6–7; root `package.json` overrides align `vite` with Vite 8.
 *
 * Merge `vite.config.mts` so test transforms match production (Babel decorators, glsl, aliases).
 */
import { defineConfig, mergeConfig } from 'vitest/config';

export default defineConfig(async () => {
    // Path aliases are not applied by Node's dynamic import(); file URL resolves at runtime.
    const { default: viteExport } = await import(new URL('./vite.config.mts', import.meta.url).href);
    const base =
        typeof viteExport === 'function'
            ? await viteExport({ command: 'serve', mode: 'test' })
            : await Promise.resolve(viteExport);
    return mergeConfig(base, {
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
}) as import('vitest/config').ViteUserConfigExport;
