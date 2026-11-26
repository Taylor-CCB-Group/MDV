Type: research
Status: resolved
Blocked by: none
Map: [Workspace in MDV - Map](../README.md)

## Question

Building on the MDV design system artifact, list where colours and the theme are defined and applied:

- the variable blocks in `src/charts/css/charts.css` and `src/catalog/catalog_index.css`, and `src/tailwind.css`
- `src/ThemeProvider.tsx` (catalog MUI theme) and `MaterialWrapper` in `src/react/react_utils.tsx` (project MUI theme)
- `ChartManager` `themes` / `setTheme` / `getPreferredColorScheme`, and canvas charts that read colours in JS (`themeChanged`)
- how many places use each colour layer (legacy variables, shadcn tokens, MUI palette) across `src`
- how the light/dark choice is (not) remembered, and how the standalone static build gets its CSS

## Context

- Design system artifact (a record of today's tokens, with known contrast failures): https://claude.ai/code/artifact/a7e9d6e3-f9a7-4fef-b42b-948024ea8649
- Feeds [13 One token set and one light-dark setting](merged/13-one-token-set-and-one-light-dark-setting.md).

## Answer

Colours live in four layers that do not reference each other: legacy underscore variables in `charts.css` and a hand-synced `catalog_index.css` (same values, but the catalog copy is missing 14 names and has no shadcn variables), shadcn variables in `tailwind.css`, two MUI themes (a full palette on catalog and login, MUI defaults on the project page), and about 470 hex literals in code (672 counting CSS).
Light or dark is a class on `<html>`. `ChartManager.setTheme` sets it on the project page and `CustomThemeProvider` sets it on catalog and login. Both start from the OS setting and reset when the OS setting changes. The choice is never stored.
Use counts in `src`: about 167 legacy variable reads, 51 shadcn token uses (47 of them as raw `hsl(var(--…))`), about 265 MUI palette references, and 212 raw Tailwind palette classes. Only two charts implement `themeChanged`, and one of those is empty.
Both static outputs load the same `all_css.js` bundle as the hosted project page, with theme from `ChartManager`. The default category palette in `DataStore.js` has 104 colours, not 30 as the design-system extraction says.

Findings: [04 Where colours and light-dark are set today - findings](../research/04-where-colours-and-light-dark-are-set-today-findings.md)
