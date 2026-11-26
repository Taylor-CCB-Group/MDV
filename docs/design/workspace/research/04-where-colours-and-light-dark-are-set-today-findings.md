Ticket: [04 Where colours and light-dark are set today](../tickets/04-where-colours-and-light-dark-are-set-today.md)
Map: [Workspace in MDV - Map](../README.md)
Source: MDV repo at `main` `b3c5a868` (2026-11-25). Paths are relative to the repo root. Line numbers are for that commit.

## Summary

- Colours come from four places that do not reference each other: the legacy underscore variables (two hand-synced CSS files), the shadcn variables in `src/tailwind.css`, two separate MUI themes, and colours written directly in components.
- Light or dark is a class on `<html>`. The project page sets it from `ChartManager.setTheme`. The catalog and login pages set it from `CustomThemeProvider`. Both start from the OS setting.
- The choice is never stored. A reload, a page change or an OS change resets it.
- Both static outputs load the same `all_css.js` bundle as the hosted project page. They have no catalog CSS and no catalog MUI theme.

## 1. Where colour variables are defined

### `src/charts/css/charts.css` (project page and static builds)

Line 1 says: sync any variable change with `catalog/catalog_index.css` "(temporary fix)".

| Block | Lines | Variables |
|---|---|---|
| `@media (prefers-color-scheme: dark) :root` | 2-29 | 25 |
| `@media (prefers-color-scheme: light) :root` | 31-59 | 25 |
| `.dark` | 61-83 | 21 |
| `.light` | 85-106 | 20 |
| `:root` (both themes) | 107-117 | 8 |

- 33 distinct names. 30 are colours, counting the `--menu-gradient` gradient. 3 are sizes: `--title_bar_height`, `--resizer-size`, `--resizer-offset`.
- The OS-dark media block (2-29) is not a dark theme. It sets `--background_color: #f5f5f5`, `--menu_bar_color: #d6d8db`, `--primary_background: #ecf0f1` (light values) next to `--text_color: white` and `--main_panel_color: black`. The `.dark` class block (61-83) holds the real dark values (`#252525`, `#404040`, `#353536`).
- The class blocks override the media blocks. They have the same specificity as `:root` and come later in the file, and the app always sets a class. Some names exist only in the media blocks: `--title_bar_height`, `--resizer-size`, `--resizer-offset` in both; `--icon_color_error` is missing from `.dark` and `.light`; `--border_menu_bar_color` is missing from `.light`. Those fall back to the media block for the OS setting. The fallback values happen to be the same in both media blocks, so nothing visible changes.
- `:root` 107-117 holds single values for both themes: `--button_color #1266f1`, the success/info/warning fills and texts, and `--scroll_bar_thumb`.
- `.dark a { color: #1266f1 }` at 118-120 is a hard-coded link colour, not the variable.
- The rest of the file uses hex literals for alerts (`.ciview-alert-info/danger/warning`, 539-555), the notice pip (`#f2c94c`, 594), the select focus glow (`#088aff`, 140) and text-box section rules (`#2d6f5e` etc., from 748). These are the same in light and dark.

### `src/catalog/catalog_index.css` (catalog page only)

Line 5 says: sync any variable change with `charts/css/charts.css` "(temporary fix)". It has 19 distinct names in the same five-block layout: OS-dark media 6-27, OS-light media 28-49, `.dark` 51-69, `.light` 71-88, `:root` 89-92. It also has `@tailwind base/components/utilities` at 1-3 and `a { color: var(--button_color) }` at 94-96.

Drift against `charts.css` (checked by script, every block):

- Every value that exists in both files is identical. The drift is only missing names.
- Missing from the catalog file: `--row_even_background`, `--row_odd_background`, `--header_menu_background`, `--scroll_bar_track` (all blocks); `--title_bar_height`, `--resizer-size`, `--resizer-offset` (media blocks); `--background_color_success/info/warning`, `--text_color_success/info/warning`, `--scroll_bar_thumb` (`:root`).
- The catalog has the same OS-dark media block quirk (light backgrounds with white text) at 6-27.
- Link colour differs: the catalog colours every `a` with `--button_color` in both themes. The project page only colours links in `.dark`, with a literal.
- The catalog page does not load `src/tailwind.css`, so the shadcn variables (`--background`, `--muted-foreground` etc.) are undefined there. Only `catalog_index.tsx:3` imports CSS for that entry.

### `src/tailwind.css` and `tailwind.config.js` (shadcn layer)

- `src/tailwind.css` `:root` 6-36 and `.dark` 38-66 define 19 shadcn "slate" colours as bare HSL triples, plus `--radius: 0.5rem`.
- It applies them globally twice: `* { border-border }` and `body { bg-background text-foreground }` at 69-76, and again with `outline-ring/50` at 114-120.
- There is no `.light` block. Light is the `:root` default.
- `tailwind.config.js:3` sets `darkMode: ["class"]`, so `dark:` variants follow the `.dark` class on `<html>`. Lines 20-54 map the Tailwind colour names to `hsl(var(--…))`. Lines 55-59 derive `lg/md/sm` radii from `--radius`.
- `components.json` declares shadcn with `baseColor: "slate"` and `cssVariables: true`. Only three shadcn components exist: `src/components/ui/accordion.tsx`, `label.tsx`, `radio-group.tsx`.
- `src/tailwind.css` is imported only by `src/modules/all_css.js:1`.

### Other CSS custom properties

- `src/utilities/css/JsonDialogStyles.css:1-50` defines `--json-*` syntax colours. Dark is the default. Light is switched by `@media (prefers-color-scheme: light)` at 37, so it follows the OS, not the app class.
- `src/react/components/chartArrayLayout.css:4-7` defines layout sizes (`--mdv-chart-array-*`), no colours.
- `src/react/components/pointerCrosshair.ts:41-42` sets `--mdv-crosshair-dx/dy` positions from JS.
- No other file in `src` defines colour custom properties.

## 2. The two MUI themes

### Catalog and login: `src/ThemeProvider.tsx`

- `getDesignTokens` (21-67) sets a full palette per mode. Light: primary `#2c3e50`, secondary `#34495e`, background `#ecf0f1` / paper `#ffffff`, text `#2c3e50` / `#7f8c8d`. Dark: primary `#ffffff`, secondary `#222`, background `#353536` / paper `#404040`, text `#ffffff` / `#ffffff`. Font Roboto, `body1` 0.9rem, flat Paper shadow, no button uppercase.
- Mode comes from `useMediaQuery("(prefers-color-scheme: dark)")` (70-71). An effect resets the mode whenever the OS value changes (73-76).
- An effect writes the mode to `<html>` with `classList.remove("light","dark")` then `add(mode)` (78-84).
- `toggleColorMode` flips the mode in React state only (86-94). The catalog toggle button is `src/catalog/Dashboard.tsx:260-270`.
- Used by `src/catalog/catalog_index.tsx:15` and `src/login/login_index.tsx:17`. Nothing else uses it.
- The light and dark backgrounds match the legacy `--primary_background` values (`#ecf0f1`, `#353536`).

### Project page: `MaterialWrapper` in `src/react/react_utils.tsx`

- `MaterialWrapper` (25-66) is a MobX observer. It reads `chartManager.theme === "dark"` (27) and passes only `palette.mode` (34-36). It sets no colours, so MUI's default palette applies. Its job is mostly `container` default props for Popper/Popover/Dialog, small TextFields, 5px radio padding, and Roboto at `fontSize: 12` (37-56).
- It renders `CssBaseline enableColorScheme` (63).
- Every React portal gets its own copy, because `createMdvPortal` (91-113) wraps each root in `MaterialWrapper`.
- So catalog dialogs and project dialogs use different MUI palettes and font sizes.

### Other MUI theme use

- `src/react/components/ChatButtons.tsx:85` re-provides the current MUI theme inside a popout window.
- Components branch on `theme.palette.mode === "dark"` for one-off rgba values, for example `src/charts/dialogs/AnnotationDialogComponent.tsx:163,395-422`, `src/react/components/SettingsDialogComponent.tsx:856`, `src/react/components/ColumnRemovalImpactDialog.tsx:123`.

## 3. `ChartManager` theme code (`src/charts/ChartManager.js`)

- `themes` (69-82) is a JS copy of four colours per theme. Dark: `title_bar_color #222`, `main_panel_color black`, `text_color white`, `background_color #333`. Light: `white`, `#f1f1f1`, `black`, `#bababa`. The two `background_color` values do not match the CSS (`#252525`, `#f5f5f5`).
- `themes` is read in two places. 179: the default `color` of each data source. No reader of `ds.color` was found. 1629-1659: the spinner colour (`t.text_color`) while a chart loads.
- `getPreferredColorScheme` (85-93) reads `matchMedia("(prefers-color-scheme: dark)")`.
- `listenPreferredColorScheme` (94-99) calls `setTheme` on every OS change. The constructor calls both (153, 156). The comment at 154-155 notes this overrides an explicit choice.
- `theme` is a MobX observable and `setTheme` an action (196-197).
- `setTheme` (523-541) sets `this.theme`, then `document.getElementsByTagName("html")[0].className = theme` (528). That replaces every class on `<html>`, unlike the catalog's `classList` call. It then calls `themeChanged()` on each chart that has it.
- `_setupThemeContextMenu` (505-513) builds `this.themeMenu` with one item per theme (`__getMenuItem`, 515-521). Nothing in `src` shows `themeMenu`. It is built and never opened.
- The toggle the user sees is `src/charts/dialogs/ToggleTheme.tsx` in the menu bar (`src/react/components/MenuBarComponent.tsx:153`). It calls `cm.setTheme` (10-13). Its comment at 9 says the theme "should be `"system" | "dark" | "light"`, and stored in local storage".

### Who reacts to a theme change

- `themeChanged()` exists on two charts. `src/charts/GenomeBrowser.js:223-225` repaints its canvas. `src/charts/TableChart.js:338-340` is an empty stub marked deprecated (called once at 137).
- React code reads `chartManager.theme` through `useTheme()` in `src/react/hooks.ts:194-207` (a MobX autorun). Users: `SelectionOverlay.tsx:79,141`, `HistogramWidget.tsx:77`, `ColorChannelComponents.tsx:240`, `histogram/BrushRangeChip.tsx:20,45`, `hooks/useSlickGridReact.ts:135,708` (passes `darkMode` to SlickGrid).
- Direct reads of `chartManager.theme`: `react_utils.tsx:27`, `CategoryHeatMap/CategoryHeatmapComponent.jsx:134` (hard-coded dark hexes at 34, 212-234), `SelectionDialogComponent.tsx:480-483` (rgba histogram colours, not reactive: plain read inside an observer), `modules/project_bootstrap.tsx:82-84` (loading overlay, falls back to `matchMedia`).
- `src/react/components/textbox/mermaidUtils.ts:6-19` reads the `<html>` class and the `--text_color`, `--main_panel_color`, `--primary_background`, `--border_menu_bar_color` variables, with hex fallbacks.

### Charts that read colours in JS for canvas, WebGL or image export

- `src/charts/RowChart.js:49` reads `--main_panel_color` once for the word cloud. Comment at 48: "won't redraw if we change the theme".
- `src/webgl/WGL2DI.js:999` reads the computed `color` for the polygon brush stroke.
- `src/charts/GenomeBrowser.js:403` and `src/browser/panel.js:658` read the computed `color` for highlight and track text.
- `src/react/components/DeckScatterReactWrapper.tsx:160` reads the panel background for image export.
- `src/react/components/dotplot/DotPlotReactWrapper.tsx:230-233` and `multiline/MultiLineReactWrapper.tsx:217-220` set `--text_color` to black during export, then restore it.
- 13 literal `fillStyle`/`strokeStyle` assignments across 7 files, mostly genome-browser tracks (`src/browser/tracks.js`, `src/browser/extra/*`) and `WGL2DI.js`.

## 4. How much each colour layer is used in `src`

Counted with a script over `.ts/.tsx/.js/.jsx/.css` files tracked in git under `src`, excluding the Font Awesome folder. Counts are occurrences, not files.

### Legacy underscore variables (`var(--text_color)` etc.)

- 167 colour reads in 27 files (211 reads minus 44 size reads, all of which are in `charts.css`).
- CSS: `charts.css` 62, `utilities/css/ContextMenu.css` 27, `browser/css/browser.css` 3, `catalog_index.css` 1.
- Code: 74 in 23 files. Top: `charts/dialogs/DebugErrorComponent.tsx` 12, `react/components/AxisComponent.tsx` 9, `ChatDialogComponent.tsx` 5, `H5MetadataPreview.tsx` 5, `ErrorComponentReactWrapper.tsx` 5, `dotplot/DotPlotPanel.tsx` 5, `SelectionOverlay.tsx` 4.
- Most read: `--text_color` 51, `--background_color` 25, `--border_menu_bar_color` 19, `--tooltip_background_color` 8.
- Defined but never read through `var()`: `--title_bar_color`, `--header_menu_background`, and all six success/info/warning variables.
- `--row_*`, `--header_menu_background` and `--scroll_bar_*` are read only inside `charts.css`.

### shadcn tokens

- Token classes (`bg-background`, `text-muted-foreground`, `border-input` and so on): 4 uses outside `tailwind.css`, all in `src/components/ui/radio-group.tsx`.
- `hsl(var(--…))` written directly, mostly as Tailwind arbitrary values like `text-[hsl(var(--muted-foreground))]`: 47 in 6 files. `ColorChannelComponents.tsx` 34, `HistogramWidget.tsx` 6, `charts/SVGChart.js` 2 (axis text fill, 61 and 74), `spatialLayers/PointsLayerPanel.tsx` 2, `PointsFeatureFilterPanel.tsx` 2, `ShapesLayerPanel.tsx` 1.
- Not tokens, but the bigger Tailwind colour layer: raw palette classes (`bg-slate-900`, `text-gray-500`…) 212 in 16 files, plus `bg-white`/`text-black` 32 in 15 files, plus `dark:` variants 103 in 14 files. Top: `file_upload/FileUploadUI.tsx` 35, `TiffMetadataTable.tsx` 32, `ChatDialogComponent.tsx` 28, `ColumnPreview.tsx` 28, `modules/project_bootstrap.tsx` 22, `SpatialLayerDialogComponent.tsx` 14.

### MUI palette

- `theme.palette.*` in code: 26 in 11 files. Top: `AnnotationDialogComponent.tsx` 6, `TabHeader.tsx` 5, `RecycleBinDialog.tsx` 3, `ColumnRemovalImpactDialog.tsx` 3.
- Palette path strings in `sx` and props (`"text.secondary"`, `"primary.main"`, `"grey.100"`…): 150 in 42 files. Top: `ColumnRemovalImpactDialog.tsx` 15, `AnnotationDialogComponent.tsx` 12, `spatialLayers/PointsFeatureFilterPanel.tsx` 10, `ChatDialogComponent.tsx` 9, `catalog/ProjectCard.tsx` 8.
- `"divider"` as a colour: 13 in 7 files.
- `color="primary|secondary|error|…"` props: 76 in 40 files.
- On the project page all of these resolve against MUI's default palette, because `MaterialWrapper` sets no colours. On the catalog they resolve against `ThemeProvider.tsx`.

### Hard-coded colours

- Hex literals: 672 in 67 files. CSS 200 (7 files): `charts.css` 98, `catalog_index.css` 52, `ContextMenu.css` 18, `JsonDialogStyles.css` 16, `table/css/slickgrid.css` 11.
- Hex in code: 472 in 60 files. `datastore/DataStore.js` 113 (the palettes), `CategoryHeatmapComponent.jsx` 36, `SingleHeatMap.js` 26, `SingleHeatMapGroup.js` 26, `HeatMap.js` 25, `tests/react/legend/ColorLegend.test.tsx` 17, `ThemeProvider.tsx` 12, `DeepToolsHeatMap.js` 12, `IGVBrowser/igvUtils.ts` 12, `dotplot/dotPlotPalettes.ts` 12, `login/AccessDenied.tsx` 10.
- `rgb()`/`rgba()` literals: 141 in 47 files. Top: `charts.css` 17, `histogram/BrushRangeChip.tsx` 12, `AnnotationDialogComponent.tsx` 10, `SplatterPlotReactWrapper.tsx` 10.
- Named colours in style keys (`color: "white"` etc.): 59 in 28 files. Top: `spatial_stats/CrossPcfMatrix.tsx` 8, `browser/tracks.js` 6, `file_upload/FileUploadViews.tsx` 6.

The regexes are approximate. A palette string match can catch a non-MUI string, and the hex match can catch a non-colour token. Spot checks found none.

## 5. How light or dark is remembered

- It is not remembered anywhere. There is no theme key in `localStorage`, `sessionStorage`, cookies, project state or the server.
- `src/lib/browserSettings.ts` is the store for per-browser preferences in `localStorage`. It has seven settings (render backend, WebGPU trace, legacy charts, feature previews, convert legacy charts, chart drag, data quality notes). None is a theme.
- `src/catalog/utils/catalogPreferences.ts` stores catalog view, sort and columns. No theme.
- Each page starts from the OS setting: `ChartManager.js:153` on the project page, `ThemeProvider.tsx:70-71` on catalog and login.
- A manual toggle lasts until the next reload, page change or OS change. Both pages reset on an OS change (`ChartManager.js:156`, `ThemeProvider.tsx:73-76`).
- The `<html>` class is set in two ways. Project page: `className = theme` (`ChartManager.js:528`). Catalog and login: `classList` remove and add (`ThemeProvider.tsx:78-84`). No server template sets a class (`python/mdvtools/templates/index.html`, `page.html`, `login.html`, `static_page.html`), so the class appears only after the JS runs.
- Popout windows copy the opener's `<html>` class and follow it with a `MutationObserver` (`src/utilities/Popout.ts:31-41`, `src/react/components/PopoutWindow.tsx:69-77`).
- Before the JS runs, only the `prefers-color-scheme` media blocks apply. In OS dark mode that is the light-background, white-text block described in section 1.

## 6. How the static builds get CSS and theme

There are two static outputs.

- `pnpm build-standalone` (`package.json:15`) runs Vite with `build=production` and writes to `dist/mdv`. The entry is `src/modules/basic_index.js`, or `basic_index_nf.js` when `nofont` is set (`vite.config.mts:234-244`). `basic_index.js:1` imports `./all_css`. The CSS is emitted as `assets/mdv{version}.css` (`vite.config.mts:199-203`). The page that loads it is not part of the build.
- `basic_index_nf.js:1-7` imports `charts.css` and the other legacy CSS directly. It does not import `src/tailwind.css`, Font Awesome, the SlickGrid material theme or Roboto. So in that variant the shadcn variables are undefined.
- `MDVProject.convert_to_static_page` (`python/mdvtools/mdvproject.py:2375-2416`) copies `python/mdvtools/static` (the Flask build, entry `project_bootstrap.tsx`, which imports `./all_css` at line 3). It fills `templates/static_page.html` with the CSS files listed by `entry_assets` (2392-2397).
- `src/modules/all_css.js` loads `tailwind.css`, `charts.css`, `ContextMenu.css`, SlickGrid CSS, Font Awesome, `wgl2di.css`, `browser.css` and Roboto 300-700.
- Neither static output loads `catalog_index.css` or `ThemeProvider.tsx`. Theme comes from `ChartManager` and the OS setting, the same as the hosted project page. `static_page.html` sets no class and no colours.

## 7. Data colours

- `defaultPalette` is in `src/datastore/DataStore.js:2092-2197`. It has 104 colours. The first eight are ColorBrewer Set1 (`#377eb8`, `#4daf4a`, `#984ea3`, `#ff7f00`, `#ffff33`, `#a65628`, `#f781bf`, `#999999`).
- `getColumnColors` (1685-1743) uses it for text columns with no `colors`, repeating it past 104 values (1722-1737).
- Numeric columns use `ov.colors || c.colors || defaultIPalette` (1691). `defaultIPalette` is an 8-stop Spectral ramp at 2199.
- Overrides, in order of where they live:
  - Column `colors` in the datasource config (`src/charts/schemas/DataSourceSchema.md:52`), changed at run time by `setColumnColors` (`DataStore.js:1665`).
  - The Color Palette dialog (`src/charts/dialogs/ColorPaletteWrapper.tsx:51-77`, `ColorPaletteComponent.tsx`, `CustomPaletteDialog.tsx`, which takes typed hex lists).
  - Synced colours across linked data sources: `ChartManager._sync_colors` (544) driven by `sync_column_colors` (`DataSourceSchema.md:179`).
  - Per-chart `overideValues` (min, max, colours, `fallbackOnZero`) built in `BaseChart.getColorFunction` (`src/charts/BaseChart.ts:757-772`) and quantile trimming at 639-658. Also `color_overlay` (1009-1016) and `color_legend`.
  - Chart-specific palettes: `src/react/components/dotplot/dotPlotPalettes.ts:6`, hard-coded scales in the heat-map charts, `d3-scale-chromatic` in `DensityPlot.js`, and image channel colours in `VivViewerMDV.js`.
- Data colours do not change with light or dark.

## 8. The colleague's extraction against the code

The extraction was taken at `ca48803` (2026-11-20). None of the theme files (`charts.css`, `catalog_index.css`, `tailwind.css`, `tailwind.config.js`, `ThemeProvider.tsx`, `react_utils.tsx`, `ChartManager.js`) changed between that commit and `b3c5a868`.

Agrees with the code:

- Every value in `tokens.json` for the legacy, shadcn and catalog MUI layers.
- The alert, notice pip, focus glow and text-box section literals.
- The OS-dark media block quirk, and that the `.dark` class block is what the app uses.
- Contrast: `#1266f1` on `#252525` is 3.06:1. `#7f8c8d` on white is 3.48:1. White on `#1266f1` is 5.0:1.
- Inter is named but not shipped. `package.json:86` has only `@fontsource/roboto`.

Disagrees with the code, or leaves something out:

- The default palette has 104 colours, not 30.
- "MUI layer: React dialogs (`src/ThemeProvider.tsx`)". That theme applies only to the catalog and login pages. Project-page dialogs use `MaterialWrapper`, which sets no colours, so they get MUI's defaults. The `mui-*` tokens describe the catalog only.
- White on light `destructive` is 3.76:1 with pure white. The token actually used, `--destructive-foreground` (`hsl(210 40% 98%)`), gives 3.6:1. The extraction says 3.8:1.
- It does not mention `catalog_index.css`, its drift, or that the catalog has no shadcn variables.
- It does not record the raw Tailwind palette layer (212 classes, 103 `dark:` variants), which is larger than shadcn token use.
- It does not record the `ChartManager.themes` JS copy, whose `background_color` values differ from the CSS.
- "Lets people switch" is true, but the switch is not stored and is undone by an OS change.
- `button_color` on `--primary_background` in dark (the canvas behind charts) is 2.45:1. The extraction lists only the `background_color` pair.

## Not verified

- Behaviour in a running browser. All of this is from reading code.
- Which rule wins for the project page `<body>` background: Tailwind's `bg-background` from `tailwind.css` or MUI `CssBaseline` from each portal. Both set it.
- The first paint on the catalog. MUI `useMediaQuery` returns `false` on the first render by default, which would paint light before switching to dark. Not checked.
- Whether anything outside `src` (for example the desktop app or embedding pages) sets the `<html>` class.
