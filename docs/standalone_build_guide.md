
## Standalone Build Guide

The standalone build (`npm run standalone`) exposes the `mdv` object (including `ChartManager` and `DataLoaders`) on the `window`. It is up to you, the developer, to create dataloaders and listeners as needed.

### How to Use

1. **Build**: Run the standalone build script:
    ```
    npm run build-standalone
    ```
2. **Deploy**: Place the resulting `dist/mdv` folder into the `/static/` directory of your web application.
3. **Include in HTML**: Add the following to your HTML page:
    ```html
    <link rel="stylesheet" href="/static/mdv/assets/mdv.css">
    <script type="module" src="/static/mdv/mdv.js"></script>
    <script>
      // div: the container for MDV
      // dataSources: loaded from dataSources.json
      // dataLoader: custom dataloader
      // state: loaded from state.json
      // listener: reacts to MDV events
      const chartManager = new mdv.ChartManager(div, dataSources, dataLoader, state, listener);
    </script>
    ```

> **Note:** If you want to place the build in a different folder, set the `asset_base` environment variable to that folder (see below).

---

### Environment Variables

- **`worker_format`**: By default, this is `iife`. You can specify `es` for certain builds (e.g. JBrowse).
- **`exclude_dir`**: By default, the `examples` directory is included in the build. Set `exclude_dir=True` to exclude it.
- **`asset_base`**: By default, assets (such as webworkers) are loaded relative to the HTML page (`./`). If you want to serve MDV from different pages, specify an absolute asset base. This should match the location where you place the build folder.
- **`nofont`**: By default, FontAwesome is included in the build. Set `nofont=True` to exclude it (you will need to include FontAwesome manually in your HTML).

---




