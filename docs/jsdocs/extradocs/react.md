# React-based charts in MDV

In order to make a new type of `Chart` that will integrate with MDV, there are a few subtle details that need to be attended to such that:

- The chart is registered with the system, so that the user can instantiate it from the ‘Add Chart’ dialog etc.
- Config settings can be used as expected:
  - Tweaking values in the standard settings dialog causes appropriate changes to be rendered.
  - Settings are able to be saved and loaded in the same way as with other charts.
- Other changes in application state - for instance data selection filtering - triggers rendering.
- When running in a dev-server, changes to React code can be (mostly) seen in real-time through HMR, without disrupting state.

For each new chart type, an entry should be added to `BaseChart.types['ChartName']`, similar to other existing charts. Often[^1] this will mean that there will be a new JS class extending `BaseReactChart`, associated with this entry. In turn, associated with this will generally be a functional React component responsible for the actual React rendering.

`BaseReactChart` has a constructor similar to `BaseChart`, but with one additional `ReactComponentFunction` parameter.


[^1]: It is possible to have multiple descriptions of chart-types that use the same `"class"`, but with different resulting configuration, which is interpreted by the same JS class to render a wide variation of actual charts. Where possible, it may be useful to try to adopt this approach more widely to avoid the need to create extra boilerplate for code that is mostly react-based, gradually re-factoring existing functionality in a way that suits the React paradigm.

## State management

`useConfig()` can be used to return `config` object that can be used to store mutable state. This is a `Proxy` that will trigger re-rendering of the component when it is mutated. This is the recommended way to store state that is specific to a particular chart instance and that will be persisted when 'save view' is called. For volatile sate, we are using `useState()` and Zustand (with the React context API).

### Implementation details

In order to make the mutable `config` object interoperable with React, we use [MobX](https://mobx.js.org/) to make it `observable`. `BaseReactChart` is responsible for this instrumentation.

This means that in the React context, we can reference members of `config`, and the relevant parts of the component will be re-rendered when they change. These changes are incorporated into state-saving by MDV in the same way as other charts.

To follow the prescribed best-practice for MobX, we make it so that the things that will cause mutations to the `config` object are wrapped in `action`s. For example, the `SettingsDialog` will modify each `func` in the entries returned from `chart.getSettings()` as such. The benefit of this would probably be more relevant in instances where several aspects of state are modified at once, but it's a good practice to follow (we get shouted at in the console otherwise, which could be avoided with `configure({enforceActions: false})`).

In some instances, we can start to use MobX reactions for effects which would otherwise arise from calls to chart methods. For example, the `useFilteredIndices()` hook will react to changes in the `DataStore`, and update automatically without the need to implement `onDataFiltered()` on the chart subclass.

### At what cost?

The amount of boilerplate with the current approach is quite small. We should take care to be conscious where possible to reduce the surface-area of our code that is specific to MobX, so that it would be easier to replace it with something else if we wanted to.

There will be non-zero overhead to the additional wrapping by MobX of various things (as well as making what were once plain-old-objects appear less plain when inspecting in debuggers, and adding some noise to stacktraces etc). The actual impact in terms of compute cost is believed to be negligible, although it is possible that there are instances in which blanket `makeAutoObservable` - while convenient in terms of reducing the amount of boilerplate we have to write and maintain - could be better replaced by a more targeted approach. It is likely that at certain points in development, subtle mistakes could lead to performance issues that may be hard to track down, but this is probably true of any approach.


## HMR

HMR - Hot Module Reloading - is a feature of Vite that allows changes to code to be reflected in the running application without the need to refresh the page. This is particularly useful when working on React components, as it allows us to see changes in the UI without losing state. However, it is not without its limitations - and it can be easy to break it in subtle ways.

In order for HMR to work without requiring a broad refactor of the code, there are some particular characteristics to be aware of.

Vite (or other bundlers) will maintain a graph of dependencies between modules, and when a module changes, it will re-evaluate it, and any other modules that depend on it. If in that process any encountered module cannot be hot-reloaded, then the page will be refreshed. That means that regressions can be easily introduced by changing the way that modules are imported and exported - for instance, `hooks.ts` had been working fine (in the sense that any changes were reflected at runtime), but at some point adding an import from it to `context.ts` meant that it could no longer be hot-reloaded, and code-paths that previously worked would now cause the page to refresh.

### Avoiding refreshes

This is generally fine, but there are some cases where it can be problematic. For example, if a module exports a class, and that class is used to create an instance, then changing the code for that class will not be reflected in the instance, as it will have already been created with the old code. It is still possible in some cases to make new instances of the class with the updated code (for example, by adding a new version of the same chart with the "Add Chart" dialog).


See `VivMDVReact` for a concrete example of how to do this.

Note that the classes referenced in `BaseChart.types` entries **must not themselves be exported**, only referenced there indirectly. Changing the code in these classes will not apply to existing instances in pages running in the dev-server, although new instances will be created with the updated code. Changing the code in the React component function will in many cases be reflected in existing instances instantly.

```jsx
// untested example illustrating the approach
import { observer } from 'mobx-react-lite';
import { BaseReactChart } from '/src/react/BaseReactChart';
import { useChart } from "./context";

// in future, we may be able to avoid applying `observer` at this point, but make it part of the `BaseReactChart` implementation, and potentially change that implementation to not use mobx
const ReactGui = observer(() => {
  const chart = useChart(); // this is a hook that returns the chart instance, see `hooks.ts` for others.
  return <div>{chart.title}</div>;
});

class MyChart extends BaseReactChart {
  // constructor signature should match that of BaseChart, so that ChartManager can instantiate it.
  constructor(dataStore, div, config) {
    super(dataStore, div, config, ReactGui);
    // console.log('uncomment this line and add a new instance of this chart to see this message in the console')
    // you could also make all kind of other changes to the class, and/or the `BaseChart.types["MyChart"]` entry, and see them reflected in new instances of the chart, but not existing ones
  }
}

BaseChart.types["MyChart"] = {
    class: MyChart,
    description: "My Chart",
    // ...
}

export default "some token that is used to drive side-effects through the system..."
```

### Why does this work?

When the code for the above module changes, HMR will re-evaluate it, mutating the `BaseChart.types` so that new instances of `MyChart` will be created with the updated code. When new module code is integrated, it cannot reasonably be expected to change the implementation of an existing instance of a class while maintaining its state - but if it was `export`ed explicitly, the system would consider any instances of that class to be invalid, requiring the page to refresh in order for their state to be re-created consistant with the code that is now in place.

The different nature of the way state works in the React context allows those changes to be reflected in existing instances, without the need to refresh the page. Components and hooks should be able to be fairly freely imported and exported in whichever way is most convenient, without needing to worry about this.

It would be possible - but probably unwise - to implement a mechanism that would allow the code for a class to be updated in-place, but this could be a lot of work, and would probably be quite fragile. It would also be somewhat explicitly tied to Vite in a way that should be avoided where possible.

## Non-Chart components

