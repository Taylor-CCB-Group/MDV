Type: grilling
Status: merged
Merged into: [18 One view, saved layout and screen sizes](../18-one-view-saved-layout-and-screen-sizes.md)
Blocked by: [02 Everything that reads or writes view layouts](../02-everything-that-reads-or-writes-view-layouts.md)
Map: [Workspace in MDV - Map](../../README.md)

## Question

What does a view store for its layout?

- grid units and row height
- tablet and phone overrides, stored only when someone rearranges them
- a format version
- how existing absolute and gridstack views load
- what the Python writers (default views, template views, chart prototypes, spatial template) must write
- whether the format stays independent of the engine chosen in [06 Choose the layout engine](06-choose-the-layout-engine.md)

## Context

Today: `gsposition` / `gssize` for gridstack, `position` / `size` in pixels for absolute, `layout` per data source.
