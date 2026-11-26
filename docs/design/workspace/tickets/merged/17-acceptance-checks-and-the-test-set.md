Type: grilling
Status: merged
Merged into: [20 Walking skeleton](../20-walking-skeleton.md)
Blocked by: [02 Everything that reads or writes view layouts](../02-everything-that-reads-or-writes-view-layouts.md)
Map: [Workspace in MDV - Map](../../README.md)

## Question

What must the build pass before each slice ships? Which real views and projects form the test set, and which automatic checks run on them: no overlap, nothing past the edge at each width, layout saved and reloaded unchanged, limits on redraws and frame rate?

## Context

The prototype's own check (8 charts after a drag and a resize: none outside the grid, none overlapping) is the starting shape.
