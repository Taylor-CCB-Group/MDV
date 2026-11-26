Type: grilling
Status: open
Blocked by: [20 Walking skeleton](20-walking-skeleton.md)
Map: [Workspace in MDV - Map](../README.md)
Merges: [13 One token set and one light-dark setting](merged/13-one-token-set-and-one-light-dark-setting.md), [15 Undo and redo](merged/15-undo-and-redo.md), [16 Rollout](merged/16-rollout.md)

## Question

**One token set and one light/dark setting**
- what replaces the colour layers (legacy CSS variables, shadcn/Tailwind, MUI), and where is the one token set defined?
- how do the legacy variables map onto it?
- do data colours stay as `DataStore`'s palette?
- is Inter shipped, and does Font Awesome stay?
- what density does the workspace use?
- where is the light/dark choice stored: the browser only, or the user's preferences on the server?
- does the design system artifact become the published source of the new tokens?

**Undo and redo**
- is it in the first release?
- what does it cover: moving and resizing charts, adding and removing charts, settings, filters?
- does it need a single model of a view's state first, and does the skeleton have one?

**Rollout**
- settled in [18 One view, saved layout and screen sizes](18-one-view-saved-layout-and-screen-sizes.md) (ADR 0026): it reaches users as the Workspace preview preference, and saves to `views_v2.json` beside an untouched `views.json`
- what has to be true before the preview becomes the default, and then before the pane workspace and `views.json` are retired?
- at retirement, how are the remaining `views.json` views converted (on open, or in one batch), and which readers move to v2 (Python writers, ChatMDV, catalog thumbnail, static build)?
- can the static build show v2 views before then?

## Context

- Findings: [04 Where colours and light-dark are set today - findings](../research/04-where-colours-and-light-dark-are-set-today-findings.md). There are four colour layers that don't reference each other, plus about 470 hard-coded hex values. Light/dark is an `<html>` class taken from the OS and never stored. The data palette has 104 colours.
- Tokens cover the catalog page and the project view, so both stay in step. The catalog redesign itself is out of scope.
- ADR 0019 rejected flags for single extractions; replacing the whole workspace may be different.
- Three topics in one session: if it runs long, finish the topic in hand, record it, and leave the rest of this ticket open.
