# About the Documentation

The main documentation is  written in reStructuredText and created with Sphinx, which requires the following python packages (pip install).

* sphinx
* sphinx-autoapi
* sphinx-rtd-theme
* spinx-js

As of this writing, `sphinx-js` has a dependency on a conflicting version of `MarkupSafe` so cannot be installed in the same environment as `mdvtools` itself.

On Unix-like systems, an npm script is provided to install the necessary packages in a virtual environment. To use it, run

```
npm run docs-env-setup
```

Also, to add the JavaScript api, jsdoc needs to be in your PATH, so it is best to install it globally:-
```
npm install -g  jsdoc
```
 (this can be more of a faff than it should be, we may review this at some point, espcially given we want to reconfigure to build typescript documentation anyway)

## Building the docs

To build all the docs, including the Python and JavaScript APIs, on Unix-like systems use
```
npm run build-docs-nix
``` 
This will activate the virtual environment and build the docs.


or 
```
sphinx-build -M html docs/maindocs docs/maindocs/_build

```
Output will be generated in `docs/maindocs/_build/html`. This output can be viewed locally with `npm run serve-docs`.

To build the JavaScript API docs alone use
```
npm run build-jsdocs
```
and they will be be placed in docs/jsdocs/build 


## Developing using VSCode
If you want a live(ish) preview of an .rst file in VSCode, then do the following:-
* install reStructuredText (lextudio.com) and reStructuredText Syntax highlighting (Trond Snekvik)
* install esbonio (pip install)

Then add the following to settings.json in the .vscode folder
```json
{
    "esbonio.sphinx.buildDir" : "${workspaceFolder}/docs/maindocs/_build/html",
    "esbonio.sphinx.confDir"  : "${workspaceFolder}/docs/maindocs"
}
```

Pressing ctrl shft K will open up a preview of the current .rst page and  every time it is saved, the preview will be updated.



