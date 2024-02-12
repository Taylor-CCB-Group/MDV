# About the Documentation

The main documentation is  written in reStructuredText and created with Sphinx, which requires the following python packages (pip install).

* sphinx
* sphinx-autoapi
* sphinx-rtd-theme
* spinx-js


Also, to to to add the JavaScript api, jsdoc needs to be in your PATH, so it is best to install it globally:-
```
npm install -g  jsdoc
```

## Building the docs

To build all the docs, including the Python and JavaScript APIs, use
```
npm run build-docs
``` 
or 
```
sphinx-build -M html docs/maindocs docs/maindocs/_build

```
and they will be generated in docs/maindocs/_build/html

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



