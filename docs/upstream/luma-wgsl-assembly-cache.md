# luma.gl: every WebGPU Model re-assembles its WGSL, even for identical inputs (draft)

Draft issue for `@luma.gl/engine` / `@luma.gl/shadertools`, found in MDV's SpatialData chart
on WebGPU. Environment: Chrome (Dawn, `metal-3`, macOS), `@luma.gl/*` 9.4.2,
`@deck.gl/core` 9.4.0. Not posted.

## Summary

`new Model()` on a WebGPU device always calls `shaderAssembler.assembleWGSLShader(...)`. For
deck's scatterplot shader with its default modules that is about 6ms of main thread
(measured hot, replaying the same props: 5.7ms; ~8ms in a busy frame), and the result is
identical every time for the same layer class on the same device. The shader and pipeline
factories after it already cache by source, so assembly is the only part repeated.

Tiled deck layers create one sublayer, and so one `Model`, per tile. Zooming a tiled
points layer created 39–51 sublayers per zoom sweep, so assembly alone was 220–290ms of
a 60-frame zoom. The GLSL path is cheap enough (~0.1ms) not to matter.

The work inside one call is also partly doubled: `assembleWGSLShader` in
`assemble-shaders.js` already computes `bindingTable` (`getShaderBindingDebugRowsFromWGSL`)
and `shaderLayout` (`scanWGSLInterface`), and `WGSLShaderAssembler.assembleWGSLShader`
discards them and computes both again on the preprocessed source, after a second
`preprocess` pass over the whole assembled source.

## Suggestion

Memoise `WGSLShaderAssembler.assembleWGSLShader` on its inputs: source, module identities
(and default modules), defines, injections, plugin inputs/varyings, entry points, hook
functions, and the platform (limits and features are per device). Return copies of
`shaderLayout`, `bindingTable` and `bindingAssignments`, since a `Model` keeps them. The
binding registry is already stable per run, so a cached result matches a fresh one.
Separately, skip the debug binding table unless something asks for it, and avoid the
second scan.

MDV does this as an interim, wrapping luma's shared WGSL assembler
(`src/react/spatialdata/wgsl_assembly_cache.ts`): after it, a zoom over a tiled points
layer made one real assembly and 46–49 cache hits averaging 0.05ms.
