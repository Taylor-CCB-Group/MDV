import { ShaderAssembler, type WGSLShaderAssembler } from "@luma.gl/shadertools";
import { recordSpatialPerf } from "./perf";

/**
 * Memoises luma's WGSL shader assembly, which it otherwise redoes for every `Model`.
 *
 * On WebGPU every new deck layer builds a `Model`, and every `Model` runs the whole WGSL
 * assembly (module injection, two preprocessor passes, two interface scans) — about 6ms
 * for a scatterplot, and identical each time for the same layer class. The shader and
 * pipeline caches after it hit; only this step is repeated. Tiled layers create a layer
 * per tile, and the points layer evicts and recreates tiles as a zoom crosses levels, so a
 * zoom on a points layer spent a third of each tile's setup here. GLSL assembly is cheap
 * enough not to matter.
 *
 * Interim: luma should cache this itself (docs/upstream/luma-wgsl-assembly-cache.md).
 *
 * The key is every input assembly reads. Modules and the device (via its `limits` and
 * `features`, both stable per device) are keyed by identity, the rest by value; an input
 * that cannot be serialised (a function) makes the call uncacheable rather than risk a
 * false hit. A hit returns copies of the layout objects, since a `Model` keeps them.
 */

type AssembleProps = Parameters<WGSLShaderAssembler["assembleWGSLShader"]>[0];
type Assembled = ReturnType<WGSLShaderAssembler["assembleWGSLShader"]>;

const MAX_ENTRIES = 256;

const objectIds = new WeakMap<object, number>();
let nextObjectId = 0;

function idOf(value: unknown): string {
    if (typeof value !== "object" || value === null) return String(value);
    let id = objectIds.get(value);
    if (id === undefined) {
        id = nextObjectId++;
        objectIds.set(value, id);
    }
    return `#${id}`;
}

class Uncacheable extends Error {}

function valueKey(value: unknown): string {
    return (
        JSON.stringify(value, (_key, item: unknown) => {
            if (typeof item === "function") throw new Uncacheable();
            return item;
        }) ?? ""
    );
}

function cacheKey(assembler: WGSLShaderAssembler, props: AssembleProps): string | null {
    if (typeof props.source !== "string") return null;
    const defaultModules: unknown = Reflect.get(assembler, "_defaultModules");
    try {
        return [
            props.platformInfo.type,
            props.platformInfo.gpu,
            idOf(props.platformInfo.limits),
            idOf(props.platformInfo.features),
            (props.modules ?? []).map(idOf).join(","),
            Array.isArray(defaultModules) ? defaultModules.map(idOf).join(",") : "",
            valueKey(Reflect.get(assembler, "_hookFunctions")),
            valueKey(props.defines),
            valueKey(props.inject),
            valueKey(props.pluginInjections),
            valueKey(props.pluginVertexInputs),
            valueKey(props.pluginVaryings),
            valueKey([props.vertexEntryPoint, props.fragmentEntryPoint, props.scanVertexAttributes, props.prologue]),
            props.source,
        ].join("\u0000");
    } catch (error) {
        if (error instanceof Uncacheable) return null;
        throw error;
    }
}

function copyOf(assembled: Assembled): Assembled {
    return {
        ...assembled,
        modules: [...assembled.modules],
        bindingAssignments: assembled.bindingAssignments.map((assignment) => ({ ...assignment })),
        bindingTable: structuredClone(assembled.bindingTable),
        shaderLayout: assembled.shaderLayout && structuredClone(assembled.shaderLayout),
    };
}

const installed = new WeakSet<WGSLShaderAssembler>();

/** Memoise `assembler.assembleWGSLShader`. Exported for tests; the app installs it on luma's shared assembler below. */
export function installWGSLAssemblyCache(assembler: WGSLShaderAssembler) {
    if (installed.has(assembler)) return;
    installed.add(assembler);
    const assemble = assembler.assembleWGSLShader.bind(assembler);
    const cache = new Map<string, Assembled>();
    assembler.assembleWGSLShader = (props) => {
        const t0 = performance.now();
        const key = cacheKey(assembler, props);
        const cached = key === null ? undefined : cache.get(key);
        if (key !== null && cached) {
            // refresh its place in the eviction order
            cache.delete(key);
            cache.set(key, cached);
            const copy = copyOf(cached);
            recordSpatialPerf("webgpu.assembleWGSL.cached", performance.now() - t0);
            return copy;
        }
        const assembled = assemble(props);
        if (key !== null) {
            cache.set(key, copyOf(assembled));
            if (cache.size > MAX_ENTRIES) {
                const oldest = cache.keys().next().value;
                if (oldest !== undefined) cache.delete(oldest);
            }
        }
        recordSpatialPerf("webgpu.assembleWGSL", performance.now() - t0);
        return assembled;
    };
}

// deck assembles every WebGPU layer's shaders with luma's shared WGSL assembler.
installWGSLAssemblyCache(ShaderAssembler.getDefaultShaderAssembler("wgsl"));
