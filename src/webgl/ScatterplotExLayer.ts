import { ScatterplotLayer } from "@deck.gl/layers";

/** don't think we can prepend '#version 300 es' in LayerExtension,
 * so we use this as a hack */

// export class ScatterplotExLayer extends ScatterplotLayer<any, any> {
//     // ts suddenly complains if we don't have this constructor
//     // todo review types after viv/deck.gl update
//     // biome-ignore lint/complexity/noUselessConstructor: we get a tsc error if we remove this - hope to fix after viv/deck.gl update
//     constructor(props: any) {
//         super(props);
//     }
//     getShaders() {
//         const shaders = super.getShaders();
//         shaders.vs = `#version 300 es\n${shaders.vs}`;
//         shaders.fs = `#version 300 es\n${shaders.fs}`;
//         return shaders;
//     }
// }
export { ScatterplotLayer as ScatterplotExLayer };