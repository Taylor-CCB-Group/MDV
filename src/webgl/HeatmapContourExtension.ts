import type { _ConstructorOf } from "@deck.gl/core/typed";
import { LayerExtension } from "deck.gl";
import { HeatmapLayer, type Layer } from "deck.gl/typed";

// original weights-fs glsl shader
/*glsl*/`
varying vec4 weightsTexture;
// Epanechnikov function, keeping for reference
// float epanechnikovKDE(float u) {
//   return 0.75 * (1.0 - u * u);
// }
float gaussianKDE(float u){
    return pow(2.71828, -u * u / 0.05555) / (1.77245385 * 0.166666);
}
void main()
{
  float dist = length(gl_PointCoord - vec2(0.5, 0.5));
    if (dist > 0.5) {
        discard;
    }
    gl_FragColor = weightsTexture * gaussianKDE(2. * dist);
    DECKGL_FILTER_COLOR(gl_FragColor, geometry);
}
`;

// shader code executed on the GPU for drawing the aggregated heatmap
/*glsl*/`
#version 100

#define SHADER_TYPE_FRAGMENT

#define DEFAULT_GPU
// Prevent driver from optimizing away the calculation necessary for emulated fp64
#define LUMA_FP64_CODE_ELIMINATION_WORKAROUND 1
// Intel's built-in 'tan' function doesn't have acceptable precision
#define LUMA_FP32_TAN_PRECISION_WORKAROUND 1
// Intel GPU doesn't have full 32 bits precision in same cases, causes overflow
#define LUMA_FP64_HIGH_BITS_OVERFLOW_WORKAROUND 1

#if (__VERSION__ > 120)

    # define FEATURE_GLSL_DERIVATIVES
    # define FEATURE_GLSL_DRAW_BUFFERS
    # define FEATURE_GLSL_FRAG_DEPTH
    # define FEATURE_GLSL_TEXTURE_LOD

    // DEPRECATED FLAGS, remove in v9
    # define FRAG_DEPTH
    # define DERIVATIVES
    # define DRAW_BUFFERS
    # define TEXTURE_LOD

#endif // __VERSION

// FRAG_DEPTH = > gl_FragDepth is available
#ifdef GL_EXT_frag_depth
    #extension GL_EXT_frag_depth : enable
    # define FEATURE_GLSL_FRAG_DEPTH
    # define FRAG_DEPTH
    # define gl_FragDepth gl_FragDepthEXT
#endif

// DERIVATIVES = > dxdF, dxdY and fwidth are available
#ifdef GL_OES_standard_derivatives
    #extension GL_OES_standard_derivatives : enable
    # define FEATURE_GLSL_DERIVATIVES
    # define DERIVATIVES
#endif

// DRAW_BUFFERS = > gl_FragData[] is available
#ifdef GL_EXT_draw_buffers
    #extension GL_EXT_draw_buffers : require
    #define FEATURE_GLSL_DRAW_BUFFERS
    #define DRAW_BUFFERS
#endif
// TEXTURE_LOD = > texture2DLod etc are available
#ifdef GL_EXT_shader_texture_lod
    #extension GL_EXT_shader_texture_lod : enable

    # define FEATURE_GLSL_TEXTURE_LOD
    # define TEXTURE_LOD

#endif



precision highp float;
#define MODULE_FP32
// END MODULE_fp32

#define MODULE_GEOMETRY

#define SMOOTH_EDGE_RADIUS 0.5

struct FragmentGeometry {
    vec2 uv;
}
geometry;
float smoothedge(float edge, float x) {
    return smoothstep(edge - SMOOTH_EDGE_RADIUS, edge + SMOOTH_EDGE_RADIUS, x);
}
// END MODULE_geometry

#define MODULE_PROJECT
// END MODULE_project

#define MODULE_PROJECT32
// END MODULE_project32





void DECKGL_FILTER_COLOR(inout vec4 color, FragmentGeometry geometry) {

}
void DECKGL_MUTATE_COLOR(inout vec4 rgba, float intensity0, float intensity1, float intensity2, float intensity3, float intensity4, float intensity5, vec2 vTexCoord) {

}
void DECKGL_PROCESS_INTENSITY(inout float intensity, vec2 contrastLimits, int channelIndex) {

}
#define SHADER_NAME triangle-layer-fragment-shader
precision highp float;
uniform float opacity;
uniform sampler2D texture;
uniform sampler2D colorTexture;
uniform float aggregationMode;
varying vec2 vTexCoords;
varying float vIntensityMin;
varying float vIntensityMax;
vec4 getLinearColor(float value) {
    float factor = clamp(value * vIntensityMax, 0., 1.);
    vec4 color = texture2D(colorTexture, vec2(factor, 0.5));
    color.a *= min(value * vIntensityMin, 1.0);
    return color;
}
void main(void) {
    vec4 weights = texture2D(texture, vTexCoords);
    float weight = weights.r;
    if (aggregationMode > 0.5) {
        weight /= max(1.0, weights.a);
    }
    // discard pixels with 0 weight.
    if (weight <= 0.) {
        discard;
    }
    vec4 linearColor = getLinearColor(weight);
    linearColor.a *= opacity;
    gl_FragColor = linearColor;
}
`


const contourDecl = /*glsl*/`
//---- HeatmapContourExtension decl

//uniforms for tweakable parameters

//function for contouring

//--------------------
`

const contourFilterColor = /*glsl*/`
//---- HeatmapContourExtension
// gl_FragColor = DECKGL_FILTER_COLOR(gl_FragColor, geometry);
// instead of getLinearColor, we want to apply a contouring function to weight...
// but for now, let's get this code injected and verify what scope we're in

//--------------------
`


export default class HeatmapContourExtension extends LayerExtension {
    static get componentName(): string {
        return 'HeatmapContourExtension';
    }
    getShaders() {
        return {
            inject: {
                'fs:#decl': contourDecl,
                'fs:DECKGL_FILTER_COLOR': contourFilterColor,
            }
        }
    }
}

//struggling to get #version 300 es to work...
//bigger picture, we want a very different version of HeatmapLayer anyway.
const myTriangleFS = /*glsl*/`//#version 300 es
#define SHADER_NAME triangle-layer-fragment-contour-shader

precision highp float;

uniform float opacity;
uniform sampler2D weightTexture; //using name 'texture' caused indirect glsl compiler errors with version 300 es(?)
uniform sampler2D colorTexture;
uniform float aggregationMode;

varying vec2 vTexCoords;
varying float vIntensityMin;
varying float vIntensityMax;

vec4 getLinearColor(float value) {
  float factor = clamp(value * vIntensityMax, 0., 1.);
  vec4 color = texture2D(colorTexture, vec2(factor, 0.5));
  color.a *= min(value * vIntensityMin, 1.0);
  return color;
}

float contour(float value) {
    //placeholder pending getting surrounding code in better shape
    return mod(value, 1.0) < 0.1 ? 1.0 : 0.0;
}
float smoothContour(float value) {
    float width = 0.1;
    float w = .1;//fwidth(value); //todo derivatives
    float f = .5;
    float wa = smoothstep(0., w * f, mod(value * f, 1.));
    wa = 1. - max(smoothstep(1.-w, 1., wa), smoothstep(w, 0., wa));
    return smoothstep(0., 1., wa*0.5);
    // return contour(value);
}

void main(void) {
  vec4 weights = texture2D(weightTexture, vTexCoords);
  float weight = weights.r;

  if (aggregationMode > 0.5) {
    weight /= max(1.0, weights.a);
  }
  // discard pixels with 0 weight.
  if (weight <= 0.) {
    discard;
  }

  vec4 linearColor = getLinearColor(weight);
  float c = smoothContour(weight);
  linearColor = mix(linearColor, vec4(1.0, 1.0, 1.0, 1.0), c);
  //   linearColor = c * vec4(1.);
  linearColor.a *= opacity;
  linearColor.a = max(linearColor.a, c*2.);
  gl_FragColor =linearColor;
}
`;

/** Original HeatmapLayer doesn't seem to apply extensions...
 * To be fair, there is some ambiguity
 * as there's more than one shader they could be applied to.
 * 
 * Anyway, this is an attempt to make HeatmapLayer work such that
 * when our HeatmapContourExtension is applied, it will be used.
 * It may well not be a complete or sustainable solution.
 * 
 * Also note we likely want to have a different version that changes
 * other aspects of how the heatmap is rendered 
 * - i.e. with a kernel in screen pixels vs coordinates.
 */
export class ExtendableHeatmapLayer extends HeatmapLayer {
    getShaders(type) {
        // type is just used as `type === 'max-weights-transform'`
        // not to distinguish e.g. between weights transform and color rendering.
        // this only lets us apply the extension to weights transform shader,
        // which is not what we want...
        // we want to replace getLinearColor() with a contouring function
        const shaders = super.getShaders(type);
        shaders.vs = `#version 300 es\n${shaders.vs}`;
        shaders._fs = `#version 300 es\n${shaders._fs}`;
        return shaders;
    }
    // biome-ignore lint/complexity/noBannedTypes: banned types are the least of our worries here
    protected getSubLayerClass<T extends Layer<{}>>(subLayerId: string, DefaultLayerClass: _ConstructorOf<T>): _ConstructorOf<T> {
        const theClass = super.getSubLayerClass(subLayerId, DefaultLayerClass);
        if (subLayerId === 'triangle') {
            if (!theClass.prototype.__originalGetShaders__) {
                console.log(">>> saving original getShaders()... this should only happen once...");
                theClass.prototype.__originalGetShaders__ = theClass.prototype.getShaders;
                const originalDraw = theClass.prototype.draw;
                //changed name of texture to weightTexture, now the shader compiles with version 300 es
                //(although it then doesn't render anything?)
                theClass.prototype.draw = function (opts) {
                    opts.uniforms.weightTexture = (this.props as any).texture;
                    return originalDraw.call(this, opts);
                }
            }
            console.log(">>> applying new getShaders() to TriangleLayer prototype...");
            const originalGetShaders = theClass.prototype.__originalGetShaders__;
            if (theClass.prototype.__lastShader !== myTriangleFS) {
                theClass.prototype.__lastShader = myTriangleFS;
                theClass.prototype.getShaders = () => {
                    console.log(">>> getShaders called...");
                    const shaders = originalGetShaders.call(this);
                    //will we get the new updated shader code without needing to mangle the prototype every time?
                    //no, probably closed on the the old module version when HMR happens
                    shaders.fs = myTriangleFS;
                    return shaders;
                }
            }
        }
        return theClass;
    }
}

