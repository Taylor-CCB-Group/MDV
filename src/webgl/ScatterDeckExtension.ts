import { type Layer, LayerExtension, type UpdateParameters } from "@deck.gl/core";
import type { ShaderModule } from "@luma.gl/shadertools";
import { editWGSL, isWebGPULayer, type WGSLLayerExtension } from "./wgslLayerExtension";

/**
 * As of this writing, used when we want to draw a scatterplot with square points.
 * WIP.
 *
 * May encorporate other rendering features in future, and the current implementation
 * is likely to change.
 */
export class ScatterSquareExtension extends LayerExtension implements WGSLLayerExtension {
    static get componentName(): string {
        return "ScatterSquareExtension";
    }
    getShaders(this: Layer) {
        if (isWebGPULayer(this)) return {};
        return {
            inject: {
                "vs:#decl": "out float _lineWidthPixels;",
                "vs:#main-end": `
                //---- ScatterSquareExtension
                // _outerRadiusPixels = outerRadiusPixels; // we might want this as an 'out', not used for now
                outerRadiusPixels = 0.; // forces all fragments to be 'inCircle'
                _lineWidthPixels = lineWidthPixels;
                ////
                
                `,
                "fs:#decl": "in float _lineWidthPixels;",
                "fs:#main-end": `
                //---- ScatterSquareExtension
                // a bit wasteful to keep the original shader for circle rendering above
                // we should arrange layers/extensions differently, but YOLO
                // (actually, this is really dodgy & it's embarassing not rendering quads correctly)
                if (scatterplot.stroked > 0.5) {
                    vec2 uv = abs(unitPosition);
                    float isLine = step(innerUnitRadius, max(uv.x, uv.y));
                    if (_lineWidthPixels <= 0.01) isLine = 0.;
                    fragColor = mix(vFillColor, vLineColor, isLine);
                }
                // make sure picking still works
                DECKGL_FILTER_COLOR(fragColor, geometry);
                //---- end ScatterSquareExtension
                `,
            },
        };
    }
    /**
     * On WebGPU, measure distance from the centre with the max norm instead of the
     * Euclidean one, so the layer's own edge smoothing and stroke draw a square.
     */
    transformWGSL(source: string) {
        return editWGSL(source, {
            label: "ScatterSquareExtension",
            anchor: "let distToCenter = length(varyings.unitPosition) * varyings.outerRadiusPixels;",
            position: "replace",
            code: `\
// ScatterSquareExtension
  let distToCenter = max(abs(varyings.unitPosition.x), abs(varyings.unitPosition.y)) * varyings.outerRadiusPixels;`,
        });
    }
}

/** WGSL twin of the GLSL block below; the variable must be named after the module. */
const scatterDensityWGSL = /* wgsl */ `\
struct ScatterDensityExensionUniforms {
  opacity: f32,
};

@group(0) @binding(auto) var<uniform> ScatterDensityExension: ScatterDensityExensionUniforms;
`;

export class ScatterDensityExension extends LayerExtension implements WGSLLayerExtension {
    static get componentName(): string {
        return "ScatterDensityExension";
    }
    getShaders(this: Layer<ScatterDensityExension>, extension: this) {
        if (isWebGPULayer(this)) {
            const module: ShaderModule = {
                name: "ScatterDensityExension",
                uniformTypes: { opacity: "f32" },
                source: scatterDensityWGSL,
            };
            return { modules: [module] };
        }
        return { ...super.getShaders(extension), modules: [
            {
                name: "ScatterDensityExension",
                uniformTypes: {
                    opacity: "f32",
                },
                inject: {
                    // todo - add a uniform for ~kernelSigma, use UBO for uniforms (previous version is broken).
                    // todo - change picking behavior
                    "fs:#decl": "uniform ScatterDensityExensionUniforms { float opacity; } scatterDensity;",
                    "fs:#main-end": `
                    //---- ScatterDensityExension
                    const float e = 2.718281828459045;
                    float d = length(unitPosition);
                    // kernalSigma relates to dst in px in muspan, but denom uniform should be pre-computed
                    // kernel = np.exp(-( dst**2 / ( 2.0 * kernelSigma**2 ) ) )
                    // denom = 2*c^2 where c is the standard deviation / kernelSigma
                    // for muspan default kernelRadius=150, kernelSigma=50 => c = 1/3
                    // 2*(1/3)**2 => 0.222...
                    float _a = exp(-(d*d)/(0.222222222));
                    // should have an instance attribute for weight
                    fragColor.a = _a * scatterDensity.opacity;
                    //---- end ScatterDensityExension
                    ////                    
                    `,
                },
            },
        ]};
    }
    /**
     * On WebGPU, replace the colour's alpha with the gaussian kernel just before deck
     * premultiplies it. Unlike the GLSL version this runs after the picking branch has
     * returned, so it does not alter picking colours.
     */
    transformWGSL(source: string) {
        return editWGSL(source, {
            label: "ScatterDensityExension",
            anchor: "fragColor = deckgl_premultiplied_alpha(fragColor);",
            position: "before",
            code: `\
  // ScatterDensityExension: see the GLSL version for the kernel constant
  let densityDistance = length(varyings.unitPosition);
  fragColor.a = exp(-(densityDistance * densityDistance) / 0.222222222) * ScatterDensityExension.opacity;`,
        });
    }
    updateState(this: Layer<ScatterDensityExension>, params: UpdateParameters<Layer<ScatterDensityExension>>) {
        const { props } = params;
        const { opacity } = props;
        for (const model of this.getModels()) {
            model.shaderInputs.setProps({
                "ScatterDensityExension": {
                    opacity,
                },
            });
        }
    }
}
