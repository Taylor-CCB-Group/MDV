import { type Layer, LayerExtension, type UpdateParameters } from "@deck.gl/core";

/**
 * As of this writing, used when we want to draw a scatterplot with square points.
 * WIP.
 *
 * May encorporate other rendering features in future, and the current implementation
 * is likely to change.
 */
export class ScatterSquareExtension extends LayerExtension {
    static get componentName(): string {
        return "ScatterSquareExtension";
    }
    getShaders() {
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
}

export class ScatterDensityExension extends LayerExtension {
    static get componentName(): string {
        return "ScatterDensityExension";
    }
    getShaders(this: Layer<ScatterDensityExension>, extension: this) {
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
