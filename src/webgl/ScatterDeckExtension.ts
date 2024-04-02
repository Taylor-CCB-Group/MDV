import { LayerExtension } from "deck.gl/typed";



/**
 * As of this writing, used when we want to draw a scatterplot with square points.
 * WIP.
 * 
 * May encorporate other rendering features in future, and the current implementation
 * is likely to change.
 */
export class ScatterSquareExtension extends LayerExtension {
    static get componentName(): string {
        return 'ScatterSquareExtension';
    }
    getShaders() {
        return {
            inject: {
                'vs:#decl': `out float _lineWidthPixels;`,
                'vs:#main-end': `
                //---- ScatterSquareExtension
                // _outerRadiusPixels = outerRadiusPixels; // we might want this as an 'out', not used for now
                outerRadiusPixels = 0.; // forces all fragments to be 'inCircle'
                _lineWidthPixels = lineWidthPixels;
                ////
                
                `,
                'fs:#decl': `uniform float opacity;
                in float _lineWidthPixels;
                uniform float lineWidthScale;
                `,
                'fs:#main-end': `
                //---- ScatterSquareExtension
                // a bit wasteful to keep the original shader for circle rendering above
                // we should arrange layers/extensions differently, but YOLO
                if (stroked > 0.5) {
                    vec2 uv = abs(unitPosition);
                    float isLine = step(innerUnitRadius, max(uv.x, uv.y));
                    if (_lineWidthPixels <= 0.01) isLine = 0.;
                    fragmentColor = mix(vFillColor, vLineColor, isLine);
                }
                // make sure picking still works
                DECKGL_FILTER_COLOR(fragmentColor, geometry);
                //---- end ScatterSquareExtension
                `,
            }
        }
    }
}

export class ScatterDensityExension extends LayerExtension {
    static get componentName(): string {
        return 'ScatterDensityExension';
    }
    getShaders() {
        return {
            inject: {
                'fs:#decl': `uniform float radiusScale; uniform float opacity;`,
                'fs:#main-end': `
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
                fragmentColor.a = _a * opacity;
                //---- end ScatterDensityExension
                ////
                
                `,
            }
        }
    }
}