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
                // denom = 2*c^2 where c is the standard deviation - should be a uniform parameter
                float _a = pow(e, -(d*d)/(0.02));
                // fragmentColor.rgb = vec3(1.);
                fragmentColor.a = _a * opacity;
                //---- end ScatterDensityExension
                ////
                
                `,
            }
        }
    }
}