import { LayerExtension } from "deck.gl/typed";



/**
 * As of this writing, used when we want to draw a scatterplot with square points.
 * WIP.
 * 
 * May encorporate other rendering features in future, and the current implementation
 * is likely to change.
 */
export class ScatterDeckExtension extends LayerExtension {
    static get componentName(): string {
        return 'ScatterDeckExtension';
    }
    getShaders() {
        return {
            inject: {
                'vs:#main-end': `
                //---- ScatterDeckExtension
                outerRadiusPixels = 0.;
                ////
                
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