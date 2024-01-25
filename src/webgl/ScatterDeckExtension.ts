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