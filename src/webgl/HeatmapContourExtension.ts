import type {
    LayersList,
    UpdateParameters,
    _ConstructorOf,
} from "@deck.gl/core/";
import { LayerExtension, type Layer } from "@deck.gl/core";
import { HeatmapLayer } from "@deck.gl/aggregation-layers";
import heatmapPresentationFS from "./shaders/heatmap-presentation-layer-fragment.glsl?raw";
import { log } from '@luma.gl/core';
// log.enable();
log.level = 0;

const contourDecl = /*glsl*/ `
//---- HeatmapContourExtension decl

//uniforms for tweakable parameters

//function for contouring

//--------------------
`;

const contourFilterColor = /*glsl*/ `
//---- HeatmapContourExtension
// gl_FragColor = DECKGL_FILTER_COLOR(gl_FragColor, geometry);
// instead of getLinearColor, we want to apply a contouring function to weight...
// but for now, let's get this code injected and verify what scope we're in

//--------------------
`;

export default class HeatmapContourExtension extends LayerExtension {
    static get componentName(): string {
        return "HeatmapContourExtension";
    }
    getShaders() {
        return {
            inject: {
                "fs:#decl": contourDecl,
                "fs:DECKGL_FILTER_COLOR": contourFilterColor,
            },
        };
    }
}
type ExtraContourProps = { contourOpacity: number };

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
export class ExtendableHeatmapLayer extends HeatmapLayer<
    Uint32Array,
    ExtraContourProps
> {
    static get componentName(): string {
        return "ExtendableHeatmapLayer";
    }
    
    protected getSubLayerProps(sublayerProps?: { [propName: string]: any; id?: string; updateTriggers?: Record<string, any>; }) {
        const mainProps = this.props;
        const props = super.getSubLayerProps(sublayerProps);
        if (sublayerProps.id === 'triangle-layer') {
            props.contourOpacity = mainProps.contourOpacity;
        }
        return props;
    }
    // biome-ignore lint/complexity/noBannedTypes: banned types are the least of our worries here
    protected getSubLayerClass<T extends Layer<{}>>(
        subLayerId: string,
        DefaultLayerClass: _ConstructorOf<T>,
    ): _ConstructorOf<T> {
        const theClass = super.getSubLayerClass(subLayerId, DefaultLayerClass);
        if (subLayerId === "triangle") { // maybe cleaner just to have our own class in this case rather than mangling the prototype
            // return TriangleLayerEx as any; // could have props._subLayerProps to override
            if (!theClass.prototype.__originalGetShaders__) {
                console.log(
                    ">>> saving original getShaders()... this should only happen once...",
                );
                theClass.prototype.__originalGetShaders__ =
                    theClass.prototype.getShaders;
                const originalDraw = theClass.prototype.draw;
                //changed name of texture to weightTexture, now the shader compiles with version 300 es
                //should we be applying binding this in updateState rather than draw?
                theClass.prototype.draw = function (opts) {
                    //opts.uniforms.weightTexture = (this.props as any).texture;
                    const { model } = this.state;
                    const { props } = this;
                    model.setUniforms({ contourOpacity: props.contourOpacity });
                    model.setBindings({ weightTexture: this.props.weightsTexture });
                    return originalDraw.call(this, opts);
                };
            }
            const originalGetShaders =
                theClass.prototype.__originalGetShaders__;
            const myTriangleFS = heatmapPresentationFS;
            //thought we could get away without doing this every time we update the shader code...
            //but moving the shader to a separate module isn't enough to get this to work
            if (theClass.prototype.__lastShader !== myTriangleFS) {
                console.log(
                    ">>> applying new getShaders() to TriangleLayer prototype...",
                );
                theClass.prototype.__lastShader = myTriangleFS;
                theClass.prototype.getShaders = () => {
                    console.log(">>> getShaders called...");
                    const shaders = originalGetShaders.call(this);
                    //will we get the new updated shader code without needing to mangle the prototype every time?
                    //no, probably closed on the the old module version when HMR happens
                    shaders.fs = myTriangleFS;
                    return shaders;
                };
            }
        }
        return theClass;
    }
}
