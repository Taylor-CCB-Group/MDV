import { Layer, LayerContext, LayerExtension, UpdateParameters } from 'deck.gl/typed';
import { ImageArray } from './ImageArray';
import { ScatterplotLayer } from 'deck.gl/typed';

/** don't think we can prepend '#version 300 es' in LayerExtension, 
 * so we use this as a hack */
export class ScatterplotExLayer extends ScatterplotLayer<any> {
    // ts suddenly complains if we don't have this constructor
    constructor(props: any) {
        super(props);
    }
    getShaders() {
        const shaders = super.getShaders();
        shaders.vs = '#version 300 es\n' + shaders.vs;
        shaders.fs = '#version 300 es\n' + shaders.fs;
        return shaders;
    }
}


type ImageArrayExtensionProps = { imageArray: ImageArray };
export class ImageArrayDeckExtension<T extends ImageArrayExtensionProps = ImageArrayExtensionProps> extends LayerExtension<T> {
    static get componentName(): string {
        return 'ImageArrayExtension';
    }
    getShaders() {
        return {
            inject: {
                'vs:#decl': `
        //---- ImageArrayExtension
        in float imageIndex;
        in float imageAspect;
        out float vImageIndex;
        out float vImageAspect;
        ////
        `,
        'vs:#main-start': `
        //---- ImageArrayExtension
        vImageIndex = imageIndex;
        vImageAspect = imageAspect;
        ////
        
        `,
        'vs:#main-end': `
        //---- ImageArrayExtension
        outerRadiusPixels = 1e-3;//HACKACK to make everything 'inCircle' for now
        ////
        
        `,
        'fs:#decl': `
        //---- ImageArrayExtension
        uniform mediump sampler2DArray imageArray;
        in float vImageIndex;
        in float vImageAspect;
        ////
        
        `,
        'fs:DECKGL_FILTER_COLOR': `
        //---- ImageArrayExtension
        vec2 uv = 0.5 * (geometry.uv + 1.0);
        // todo fix aspect ratio
        // uv.x *= max(1., vImageAspect);
        // uv.y /= vImageAspect;
        vec3 uvw = vec3(uv, vImageIndex);
        color *= texture(imageArray, uvw);
        // vec3 s = vec3(textureSize(imageArray, 0));
        // uvw.z /= s.z;
        // color.rgb = uvw;
        // color.r = vImageIndex / s.z;
        ////

        `
            }
        };
    }
    initializeState(this: Layer<{}>, context: LayerContext, extension: this): void {
        this.getAttributeManager()?.addInstanced({
            imageIndex: { size: 1, accessor: 'getImageIndex', defaultValue: 0 },
            imageAspect: { size: 1, accessor: 'getImageAspect', defaultValue: 1 }
        });
    }
    updateState(this: Layer<{}>, params: UpdateParameters<Layer<ImageArrayExtensionProps>>) {
        const { texture } = params.props.imageArray;
        for (const model of this.getModels()) {
            const gl = model.gl as WebGL2RenderingContext;
            gl.activeTexture(gl.TEXTURE10);
            gl.bindTexture(gl.TEXTURE_2D_ARRAY, texture);
            model.setUniforms({ imageArray: 10 });
        }
    }
}