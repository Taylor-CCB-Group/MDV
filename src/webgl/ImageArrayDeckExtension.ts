import {
    type Layer,
    type LayerContext,
    LayerExtension,
    type UpdateParameters,
} from "@deck.gl/core";
import type { ImageArray } from "./ImageArray";

export type ImageArrayExtensionProps = { imageArray: ImageArray; saturation: number };
export class ImageArrayDeckExtension<
    T extends ImageArrayExtensionProps = ImageArrayExtensionProps,
> extends LayerExtension<T> {
    static get componentName(): string {
        return "ImageArrayExtension";
    }
    getShaders() {
        // todo, make tests / test-project for this and review.
        // clearly will need to change for UBO but not going to do that without checking the results.
        return {
            inject: {
                "vs:#decl": `
        //---- ImageArrayExtension
        in float imageIndex;
        in float imageAspect;
        out float vImageIndex;
        out float vImageAspect;
        ////
        `,
                "vs:#main-start": `
        //---- ImageArrayExtension
        vImageIndex = imageIndex;
        vImageAspect = imageAspect;
        ////
        
        `,
                "vs:#main-end": `
        //---- ImageArrayExtension
        outerRadiusPixels = 1e-3;//HACKACK to make everything 'inCircle' for now
        // - but it breaks opacity, and generally we can do better...
        // it would be good to have border lines back, too
        ////
        
        `,
                "fs:#decl": `
        //---- ImageArrayExtension
        uniform mediump sampler2DArray imageArray;
        in float vImageIndex;
        in float vImageAspect;
        uniform float opacity;
        uniform float saturation;

        vec3 rgb2hsv(vec3 c){
            vec4 K = vec4(0., -1./3., 2./3., -1.);
            vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
            vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));
            float d = q.x - min(q.w, q.y);
            float e = 1.0e-10;
            return vec3(abs(q.z + (q.w - q.y) / (6. * d + e)), d / (q.x + e), q.x);
        }
        vec3 hsv2rgb(vec3 c){
            vec4 K = vec4(1., 2./3., 1./3., 3.);
            vec3 p = abs(fract(c.xxx + K.xyz) * 6. - K.www);
            return c.z * mix(K.xxx, clamp(p - K.xxx, 0., 1.), c.y);
        }
        
        `,
                "fs:DECKGL_FILTER_COLOR": `
        //---- ImageArrayExtension
        vec2 uv = 0.5 * (geometry.uv + 1.0);
        uv.y = 1. - uv.y; // flip y
        // todo fix aspect ratio in vertex shader
        uv.x /= min(1., vImageAspect);
        uv.y *= max(1., vImageAspect);
        if (uv.y > 1. || uv.x > 1.) discard;
        vec3 uvw = vec3(uv, vImageIndex);
        vec4 t = texture(imageArray, uvw);
        vec3 c = rgb2hsv(t.rgb);
        c.y *= saturation;
        t.rgb = hsv2rgb(c);
        color *= t;
        ///--- opacity may well not be correct gamma etc
        color.a = t.a * opacity; //HACK so broken inCircle doesn't break opacity
        // color.r = vImageAspect - 0.5;
        // vec3 s = vec3(textureSize(imageArray, 0));
        // uvw.z /= s.z;
        // color.rgb = uvw;
        // color.r = vImageIndex / s.z;
        ////

        `,
            },
        };
    }
    initializeState(
        this: Layer<ImageArrayExtensionProps>,
        context: LayerContext,
        extension: this,
    ): void {
        this.getAttributeManager()?.addInstanced({
            imageIndex: { size: 1, accessor: "getImageIndex", defaultValue: 0 },
            imageAspect: {
                size: 1,
                accessor: "getImageAspect",
                defaultValue: 1,
            },
        });
        this.props.saturation
    }
    updateState(
        this: Layer<ImageArrayExtensionProps>,
        params: UpdateParameters<Layer<ImageArrayExtensionProps>>,
    ) {
        // const { texture } = params.props.imageArray;
        const { saturation } = params.props;
        for (const model of this.getModels()) {
            const device = model.device;
            // todo clean this up
            if (!(device.type === "webgl")) throw new Error("non-webgl device not implemented");
            params.props.imageArray.wrapLumaTexture(device);
            const { lumaTexture } = params.props.imageArray;
            if (!lumaTexture) throw new Error("no luma texture");
            // const gl = (device as any).gl as WebGL2RenderingContext;
            // gl.activeTexture(gl.TEXTURE10);
            // gl.bindTexture(gl.TEXTURE_2D_ARRAY, texture);
            model.setUniforms({ saturation });
            model.setBindings({
                imageArray: lumaTexture,
            })
        }
    }
}
