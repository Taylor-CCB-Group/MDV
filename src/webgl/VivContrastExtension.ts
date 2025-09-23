import { type Layer, LayerExtension, type UpdateParameters } from "@deck.gl/core";
import type { ShaderModule } from "@luma.gl/shadertools";

export type ContrastProps = {
    contrast: number[];
    brightness: number[];
};

const contrastModule = {
    name: "VivContrastExtension",
    inject: {
        "fs:#decl": /*glsl*/ `///----- VivContrastExtension decl
                // todo - UBO version.
                // may need future revision if we change the number of channels.
                // layout(std140) uniform ChannelData {
                //     float contrast[6];
                //     float brightness[6];
                // };
                uniform float contrast[6];
                uniform float brightness[6];

                //https://cis700-procedural-graphics.github.io/files/toolbox_functions.pdf
                float bias(float b, float t) {
                    return pow(t, log(b) / log(0.5));
                }
                float gain(float g, float t) {
                    if (t < 0.5) {
                        return bias(1. - g, 2. * t) / 2.;
                    } else {
                        return 1. - bias(1. - g, 2. - 2. * t) / 2.;
                    }
                }
                float applyBrightnessContrast(float intensity, int channelIndex) {
                    float c = contrast[channelIndex];
                    float b = brightness[channelIndex];
                    return bias(b, gain(c, intensity));
                }
                ///---- end VivContrastExtension
                `,
        "fs:DECKGL_PROCESS_INTENSITY": /*glsl*/ `///----- VivContrastExtension DECKGL_PROCESS_INTENSITY
                intensity = apply_contrast_limits(intensity, contrastLimits);
                intensity = clamp(intensity, 0., 1.);
                intensity = applyBrightnessContrast(intensity, channelIndex);
                ///---- end VivContrastExtension
                `,
        // 'fs:DECKGL_MUTATE_COLOR': /*glsl*/`
        // //---- VivContrastExtension DECKGL_MUTATE_COLOR
        // // contrast adjustment
        // `
    },
} as const satisfies ShaderModule<ContrastProps>;

export default class VivContrastExtension extends LayerExtension<ContrastProps> {
    static get componentName(): string {
        return "VivContrastExtension";
    }
    //"only called if attached to a primitive layer"
    getShaders(this: Layer<ContrastProps>, extension: this) {
        // we may need something a bit like this - how does `ShaderModule` work in relation to extensions?
        // return super.getShaders({..., modules: [contrastModule]}) //?
        return { ...super.getShaders(extension), modules: [contrastModule] };
    }
    updateState(this: Layer<ContrastProps>, params: UpdateParameters<Layer<ContrastProps>>, extension: this): void {
        const { props } = params;
        const { contrast, brightness } = props;
        // should be arrays of 6... probably in imageSettingsStore.
        // --- would be nice if viv etc were better about managing variable number of channels ---
        if (contrast.length !== 6 || brightness.length !== 6) {
            // throw new Error('contrast and brightness must be arrays of length 6');
            const contrastA = Array(6);
            const brightnessA = Array(6);
            for (let i = 0; i < 6; i++) {
                //don't really care about default values for missing entries...
                contrastA[i] = contrast[i] ?? 1;
                brightnessA[i] = brightness[i] ?? 0;
            }
            //(method) Layer<ContrastProps>.getModels(): Model[]
            for (const model of this.getModels()) {
                // todo - update this to fit current deck.gl syntax for setting buffer...
                //Instead of calling `model.setUniforms` (or `model.setBindings`) use `model.shaderInputs.setProps` to update the UBO with props
                // model.setUniforms({
                //     contrast: contrastA,
                //     brightness: brightnessA,
                // });
                //We've tried to make sure that this extension has an appropriate type parameter, so ShaderPropsT should be ContrastProps
                //However, although `this` has a generic parameter, `getModels` is returning `Model[]`
                // so I don't understand how `model.shaderInputs` is supposed to have the correct type.
                //definition in deck.gl: setProps(props: Partial<{[P in keyof ShaderPropsT]?: Partial<ShaderPropsT[P]>}>): void
                //LSP says
                //(method) ShaderInputs<Partial<Record<string, Record<string, unknown>>>>.setProps(props: Partial<{
                //[x: string]: Partial<Record<string, unknown> | undefined>;
                //}>): void

                //... but then again, the LSP seems a bit confused about basic indentation, let alone type generics.

                //model.shaderInputs.setProps({custom: CustomProps});

                model.shaderInputs.setProps({
                    // contrast: contrastA,
                    // brightness: brightnessA
                });
            }
            return;
        }
        for (const model of this.getModels()) {
            model.setUniforms({ contrast, brightness });
        }
    }
}
