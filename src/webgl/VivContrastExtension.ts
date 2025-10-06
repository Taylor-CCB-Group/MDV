import { type Layer, LayerExtension, type UpdateParameters } from "@deck.gl/core";
import type { ShaderModule } from "@luma.gl/shadertools";
const N_CHANNELS = 6;
export type ContrastProps = {
    // order may matter here(?)
    redTest: number;
    greenTest: number;
    contrast: [number, number, number, number, number, number];
    brightness: [number, number, number, number, number, number];
};

const vivDeclCode = /*glsl*/ `
///----- VivContrastExtension decl
// in future we might have a #define for N_CHANNELS, at the viv level (xr-layer-fragment).
uniform VivContrastExtensionUniforms {
    float redTest;
    float greenTest;
    float contrast[${N_CHANNELS}];
    float brightness[${N_CHANNELS}];
} vivContrast;

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
    float c = vivContrast.contrast[channelIndex];
    float b = vivContrast.brightness[channelIndex];
    // return bias(b, gain(c, intensity));
    return intensity;
}
///---- end VivContrastExtension
`;

const vivProcessIntensityCode = /*glsl*/ `
///----- VivContrastExtension DECKGL_PROCESS_INTENSITY
intensity = apply_contrast_limits(intensity, contrastLimits);
intensity = clamp(intensity, 0., 1.);
intensity = applyBrightnessContrast(intensity, channelIndex);
///---- end VivContrastExtension
`;

const contrastModule = {
    name: "VivContrastExtension",
    uniformTypes: {
        // this is not how we type these arrays
        // we should probably start with a simpler extension that doesn't have arrays
        redTest: "f32",
        greenTest: "f32",
        contrast: "mat3x2<f32>", //will decode to `{ type: 'f32', components: 6 }`
        brightness: "mat3x2<f32>",
    },
    inject: {
        "fs:#decl": vivDeclCode,
        "fs:DECKGL_PROCESS_INTENSITY": vivProcessIntensityCode,
        'fs:DECKGL_MUTATE_COLOR': /*glsl*/`
        //---- VivContrastExtension DECKGL_MUTATE_COLOR
        rgba.r *= vivContrast.redTest;
        rgba.g *= vivContrast.greenTest;
        ///---- end VivContrastExtension
        `
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
        // should be arrays of N_CHANNELS... probably in imageSettingsStore.
        // --- would be nice if viv etc were better about managing variable number of channels ---
        if (contrast.length !== N_CHANNELS || brightness.length !== N_CHANNELS) {
            // throw new Error('contrast and brightness must be arrays of length 6');
            const contrastA = Array(N_CHANNELS);
            const brightnessA = Array(N_CHANNELS);
            for (let i = 0; i < N_CHANNELS; i++) {
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
                    "VivContrastExtension": {
                        // vivContrast: {
                        //     contrast: contrastA,
                        //     brightness: brightnessA,
                        // },
                        contrast: contrastA,
                        brightness: brightnessA,
                        redTest: 0.5,
                        greenTest: 0.5,
                    },
                });
            }
            return;
        }
        for (const model of this.getModels()) {
            // model.setUniforms({ contrast, brightness });
            model.shaderInputs.setProps({
                "VivContrastExtension": {
                    // vivContrast: {
                    //     contrast,
                    //     brightness,
                    // },
                    // re-using the arrays - is there a possibility that luma/deck.gl will not update the mutated contents?
                    contrast,
                    brightness,
                    redTest: 1.5,
                    greenTest: 0.9,
                },
            });
        }
    }
}
