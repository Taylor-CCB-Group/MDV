import {
    type Layer,
    LayerExtension,
    type UpdateParameters,
} from "@deck.gl/core";

export type ContrastProps = {
    contrast: number[];
    brightness: number[];
};

export default class VivContrastExtension extends LayerExtension<ContrastProps> {
    static get componentName(): string {
        return "VivContrastExtension";
    }
    getShaders(this: Layer<ContrastProps>, extension: this) {
        return {
            inject: {
                "fs:#decl": /*glsl*/ `///----- VivContrastExtension decl
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
        };
    }
    updateState(
        this: Layer<ContrastProps>,
        params: UpdateParameters<Layer<ContrastProps>>,
        extension: this,
    ): void {
        const { props } = params;
        const { contrast, brightness } = props;
        // should be arrays of 6... probably in imageSettingsStore.
        // --- would be nice if viv etc were better about managing variable number of channels ---
        if (contrast.length !== 6 || brightness.length !== 6) {
            // throw new Error('contrast and brightness must be arrays of length 6');
            const contrastA = Array(6);
            const brightnessA = Array(6);
            for (let i = 0; i < 6; i++) {
                //don't really car about default values for missing entries...
                contrastA[i] = contrast[i] ?? 1;
                brightnessA[i] = brightness[i] ?? 0;
            }
            for (const model of this.getModels())
                model.setUniforms({
                    contrast: contrastA,
                    brightness: brightnessA,
                });
            return;
        }
        for (const model of this.getModels())
            model.setUniforms({ contrast, brightness });
    }
}
