import { type Layer, LayerExtension, type UpdateParameters } from "deck.gl/typed";

export default class VivContrastExtension extends LayerExtension {
    static get componentName(): string {
        return 'VivContrastExtension';
    }
    // biome-ignore lint/complexity/noBannedTypes: deck.gl idiosyncracy
    getShaders(this: Layer<{}>, extension: this) {
        return {
            inject: {
                'fs:#decl': /*glsl*/`///----- VivContrastExtension decl
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
                'fs:DECKGL_PROCESS_INTENSITY': /*glsl*/`///----- VivContrastExtension DECKGL_PROCESS_INTENSITY
                intensity = apply_contrast_limits(intensity, contrastLimits);
                intensity = applyBrightnessContrast(intensity, channelIndex);
                ///---- end VivContrastExtension
                `,
                // 'fs:DECKGL_MUTATE_COLOR': /*glsl*/`
                // //---- VivContrastExtension DECKGL_MUTATE_COLOR
                // // contrast adjustment
                // `
            }
        }
    }
    // biome-ignore lint/complexity/noBannedTypes: deck.gl idiosyncracy
    updateState(this: Layer<{}>, params: UpdateParameters<Layer<{}>>, extension: this): void {
        // should be arrays of 6... probably in imageSettingsStore.
        // const { contrast, brightness } = this.props;
        const contrast = [0.75, 0.75, 0.75, 0.75, 0.75, 0.75];
        const brightness = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
        for (const model of this.getModels()) model.setUniforms({ contrast, brightness });
    }
}