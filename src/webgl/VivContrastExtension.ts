import type { Layer, UpdateParameters } from "@deck.gl/core";
import { VivLayerExtension } from "@hms-dbmi/viv";
import { VIV_CHANNEL_INDEX_PLACEHOLDER as I } from "@vivjs/constants";

export type ContrastProps = {
    contrast: number[];
    brightness: number[];
};

const M = "vivContrast";
const DEFAULT_BRIGHTNESS = 0.5;
const DEFAULT_CONTRAST = 0.5;

export default class VivContrastExtension extends VivLayerExtension<ContrastProps> {
    static get componentName(): string {
        return "VivContrastExtension";
    }

    getVivShaderTemplates() {
        return {
            modules: [
                {
                    name: M,
                    uniformTypes: {
                        [`contrast${I}`]: "f32",
                        [`brightness${I}`]: "f32",
                    },
                    fs: /*glsl*/ `///----- VivContrastExtension decl
                uniform ${M}Uniforms {
                    float contrast${I};
                    float brightness${I};
                } ${M};

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
                    // could we change the signature here rather than expanding to all channels internally?
                    // hopefully the glsl stack will make this not as bad as it looks
                    // maybe one day luma.gl will understand arrays in UBOs so we can do
                    // "${M}.contrast[channelIndex]"
                    float contrast[NUM_CHANNELS] = float[NUM_CHANNELS](
                        ${M}.contrast${I},
                    );
                    float brightness[NUM_CHANNELS] = float[NUM_CHANNELS](
                        ${M}.brightness${I},
                    );
                    float c = contrast[channelIndex];
                    float b = brightness[channelIndex];
                    return bias(b, gain(c, intensity));
                }
                ///---- end VivContrastExtension
                `,
                    inject: {
                        "fs:DECKGL_PROCESS_INTENSITY": /*glsl*/ `///----- VivContrastExtension DECKGL_PROCESS_INTENSITY
                intensity = apply_contrast_limits(intensity, contrastLimits);
                intensity = clamp(intensity, 0., 1.);
                intensity = applyBrightnessContrast(intensity, channelIndex);
                ///---- end VivContrastExtension
                `,
                    },
                },
            ],
        };
    }

    updateState(
        this: Layer<ContrastProps>,
        params: UpdateParameters<Layer<ContrastProps>>,
        extension: this,
    ): void {
        super.updateState(params, extension);
        const { props } = params;
        const { contrast = [], brightness = [] } = props;
        const getNumChannels = "getNumChannels" in this && typeof this.getNumChannels === "function"
            ? this.getNumChannels.bind(this)
            : undefined;
        const numChannels = getNumChannels?.()
            ?? ("selections" in props && Array.isArray(props.selections) ? props.selections.length : 0);
        const vivContrastUniforms: Record<string, number> = {};

        for (let i = 0; i < numChannels; i++) {
            vivContrastUniforms[`contrast${i}`] = contrast[i] ?? DEFAULT_CONTRAST;
            vivContrastUniforms[`brightness${i}`] = brightness[i] ?? DEFAULT_BRIGHTNESS;
        }

        for (const model of this.getModels()) {
            model.shaderInputs.setProps({ [M]: vivContrastUniforms });
        }
    }
}
