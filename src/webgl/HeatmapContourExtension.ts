import { LayerExtension, Layer } from "@deck.gl/core";
import type { Buffer, Device, Texture } from '@luma.gl/core';
import { Model } from '@luma.gl/engine';
import { type LayerContext, project32 } from '@deck.gl/core';
import type { _ConstructorOf } from "@deck.gl/core/";

import heatmapPresentationFS from "./shaders/heatmap-presentation-layer-fragment.glsl?raw";

import { log } from '@luma.gl/core';
// log.enable();
log.level = 0;

export type ExtraContourProps = { contourOpacity: number };

// ---------- copied from deck.gl triangle-layer
const triangleVs = /*glsl*/ `#version 300 es
#define SHADER_NAME heatp-map-layer-vertex-shader
uniform sampler2D maxTexture;
uniform float intensity;
uniform vec2 colorDomain;
uniform float threshold;
uniform float aggregationMode;
in vec3 positions;
in vec2 texCoords;
out vec2 vTexCoords;
out float vIntensityMin;
out float vIntensityMax;
void main(void) {
    gl_Position = project_position_to_clipspace(positions, vec3(0.0), vec3(0.0));
    vTexCoords = texCoords;
    vec4 maxTexture = texture(maxTexture, vec2(0.5));
    float maxValue = aggregationMode < 0.5 ? maxTexture.r : maxTexture.g;
    float minValue = maxValue * threshold;
    if (colorDomain[1] > 0.) {
        maxValue = colorDomain[1];
        minValue = colorDomain[0];
    }
    vIntensityMax = intensity / maxValue;
    vIntensityMin = intensity / minValue;
}`;


export type _TriangleLayerProps = {
    data: { attributes: { positions: Buffer; texCoords: Buffer } };
    colorDomain: number[];
    aggregationMode: string;
    threshold: number;
    intensity: number;
    vertexCount: number;
    colorTexture: Texture;
    maxTexture: Texture;
    weightsTexture: Texture;
};

export class TriangleLayerContours extends Layer<_TriangleLayerProps & ExtraContourProps> {
    static layerName = 'TriangleLayerContours';

    declare state: {
        model: Model;
        positions: Buffer;
        texCoords: Buffer;
    };

    getShaders() {
        return { vs: triangleVs, fs: heatmapPresentationFS, modules: [project32] };
    }

    initializeState({ device }: LayerContext): void {
        this.setState({ model: this._getModel(device) });
    }

    _getModel(device: Device): Model {
        const { vertexCount, data, weightsTexture, maxTexture, colorTexture } = this.props;

        return new Model(device, {
            ...this.getShaders(),
            id: this.props.id,
            bindings: { weightsTexture, maxTexture, colorTexture },
            attributes: data.attributes,
            bufferLayout: [
                { name: 'positions', format: 'float32x3' },
                { name: 'texCoords', format: 'float32x2' }
            ],
            topology: 'triangle-fan-webgl',
            vertexCount
        });
    }

    //@ts-expect-error `draw({uniforms})`
    draw({ uniforms }): void {
        const { model } = this.state;
        const { intensity, threshold, aggregationMode, colorDomain, contourOpacity } = this.props;
        model.setUniforms({
            ...uniforms,
            intensity,
            threshold,
            aggregationMode,
            colorDomain,
            contourOpacity
        });
        model.draw(this.context.renderPass);
    }
}

// ---------- copied from deck.gl triangle-layer
