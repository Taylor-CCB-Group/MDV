import { Layer } from "@deck.gl/core";
import type { Buffer, Device, Texture } from '@luma.gl/core';
import { Model } from '@luma.gl/engine';
import { type LayerContext, project32 } from '@deck.gl/core';

import heatmapPresentationFS from "./shaders/heatmap-presentation-layer-fragment.glsl";

import { log } from '@luma.gl/core';
// log.enable();
log.level = 0;

export type ExtraContourProps = { 
    contourOpacity: number;
    contourFill: number;
    fillOpacity: number;
};

// ---------- copied from deck.gl triangle-layer
// (updated for UBO, not compared with new deck.gl implementation)
const triangleVs = /*glsl*/ `#version 300 es
#define SHADER_NAME heatp-map-layer-vertex-shader
uniform sampler2D maxTexture;
in vec3 positions;
in vec2 texCoords;
out vec2 vTexCoords;
out float vIntensityMin;
out float vIntensityMax;
void main(void) {
    gl_Position = project_position_to_clipspace(positions, vec3(0.0), vec3(0.0));
    vTexCoords = texCoords;
    vec4 maxTexture = texture(maxTexture, vec2(0.5));
    float maxValue = triangle.aggregationMode < 0.5 ? maxTexture.r : maxTexture.g;
    float minValue = maxValue * triangle.threshold;
    if (triangle.colorDomain[1] > 0.) {
        maxValue = triangle.colorDomain[1];
        minValue = triangle.colorDomain[0];
    }
    vIntensityMax = triangle.intensity / maxValue;
    vIntensityMin = triangle.intensity / minValue;
}`;

const contourUniformCode = /*glsl*/ `uniform triangleUniforms {
    float aggregationMode;
    vec2 colorDomain;
    float intensity;
    float threshold;
    float contourOpacity;
    float contourFill;
    float fillOpacity;
} triangle;
`

const triangleContourUniforms = {
    name: "triangle",
    vs: contourUniformCode,
    fs: contourUniformCode,
    uniformTypes: {
        aggregationMode: "f32",
        colorDomain: "vec2<f32>",
        intensity: "f32",
        threshold: "f32",
        contourOpacity: "f32",
        contourFill: "f32",
        fillOpacity: "f32",
    },
};


export type _TriangleLayerProps = {
    data: { attributes: { positions: Buffer; texCoords: Buffer } };
    colorDomain: number[];
    aggregationMode: number;
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
        // I see warnings about various moules not being found 
        // (picking - also layer, shadow, lighting, phongMaterial, gouraudMaterial), as well as
        // > Ignoring constant supplied for unknown attribute "instancePickingColors"
        return super.getShaders({ vs: triangleVs, fs: heatmapPresentationFS, modules: [project32, triangleContourUniforms] });
    }

    initializeState({ device }: LayerContext): void {
        this.setState({ model: this._getModel(device) });
    }

    _getModel(device: Device): Model {
        const { vertexCount, data } = this.props;

        return new Model(device, {
            ...this.getShaders(),
            id: this.props.id,
            attributes: data.attributes,
            bufferLayout: [
                { name: 'positions', format: 'float32x3' },
                { name: 'texCoords', format: 'float32x2' }
            ],
            topology: 'triangle-strip',
            vertexCount,
            
            // if we include this, we get 
            // """
            // drawing TriangleLayerContours({id: 'scatter_-#srslt0detail-react#-contour1-triangle-layer'}) to screen: 
            // Cannot read properties of undefined (reading 'setUniforms') 
            // TypeError: Cannot read properties of undefined (reading 'setUniforms')
            // """
            // indexBuffer: device.createBuffer(new Uint16Array([0, 1, 3, 1, 2, 3]))
        });
    }

    draw(): void {
        const { model } = this.state;
        const {
            intensity,
            threshold,
            aggregationMode,
            colorDomain,
            contourOpacity,
            contourFill,
            fillOpacity,
            colorTexture,
            maxTexture,
            weightsTexture,
        } = this.props;
        model.shaderInputs.setProps({
            triangle: {
                aggregationMode,
                colorDomain,
                colorTexture,
                contourFill,
                contourOpacity,
                fillOpacity,
                maxTexture,
                weightsTexture,
                intensity,
                threshold,
            },
        });
        model.draw(this.context.renderPass);
    }
}

// ---------- copied from deck.gl triangle-layer
