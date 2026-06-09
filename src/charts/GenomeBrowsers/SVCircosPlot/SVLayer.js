import { project, picking } from "@deck.gl/core";
import { Layer } from "@deck.gl/core";
import { Model, Geometry } from "@luma.gl/engine";

/* i could use deck.gl's filtering extension but it 
 it appears to work in a similar way shader-modules.ts line 241
    if (dataFilter_value == 0.0) {
       gl_Position = vec4(0.);
    }
 at the end of the vertex shadwer 
 so if we do it at the beginning and return early
 it should be more efficient
 keeps the compiler happy
*/
const vs = `#version 300 es
precision highp float;
in vec3 instanceColor;
in float instanceFilter;
in vec3 instancePosition; // [start, length] in 0-100 units
in float t;
in vec3 instancePickingColors;
uniform svLayerUniforms {
  float radius;
  float rotation;
  vec2 center;
  float totalLength;
} svLayer;
out vec3 vColor;

vec2 getCirclePoint(float angle, float r, vec2 c) {
  float rad = radians(angle);
  return c + vec2(cos(rad), sin(rad)) * r;
}

void main(void) {
  if (instanceFilter != 0.0){
    gl_Position = vec4(0.0);
    return;
  }
  picking_setPickingColor(instancePickingColors);
  float startAngle = instancePosition.x * 360.0 / svLayer.totalLength - 90.0 + svLayer.rotation;
  float arcLength = instancePosition.y * 360.0 / svLayer.totalLength;
  vec2 pos;
  // DEL (inner), DUP (middle), INS (outer) as separate non-overlapping concentric rings.
  // DEL, DUP, INS, UNKNOWN/ND as concentric rings, all thick and spaced out
  if (instancePosition.z >= 2.0 && instancePosition.z <= 5.0) {
    float ringBase;
    float thickness = 28.0;
    if (abs(instancePosition.z - 2.0) < 0.1) {
      ringBase = svLayer.radius - 80.0; // DEL innermost
    } else if (abs(instancePosition.z - 4.0) < 0.1) {
      ringBase = svLayer.radius - 54.0; // DUP
    } else if (abs(instancePosition.z - 3.0) < 0.1) {
      ringBase = svLayer.radius - 28.0; // INS
    } else {
      ringBase = svLayer.radius - 2.0; // UNKNOWN/ND outermost
    }
    pos = getCirclePoint(startAngle, ringBase + (t * thickness), svLayer.center);
  } else {
    // Inversions (z=1) and TRA/BND (z=0) are now further inward to avoid overlap
    float invBase, traBase;
    float thickness = 20.0;
    if (abs(instancePosition.z - 1.0) < 0.1) {
      invBase = svLayer.radius - 108.0;
      float traBase = svLayer.radius - 160.0;
      if (abs(instancePosition.y) > 10000.0) {
        float endAngle = startAngle + arcLength;
        vec2 p0 = getCirclePoint(startAngle, invBase + 10.0, svLayer.center);
        vec2 p2 = getCirclePoint(endAngle, invBase + 10.0, svLayer.center);
        float midAngle = (startAngle + endAngle) / 2.0;
        // Control point: not below TRA/BND ring, and slightly taller
        float controlRadius = invBase - 48.0;
        vec2 p1 = getCirclePoint(midAngle, controlRadius, svLayer.center);
        pos = mix(mix(p0, p1, t), mix(p1, p2, t), t);
      } else {
        pos = getCirclePoint(startAngle, invBase + (t * thickness), svLayer.center);
      }
    } else if (abs(instancePosition.z - 0.0) < 0.1) {
      traBase = svLayer.radius - 140.0;
      if (abs(instancePosition.y) > 10000.0) {
        float endAngle = startAngle + arcLength;
        vec2 p0 = getCirclePoint(startAngle, traBase + 10.0, svLayer.center);
        vec2 p2 = getCirclePoint(endAngle, traBase + 10.0, svLayer.center);
        float midAngle = (startAngle + endAngle) / 2.0;
        float controlRadius = traBase * 0.7;
        vec2 p1 = getCirclePoint(midAngle, controlRadius, svLayer.center);
        pos = mix(mix(p0, p1, t), mix(p1, p2, t), t);
      } else {
        pos = getCirclePoint(startAngle, traBase + (t * thickness), svLayer.center);
      }
    } else {
      pos = getCirclePoint(startAngle, svLayer.radius, svLayer.center);
    }
  }
  
    //draw as an arc problem is that for short svs the arc is very small
    // plus they all overlap
    //else{
    //    float angle = startAngle + t * arcLength;
    //    pos = getCirclePoint(angle, radius-offset, center);
    //}

  pos= project_position(pos);
  gl_Position = project_common_position_to_clipspace(vec4(pos, 0, 1.0));
  vColor = instanceColor / 255.0;
}
`;

const fs = `#version 300 es
precision highp float;
in vec3 vColor;
out vec4 fragColor;
void main(void) {
    //highlight color is not working
    fragColor = picking_filterHighlightColor(vec4(vColor, 1.0));
    fragColor = picking_filterPickingColor(fragColor);
}
`;

const svLayerModule = {
  name: 'svLayer',
  uniformTypes: {
    radius: 'f32',
    rotation: 'f32',
    center: 'vec2<f32>',
    totalLength: 'f32'
  },
  defaultUniforms: {
    radius: 150.0,
    rotation: 0.0,
    center: [0.0, 0.0],
    totalLength: 3088269832.0
  }
};

class SVLayer extends Layer {
  static layerName = 'SVLayer';

  initializeState() {
    this.getAttributeManager().addInstanced({
      instancePosition: {
        size: 3,
        accessor: 'getSV'
      },
      instanceFilter:{
        size:1,
        type:"uint8",
        accessor:"getFilter"
      },
    
      instanceColor: {
        size: 3,
        type: "uint8",
        accessor: 'getColor'
      }
    });

    this.setState({ model: this._getModel() });
    this._updateSVUniforms();
  }

  getShaders() {
    return {
      vs,
      fs,
      modules: [project, picking, svLayerModule]
    };
  }

  updateState(params) {
    super.updateState(params);

    const { changeFlags } = params;
    if (!this.state.model) {
      this.setState({ model: this._getModel() });
    }

    if (changeFlags.propsChanged || changeFlags.dataChanged) {
      this._updateSVUniforms();
    }

    const trig = changeFlags.updateTriggersChanged;
    if (trig) {
      const instance = trig.getColor ? 'instanceColor' : 'instanceFilter';
      const am = this.getAttributeManager();
      const attribute = am.attributes[instance];
      attribute.updateSubBuffer();
    }
  }

  finalizeState() {
    const { model } = this.state;
    if (model) {
      model.destroy();
    }
    super.finalizeState();
  }

  draw() {
    const { model } = this.state;
    if (model) {
      model.draw(this.context.renderPass);
    }
  }

  getNumInstances() {
    const sv = this.props.data?.attributes?.getSV;
    if (sv?.length) return sv.length / 3;
    return 0;
  }

  _getModel() {
    const NUM_SEGMENTS = 12;
    const geometry = new Geometry({
      topology: 'line-strip',
      vertexCount: NUM_SEGMENTS,
      isInstanced:true,
      attributes: {
          t: { size: 1, value: new Float32Array([...Array(NUM_SEGMENTS).keys()].map(i => i / (NUM_SEGMENTS - 1))) }
      },
    });

    return new Model(this.context.device, {
      ...this.getShaders(),
      geometry,
      isInstanced: true,
      bufferLayout: this.getAttributeManager().getBufferLayouts()
    });
  }

  _updateSVUniforms() {
    this.setShaderModuleProps({
      svLayer: {
        radius: this.props.radius || 150.0,
        center: this.props.center || [0.0, 0.0],
        rotation: this.props.rotation || 0.0,
        totalLength: this.props.totalLength || 3088269832.0
      }
    });
  }
}

export default SVLayer;
