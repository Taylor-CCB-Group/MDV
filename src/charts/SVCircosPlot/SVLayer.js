import { Layer } from '@deck.gl/core';
import { Model, Geometry } from '@luma.gl/engine';
import { project,picking } from '@deck.gl/core';



const vs = `#version 300 es
precision highp float;
in vec3 instanceColor;
in float instanceFilter;
in vec3 instancePosition; // [start, length] in 0-100 units
in float t;
in vec3 instancePickingColors;
uniform float radius;
uniform float rotation;
uniform vec2 center;
uniform float totalLength;
out vec3 vColor;

vec2 getCirclePoint(float angle, float r, vec2 c) {
  float rad = radians(angle);
  return c + vec2(cos(rad), sin(rad)) * r;
}

void main(void) {
  if (instanceFilter != 0.0){
    return;
  }
  picking_setPickingColor(instancePickingColors);
  float startAngle = instancePosition.x * 360.0 / totalLength - 90.0 + rotation;
  float arcLength = instancePosition.y * 360.0 / totalLength;
  vec2 pos;
  float offset;
  //any insertions or unknown svs (z=3 or 5) are drawn on the outer edge
  if (instancePosition.z > 2.5){
    offset = 0.0;
  }
  else {
    offset = (3.0-instancePosition.z)*20.0;
  }
  // for translocation , inversions and deletions greater than 10 000 draw bezier curves
  // linking start and end of sv
  if (instancePosition.z < 3.0 && instancePosition.y > 10000.0) {
    float endAngle = startAngle + arcLength;
    // Endpoints on the circle
    vec2 p0 = getCirclePoint(startAngle, radius-offset+20.0, center);
    vec2 p2 = getCirclePoint(endAngle, radius-offset+20.0, center);
    // Control point: inward from the midpoint
    float midAngle = (startAngle + endAngle) / 2.0;
    float controlRadius = radius * 0.7; // inward, adjust as needed
    vec2 p1 = getCirclePoint(midAngle, controlRadius, center);
    // Quadratic Bezier interpolation
    pos = mix(mix(p0, p1, t), mix(p1, p2, t), t);
  }

  //draw as straight line if insertion or ends of SV are close together
  else {
      pos = getCirclePoint(startAngle, radius + (t*20.0)-offset, center);
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

class SVLayer extends Layer {
  static layerName = 'SVLayer';

  initializeState() {
    this.getAttributeManager().addInstanced({
      instancePosition: {
        size: 3,
        type: "uint32",
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
  }



    updateState(params) {
  
        super.updateState(params);
        if (params.changeFlags.dataChanged) {
          if(params.changeFlags.dataChanged == "init"){
            this.state.model = this._getModel();
          }
        }

        else{
            const trig = params.changeFlags.updateTriggersChanged;
        
          if (trig){
            const instance = trig.getColor?"instanceColor":"instanceFilter";
            const am = this.getAttributeManager();
            const attribute = am.attributes[instance];
            attribute.setData(attribute.state.binaryValue);
        }
        
    }
  }

  draw({uniforms}) {
    const { model } = this.state;
    if (model) {
      model.setUniforms({
        ...uniforms,
        radius: this.props.radius || 150.0,
        center: this.props.center || [0.0, 0.0],
        rotation: this.props.rotation || 0.0,
        totalLength: this.props.totalLength|| 3088269832.0
      });
      model.draw(this.context.renderPass);
    }
  }

  getNumInstances() {
    return 400000;
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
      vs,
      fs,
      geometry,
      modules: [project,picking],
      isInstanced: true,
      bufferLayout: this.getAttributeManager().getBufferLayouts()
    });
  }
}

export default SVLayer;