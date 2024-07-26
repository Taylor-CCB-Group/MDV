#version 300 es
#define SHADER_NAME triangle-layer-fragment-contour-shader

precision highp float;

uniform float opacity;
uniform sampler2D weightTexture; //using name 'texture' caused indirect glsl compiler errors with version 300 es(?)
uniform sampler2D colorTexture;
uniform float aggregationMode;

varying vec2 vTexCoords;
varying float vIntensityMin;
varying float vIntensityMax;

vec4 getLinearColor(float value) {
  float factor = clamp(value * vIntensityMax, 0., 1.);
  vec4 color = texture2D(colorTexture, vec2(factor, 0.5));
  color.a *= min(value * vIntensityMin, 1.0);
  return color;
}

float contour(float value) {
    //placeholder pending getting surrounding code in better shape
    return mod(value, 1.0) < 0.1 ? 1.0 : 0.0;
}
float smoothContour(float value) {
    float width = 0.1; //todo better control
    float w = fwidth(value);
    float f = .5;
    if (value < f) return 0.; //todo something better
    float wa = smoothstep(0., w * f, mod(value * f, 1.));
    wa = 1. - max(smoothstep(1.-w, 1., wa), smoothstep(w, 0., wa));
    return smoothstep(0., 1., wa*0.5);
    // return contour(value);
}

void main(void) {
  vec4 weights = texture2D(weightTexture, vTexCoords);
  float weight = weights.r;

  if (aggregationMode > 0.5) {
    weight /= max(1.0, weights.a);
  }
  // discard pixels with 0 weight.
  if (weight <= 0.) {
    discard;
  }

  vec4 linearColor = getLinearColor(weight);
  // todo allow for multiple contours, with different properties
  float c = smoothContour(weight);
  linearColor = mix(linearColor, vec4(1.0, 1.0, 1.0, 1.0), c);
  //   linearColor = c * vec4(1.);
  linearColor.a *= opacity;
  linearColor.a = max(linearColor.a, c*2.);
  gl_FragColor =linearColor;
}