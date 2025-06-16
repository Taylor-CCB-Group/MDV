#version 300 es
#define SHADER_NAME triangle-layer-fragment-contour-shader

precision highp float;
out vec4 fragColor;
uniform float opacity;
uniform sampler2D weightsTexture; //using name 'texture' caused indirect glsl compiler errors with version 300 es(?)
uniform sampler2D colorTexture;
uniform float aggregationMode;

in vec2 vTexCoords;
in float vIntensityMin;
in float vIntensityMax;

// todo more structured uniforms, array of structs...
uniform float contourOpacity;
// struct ContourProps {
//     vec2 increment;
//     float lineWidth;
//     vec3 color;
// };
// uniform ContourProps contourProps;

float smoothContour(float value) {
    // value = clamp(value * vIntensityMax, 0., 1.);
    float width = 1.5; //todo better control, more coherent maths
    float w = width * fwidth(value);
    float f = 0.5; // reciprocol bandwidth
    if (value < f) return 0.; //todo something better
    if (value > 2. * f) return 1.; //metaballs - should be controllable parameter with nice animation
    float wa = smoothstep(0., (w * f), mod(value * f, 1.)-0.5); //nb -0.5 is to center the contour
    wa = 1. - max(smoothstep(1.-w, 1., wa), smoothstep(w, 0., wa));
    // return 0.0;
    return smoothstep(0., 1., wa*0.5);
    // return contour(value);
}


vec4 getLinearColor(float value) {
  float factor = clamp(value * vIntensityMax, 0., 1.);
  vec4 color = texture(colorTexture, vec2(factor, 0.5));
  color.a *= min(value * vIntensityMin, 1.0);
  return color;
}


void main(void) {
  vec4 weights = texture(weightsTexture, vTexCoords);
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
  float c = smoothContour(weight) * contourOpacity;
  vec4 fullColor = texture(colorTexture, vec2(1.0, 0.5)); //should be a uniform
  //fullColor = vec4(1.);
  linearColor = mix(linearColor, fullColor, c);
  //   linearColor = c * vec4(1.);
  linearColor.a *= opacity;
  linearColor.a = max(linearColor.a, c*2.);
  fragColor = linearColor;
}