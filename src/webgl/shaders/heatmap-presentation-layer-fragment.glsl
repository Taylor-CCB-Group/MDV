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
uniform float contourFill;
uniform float fillOpacity;
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
    // if (value > contourFill * f) return 1.; //metaballs - should be controllable parameter with nice animation //no bool type for uniforms in new luma.gl afaik
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
  // Calculate contour shape (0-1 based on weight position relative to contour lines)
  float contourShape = smoothContour(weight);
  // float f = 0.5; // reciprocol bandwidth
  // if (value > contourFill * f) return 1.; //metaballs - should be controllable parameter with nice animation //no bool type for uniforms in new luma.gl afaik
  // Determine fill strength: fillOpacity controls the overall strength of the fill contribution
  float fillStrength = (weight > contourFill) ? fillOpacity : 0.0;
  vec4 fullColor = texture(colorTexture, vec2(1.0, 0.5)); //should be a uniform
  
  // Contour contribution: contourOpacity controls the overall strength
  // Only show contour when contourShape > 0 (on a contour line)
  // Don't mix with linearColor to avoid gradient - show fullColor directly
  float c = contourShape * contourOpacity;
  vec4 contourColor = fullColor;
  // Apply contourOpacity and boost alpha for visibility (original logic used c*2.)
  contourColor.a = c * 2.0;
  
  // Fill contribution: show linearColor when weight exceeds threshold
  // fillOpacity controls the overall strength of the fill contribution
  vec4 fillColor = linearColor;
  fillColor.a *= fillStrength;
  fillColor.rgb *= fillStrength; // Apply fillOpacity to RGB for proper blending
  
  // Combine: fill is base, contour overlays on top
  // When contourShape is 0, contourColor.a is 0, so we only see fillColor
  // When contourShape > 0, contourColor overlays on fillColor
  fragColor.rgb = fillColor.rgb * (1.0 - contourColor.a) + contourColor.rgb * contourColor.a;
  fragColor.a = fillColor.a * (1.0 - contourColor.a) + contourColor.a;
}