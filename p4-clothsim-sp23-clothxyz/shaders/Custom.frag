#version 330

// (Every uniform is available here.)

uniform mat4 u_view_projection;
uniform mat4 u_model;

uniform float u_normal_scaling;
uniform float u_height_scaling;

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

// Feel free to add your own textures. If you need more than 4,
// you will need to modify the skeleton.
uniform sampler2D u_texture_1;
uniform sampler2D u_texture_2;
uniform sampler2D u_texture_3;
uniform sampler2D u_texture_4;
uniform sampler2D u_texture_6;
uniform sampler2D u_texture_7;

uniform vec2 u_texture_6_size;
// Environment map! Take a look at GLSL documentation to see how to
// sample from this.
uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  return texture(u_texture_6, uv).r;
}

void main() {
  // Your awesome shader here!
  vec3 t = normalize(v_tangent.xyz);
  vec3 n = normalize(v_normal.xyz);
  vec3 b = cross(n, t);
  mat3 tbn = mat3(t, b, n);

  float w = u_texture_6_size[0];
  float hi = u_texture_6_size[1];
  float k = u_normal_scaling * u_height_scaling;

  float dU = (h(v_uv + vec2(1. / w, 0.)) - h(v_uv)) * k;
  float dV = (h(v_uv + vec2(0., 1. / hi)) - h(v_uv)) * k;

  vec3 no = vec3(-dU, -dV, 1.);
  vec3 nd = tbn * normalize(no);
  
  vec3 ka = vec3(0.1, 0.1, 0.1);
  vec3 kd = vec3(0.5, 0.5, 0.5);
  vec3 ks = vec3(0.7, 0.7, 0.7);
  vec3 Ia = vec3(0.1, 0.1, 0.1);
  float p = 64;

  vec3 light = u_light_pos - v_position.xyz;
  vec3 l = normalize(light);
  vec3 camdir = normalize(u_cam_pos - v_position.xyz);
  vec3 h = normalize(camdir + l);
  float radius = length(light);
  float rad2 = radius * radius;
  vec3 color = ka * Ia + kd * u_light_intensity / rad2 * max(0., dot(nd, l)) + ks * u_light_intensity / rad2 * pow(max(0, dot(nd, h)), p) + texture(u_texture_7, v_uv).xyz;
  out_color = vec4(color, 0.0);
  out_color.a = 1.0;
  // out_color = texture(u_texture_7, v_uv);
}
