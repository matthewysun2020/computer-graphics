#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  
  vec3 ka = vec3(0.5, 0.5, 0.5);
  vec3 kd = vec3(1., 1., 1.);
  vec3 ks = vec3(0.7, 0.7, 0.7);
  vec3 Ia = vec3(0.1, 0.1, 0.1);
  float p = 100;

  vec3 light = u_light_pos - v_position.xyz;
  vec3 l = normalize(light);
  vec3 camdir = normalize(u_cam_pos - v_position.xyz);
  vec3 h = normalize(camdir + l);
  float radius = length(light);
  float rad2 = radius * radius;
  vec3 color = ka * Ia + kd * u_light_intensity / rad2 * max(0., dot(v_normal.xyz, l)) + ks * u_light_intensity / rad2 * pow(max(0, dot(v_normal.xyz, h)), p);
  out_color = vec4(color, 0.0);
  out_color.a = 1.0;
}

