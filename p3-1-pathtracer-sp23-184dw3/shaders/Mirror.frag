#version 330


uniform vec3 u_cam_pos;

uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
  
    vec3 wo = u_cam_pos-vec3(v_position);
    wo /= length(wo);
    vec3 n = vec3(v_normal);
    
    vec3 wi = -wo - 2*dot(-wo, n)*n;
    wi /= length(wi);
    
    out_color = texture(u_texture_cubemap, wi);
  out_color.a = 1;
}
