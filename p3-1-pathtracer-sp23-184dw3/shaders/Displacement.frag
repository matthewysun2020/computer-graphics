#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
    return texture(u_texture_2, uv).x;
}

void main() {
    vec3 t = vec3(v_tangent);
    vec3 n = vec3(v_normal);
    vec3 b = cross(t, n);
    mat3 TBN = mat3(t, b, n);
    float dU = (h(vec2(v_uv.x + (1 / u_texture_2_size.x), v_uv.y)) - h(v_uv)) * u_height_scaling * u_normal_scaling;
    float dV = (h(vec2(v_uv.x, v_uv.y+(1 / u_texture_2_size.y))) - h(v_uv)) * u_height_scaling * u_normal_scaling;
    
    vec3 n0 = vec3(-dU, -dV, 1.0);
    vec3 nd = TBN * n0;
        
    vec3 illumination = u_light_intensity / pow(length(u_light_pos - vec3(v_position)), 2.0);
    
    vec3 l = normalize(u_light_pos - vec3(v_position));
    vec3 v = normalize(u_cam_pos - vec3(v_position));
    vec3 bisect = normalize(v + l);
    
    float ka = 0.1;
    vec3 Ia = vec3(1.0);
    float kd = 0.8;
    float ks = 0.5;
    
    vec3 p0 = ka * Ia;
    vec3 p1 = kd * illumination * max(0, dot(nd, l));
    vec3 p2 = ks * illumination * pow(max(0, dot(nd, bisect)), 100.0);
    
    out_color = vec4(p0 + p1 + p2, 1);
}

