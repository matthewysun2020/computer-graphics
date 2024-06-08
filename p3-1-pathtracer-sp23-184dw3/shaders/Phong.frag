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
    vec3 illumination = u_light_intensity / pow(length(u_light_pos - vec3(v_position)), 2.0);
    
    vec3 l = (u_light_pos - vec3(v_position));
    l = l / length(l);
    
    vec3 v = u_cam_pos - vec3(v_position);
    v = v / length(v);
    
    vec3 bisect = v + l;
    bisect = bisect / length(bisect);
    
    float ka = 0.1;
    vec3 Ia = vec3(1);
    vec3 kd = vec3(u_color);
    float ks = 0.5;
    
    vec3 p0 = ka * Ia;
    vec3 p1 = kd * illumination * max(0, dot(vec3(v_normal), l));
    vec3 p2 = ks * illumination * pow(max(0, dot(vec3(v_normal), bisect)), 100.0);
    
    out_color = vec4(p0+p1+p2, 1);
}

