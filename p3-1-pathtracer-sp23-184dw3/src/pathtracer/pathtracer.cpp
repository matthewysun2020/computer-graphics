#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"

#include "scene/raymarch.h"



using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
    gridSampler = new UniformGridSampler2D();
    hemisphereSampler = new UniformHemisphereSampler3D();
    
    tm_gamma = 2.2f;
    tm_level = 1.0f;
    tm_key = 0.18;
    tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
    delete gridSampler;
    delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
    sampleBuffer.resize(width, height);
    sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
    bvh = NULL;
    scene = NULL;
    camera = NULL;
    sampleBuffer.clear();
    sampleCountBuffer.clear();
    sampleBuffer.resize(0, 0);
    sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
    sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
    // Estimate the lighting from this intersection coming directly from a light.
    // For this function, sample uniformly in a hemisphere.
    
    // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
    // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.
    
    // make a coordinate system for a hit point
    // with N aligned with the Z direction.
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();
    
    // w_out points towards the source of the ray (e.g.,
    // toward the camera if this is a primary ray)
    const Vector3D hit_p = r.o + r.d * isect.t;
    const Vector3D w_out = w2o * (-r.d);
    
    // This is the same number of total samples as
    // estimate_direct_lighting_importance (outside of delta lights). We keep the
    // same number of samples for clarity of comparison.
    int num_samples = scene->lights.size() * ns_area_light;
    Vector3D L_out;
    
    // TODO (Part 3): Write your sampling loop here
    // TODO BEFORE YOU BEGIN
    // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading
    
    Vector3D total = Vector3D();
    
    for (int i = 0; i < num_samples; i++) {
        Vector3D wi = hemisphereSampler->get_sample(); // in obj coords
        Vector3D f = isect.bsdf->f(w_out, wi);
        
        Vector3D wi_world = o2w * wi;
        Ray ray(hit_p, wi_world, 1);
        ray.min_t = EPS_F;
        
        Intersection intersection;
        if (bvh->intersect(ray, &intersection) && cos_theta(wi) > EPS_F) {
            total += intersection.bsdf->get_emission() * f * cos_theta(wi) / (1.0/PI);
            //            cout << intersection.bsdf->get_emission() << f << cos_theta(wi) << (1.0/PI) << "\n";
            //            cout <<intersection.bsdf->get_emission() * f * cos_theta(wi) / (1.0/PI) << "\n";
            //            if (!(n == Vector3D())) {
            //                cout << n << "\n";
            //            }
            //            total += n;
        }
    }
    
    //    cout << total << "\n";
    
    return total / (double) num_samples;
    
}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
    // Estimate the lighting from this intersection coming directly from a light.
    // To implement importance sampling, sample only from lights, not uniformly in
    // a hemisphere.
    
    // make a coordinate system for a hit point
    // with N aligned with the Z direction.
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();
    
    // w_out points towards the source of the ray (e.g.,
    // toward the camera if this is a primary ray)
    const Vector3D hit_p = r.o + r.d * isect.t;
    const Vector3D w_out = w2o * (-r.d);
    Vector3D L_out;
    
    for (auto light : scene->lights) {
        if (light->is_delta_light()) {
            double distToLight, pdf;
            Vector3D wi;
            Intersection i;
            Vector3D rad = light->sample_L(hit_p, &wi, &distToLight, &pdf);
            Ray shadow = Ray(hit_p, wi, 1);
            shadow.min_t = EPS_F;
            shadow.max_t = distToLight - EPS_F;
            if(!bvh->intersect(shadow, &i) && dot(wi, isect.n) > EPS_F)
                L_out += rad * isect.bsdf->f(w_out, wi) * dot(wi, isect.n) / pdf;
        } else {
            for (int j = 0; j < ns_area_light; j++) {
                double distToLight, pdf;
                Vector3D wi;
                Intersection i;
                Vector3D rad = light->sample_L(hit_p, &wi, &distToLight, &pdf);
                Ray shadow = Ray(hit_p, wi);
                shadow.min_t = EPS_F;
                shadow.max_t = distToLight - EPS_F;
                if(!bvh->intersect(shadow, &i) && dot(wi, isect.n) > EPS_F)
                    L_out += rad * isect.bsdf->f(w_out, wi) * dot(wi, isect.n) / pdf / ns_area_light;
            }
        }
    }
    
    return L_out;
    
}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
    // TODO: Part 3, Task 2
    // Returns the light that results from no bounces of light
    
    return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
    // TODO: Part 3, Task 3
    // Returns either the direct illumination by hemisphere or importance sampling
    // depending on `direct_hemisphere_sample`
    if (direct_hemisphere_sample)
        return estimate_direct_lighting_hemisphere(r, isect);
    
    return estimate_direct_lighting_importance(r, isect);
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();
    
    Vector3D hit_p = r.o + r.d * isect.t;
    Vector3D w_out = w2o * (-r.d);
    
    Vector3D L_out(0, 0, 0);
    
    // TODO: Part 4, Task 2
    // Returns the one bounce radiance + radiance from extra bounces at this point.
    // Should be called recursively to simulate extra bounces.
    if (r.depth < 1)
        return L_out;
    
    L_out = one_bounce_radiance(r, isect);
    
    double prob = 0.35;
    Vector3D wi;
    double pdf;
    Vector3D f = isect.bsdf->sample_f(w_out, &wi, &pdf);
    Vector3D new_dir = o2w * wi;
    Intersection i = Intersection(isect);
    double cos = cos_theta(wi);

    if (coin_flip(1.0-prob)) {
        Ray bounce = Ray(hit_p, new_dir, (int) (r.depth - 1));
        bounce.min_t = EPS_F;
        if (bvh->intersect(bounce, &i) && cos > EPS_F) {
            L_out += (at_least_one_bounce_radiance(bounce, i) * f * cos) / pdf / (1.0-prob);
//            cout << "im in.\n";
//            cout << at_least_one_bounce_radiance(bounce, isect) << "\n";
//            cout << bsdf << "\n";
//            cout << cos_theta(wi) << "\n";
//            cout << pdf << "\n";
//            cout << 1.0-prob << "\n";

//            L_out += diff;
        }
    }
    return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
    Intersection isect;
    Vector3D L_out;
    
//    return Vector3D(1);
    
    Vector3D color = RayMarch::trace(r, &isect);
//    cout << color;
    return color;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
    // TODO (Part 1.2):
    // Make a loop that generates num_samples camera rays and traces them
    // through the scene. Return the average Vector3D.
    // You should call est_radiance_global_illumination in this function.
    
    // TODO (Part 5):
    // Modify your implementation to include adaptive sampling.
    // Use the command line parameters "samplesPerBatch" and "maxTolerance"
    
    double num_samples = 0;          // total samples to evaluate
    Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel
    double xk = 0, xk2 = 0;
    
    Vector3D total = Vector3D(0);
    
    bool br = false;
    
    for(int i = 0; i < ns_aa && !br; i++) {
        if (i > 0 && i % samplesPerBatch == 0) {
            double var = (1 / (num_samples - 1)) * (xk2 - (xk * xk) / num_samples);
            if(1.96 * (sqrt(var) / sqrt(num_samples)) <= maxTolerance * (xk / num_samples))
                br = true;
        }
        
        Vector2D sample = gridSampler->get_sample();
        Ray ray = camera->generate_ray((x+sample.x)/sampleBuffer.w, (y+sample.y)/sampleBuffer.h);
        ray.depth = max_ray_depth;
        Vector3D color = est_radiance_global_illumination(ray);
        auto illum = color.illum();
        xk += illum;
        xk2 += illum * illum;
        total += color;
        num_samples++;
    }
    
    total /= num_samples;
    sampleBuffer.update_pixel(total, x, y);
    sampleCountBuffer[x + y * sampleBuffer.w] = num_samples;
}

void PathTracer::autofocus(Vector2D loc) {
    Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
    Intersection isect;
    
    bvh->intersect(r, &isect);
    
    camera->focalDistance = isect.t;
}

} // namespace CGL
