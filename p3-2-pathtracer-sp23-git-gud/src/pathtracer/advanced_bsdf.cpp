#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Implement MirrorBSDF
        // return Vector3D();
        reflect(wo, wi);
        *pdf = 1;
        return reflectance / abs_cos_theta(*wi);
    }

    void MirrorBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Mirror BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Project 3-2: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        
        return exp(-1.0 * pow(tan(acos(h.z)) / alpha, 2)) / (PI * pow(alpha, 2) * pow(h.z, 4)); //lmao
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.

        double cti = abs_cos_theta(wi);
        Vector3D espks = eta * eta + k * k;
        Vector3D tecti = 2 * eta * cti;
        double cti_squared = cti * cti;

        Vector3D Rs = (espks - tecti + cti_squared) / (espks + tecti + cti_squared);
        Vector3D Rp = (espks * cti_squared - tecti + 1) / (espks * cti_squared + tecti + 1);
        return (Rs + Rp) / 2.0;
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Implement microfacet model here.
        Vector3D h = (wo + wi).unit();
        Vector3D n(0, 0, 1);

        // potentially need to check for wo and wi coming from outside
        return (wo.z > 0 && wi.z > 0) ? F(wi) * G(wo, wi) * D(h) / (4 * dot(n, wi) * dot(n, wo)) : Vector3D();
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        Vector2D sample = sampler.get_sample();
        double r1 = sample.x;
        double r2 = sample.y;

        double theta_h = atan(sqrt(-pow(alpha, 2) * log(1.0 - r1)));
        double phi_h = 2 * PI * r2;
        Vector3D h(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));

        *wi = (2.0 * dot(h, wo)) * h - wo;
        wi->normalize();

        if ((*wi).z < 0) {
            *pdf = 0;
            return Vector3D();
        }
        
        double p_theta = (2 * sin(theta_h)) / (pow(alpha, 2) * pow(cos(theta_h), 3)) * exp(-pow(tan(theta_h), 2) / pow(alpha, 2));
        double p_phi = 1 / (2 * PI);
        double p_w = (p_theta * p_phi) / sin(theta_h);

        *pdf = p_w / (4 * (dot(*wi, h)));
        
        // *wi = cosineHemisphereSampler.get_sample(pdf);

        return MicrofacetBSDF::f(wo, *wi);
    }

    void MicrofacetBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Micofacet BSDF"))
        {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 1
        // Implement RefractionBSDF
        // return Vector3D();
        if (refract(wo, wi, ior)) {
            *pdf = 1.0;
            double eta = wo.z > 0 ? 1.0 / ior : ior;
            return transmittance / abs_cos_theta(*wi) / pow(eta, 2);
        }
        return Vector3D();
    }

    void RefractionBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        // compute Fresnel coefficient and use it as the probability of reflection
        // - Fundamentals of Computer Graphics page 305
        if (!refract(wo, wi, ior)) {
            reflect(wo, wi);
            *pdf = 1.0;
            return reflectance / abs_cos_theta(*wi);
        } else {
            double r0 = pow((1 - ior) / (1 + ior), 2);
            double schlick = r0 + (1 - r0) * pow((1 - abs(wo.z)), 5);
            
            if (coin_flip(schlick)) {
                reflect(wo, wi);
                *pdf = schlick;
                return schlick * reflectance / abs_cos_theta(*wi);
            } else {
                *pdf = 1 - schlick;
                refract(wo, wi, ior);
                double eta = wo.z > 0 ? 1.0 / ior : ior;
                return (1 - schlick) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
            }
        }

    }

    void GlassBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

        // TODO Project 3-2: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.
        Matrix3x3 reflector(-1., 0., 0.,
                            0., -1., 0.,
                            0., 0.,  1.);
        *wi = reflector * wo;
    }

    bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

        // TODO Project 3-2: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.
        double wo_z = wo.z;

        // assign eta
        double eta = wo_z > 0 ? 1.0 / ior : ior;
        double snell = 1 - pow(eta, 2) * (1 - pow(wo_z, 2));

        bool cond = snell < 0;

        if (!cond) {
            wi->x = -eta * wo.x;
            wi->y = -eta * wo.y;
            wi->z = wo.z > 0 ? -sqrt(snell) : sqrt(snell);
            wi->normalize();
            return true;
        }

        return false;
        
    }

} // namespace CGL
