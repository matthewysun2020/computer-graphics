#ifndef CGL_RAYMARCH_H
#define CGL_RAYMARCH_H

#include "scene.h"
#include "aggregate.h"

#include <vector>

namespace CGL { namespace SceneObjects {

/**
 * Bounding Volume Hierarchy for fast Ray - Primitive intersection.
 * Note that the BVHAccel is an Aggregate (A Primitive itself) that contains
 * all the primitives it was built from. Therefore once a BVHAccel Aggregate
 * is created, the original input primitives can be ignored from the scene
 * during ray intersection tests as they are contained in the aggregate.
 */
class RayMarch {
 public:

    RayMarch () { }

  /**
   * Ray - Aggregate intersection 2.
   * Check if the given ray intersects with the aggregate (any primitive in
   * the aggregate). If so, the input intersection data is updated to contain
   * intersection information for the point of intersection. Note that the
   * intersected primitive entry in the intersection should be updated to
   * the actual primitive in the aggregate that the ray intersected with and
   * not the aggregate itself.
   * \param r ray to test intersection with
   * \param i address to store intersection info
   * \return true if the given ray intersects with the aggregate,
             false otherwise
   */
    
    static float sharpenScale(float x){
        for(int i = 0; i < 6; i++){
            x = pow(2.0, x) - 1;
        }
        return x;
    }
    static Vector3D trace(const Ray& r, Intersection* i) {
//        total_rays++;

        float totalDistance = 0.0;
        int steps;
        for (steps=0; steps < MaximumRaySteps; steps++) {
            Vector3D p = r.o + totalDistance * r.d;
            float distance = fieldDE(p);
            totalDistance += distance;
            if (distance < MinimumDistance) {
                float multip = (sharpenScale(1.0-float(steps)/float(MaximumRaySteps)));
//                cout << "REE: " << multip << "\n";
                return multip * Vector3D(1) ;
            }
        }
        return Vector3D(0);
    }
    
    static float fieldDE(Vector3D z) {
        Vector3D v = Vector3D(z.x, z.y, z.z);
        v.x = z.x - floor(z.x) - 0.5;
        v.y = z.y - floor(z.y) - 0.5;
        v.z = z.z - floor(z.z) - 0.5;
        return (float) v.norm() - 0.1;
    }
    
    static float sphereDE(Vector3D z) {
        return ((float) (z).norm() - 0.2);
    }
    
    static float slirpenskiDE(Vector3D z)
    {
        int Iterations = 10;
        float Scale = 2;
    
        Vector3D a1 = Vector3D(1,1,1);
        Vector3D a2 = Vector3D(-1,-1,1);
        Vector3D a3 = Vector3D(1,-1,-1);
        Vector3D a4 = Vector3D(-1,1,-1);
        Vector3D c;
        int n = 0;
        float dist, d;
        while (n < Iterations) {
             c = a1;
            dist = (z-a1).norm();
                d = (z-a2).norm();
            if (d < dist) {
                c = a2;
                dist=d;
            }
             d = (z-a3).norm();
            if (d < dist) {
                c = a3;
                dist=d;
            }
             d = (z-a4).norm();
            if (d < dist) {
                c = a4;
                dist=d;
            }
            z = Scale*z-c*(Scale-1.0);
            n++;
        }

        return (z).norm() * pow(Scale, float(-n));
    }
    
    static float mandelbulbDE(Vector3D pos) {
        int Iterations = 10;
        float Bailout = 2;
        float Power = 8;
        
        
        Vector3D z = pos;
        float dr = 1.0;
        float r = 0.0;
        for (int i = 0; i < Iterations ; i++) {
            r = z.norm();
            if (r>Bailout) break;
            
            // convert to polar coordinates
            float theta = acos(z.z/r);
            float phi = atan2(z.y,z.x);
            dr =  pow( r, Power-1.0)*Power*dr + 1.0;
            
            // scale and rotate the point
            float zr = pow( r,Power);
            theta = theta*Power;
            phi = phi*Power;
            
            // convert back to cartesian coordinates
            z = zr*Vector3D(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
            z+=pos;
        }
        return 0.5*log(r)*r/dr;
    }

//  static unsigned long long total_rays, total_isects;

private:
    static const int MaximumRaySteps = 100;
    static constexpr float MinimumDistance = 0.001;
};

} // namespace SceneObjects
} // namespace CGL

#endif // CGL_BVH_H
