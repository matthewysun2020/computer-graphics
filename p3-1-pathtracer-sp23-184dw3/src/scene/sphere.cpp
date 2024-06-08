#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {
    
    // TODO (Part 1.4):
    // Implement ray - sphere intersection test.
    // Return true if there are intersections and writing the
    // smaller of the two intersection times in t1 and the larger in t2.
    
    
    return true;
    
}

bool Sphere::has_intersection(const Ray &r) const {
    
    // TODO (Part 1.4):
    // Implement ray - sphere intersection.
    // Note that you might want to use the the Sphere::test helper here.
    
    double a = dot(r.d,r.d);
    double b = dot(2*(r.o - o), r.d);
    double c = dot(r.o - o, r.o - o) - r2;
    
    double det = b*b - 4 * a * c;
    if (det == 0) {
        double t = -b / 2 / a;
        if (t >= r.min_t && t <= r.max_t) {
            r.max_t = t;
            return true;
        }
    } else if (det > 0) {
        double t1 = (-b + sqrt(det)) / 2 / a;
        double t2 = (-b - sqrt(det)) / 2 / a;
        
        if (t1 >= r.min_t && t1 <= r.max_t) {
            r.max_t = min(r.max_t, t1);
        }
        if (t2 >= r.min_t && t2 <= r.max_t) {
            r.max_t = min(r.max_t, t2);
        }
        
        if ((t1 >= r.min_t && t1 <= r.max_t) || (t2 >= r.min_t && t2 <= r.max_t)) {
            return true;
        }
    }
    
    return false;
}

bool Sphere::intersect(const Ray &r, Intersection *i) const {
    
    // TODO (Part 1.4):
    // Implement ray - sphere intersection.
    // Note again that you might want to use the the Sphere::test helper here.
    // When an intersection takes place, the Intersection data should be updated
    // correspondingly.
    
    double a = dot(r.d,r.d);
    double b = dot(2*(r.o - o), r.d);
    double c = dot(r.o - o, r.o - o) - r2;
    
    double det = b*b - 4 * a * c;
    if (det == 0) {
        double t = -b / 2 / a;
        if (t >= r.min_t && t <= r.max_t) {
            r.max_t = t;
            i->t = t;
            i->n = (r.o+t*r.d - o).unit();
            i->primitive = this;
            i->bsdf = get_bsdf();
            return true;
        }
    } else if (det > 0) {
        double t;
        double t1 = (-b + sqrt(det)) / 2 / a;
        double t2 = (-b - sqrt(det)) / 2 / a;
        
        if (t1 >= r.min_t && t1 <= r.max_t) {
            if (t2 >= r.min_t && t2 <= r.max_t) {
                t = min(t1, t2);
            } else {
                t = t1;
            }
        } else {
            if (t2 >= r.min_t && t2 <= r.max_t) {
                t = t2;
            } else {
                return false;
            }
        }
        r.max_t = t;
        i->t = t;
        i->n = (r.o+t*r.d - o).unit();
        i->primitive = this;
        i->bsdf = get_bsdf();
        return true;
    }
    
    return false;
}

void Sphere::draw(const Color &c, float alpha) const {
    Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
    // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
