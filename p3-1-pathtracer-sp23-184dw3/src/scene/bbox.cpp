#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>
#include <vector>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {
    // TODO (Part 2.2):
    // Implement ray - bounding box intersection test
    // If the ray intersected the bouding box within the range given by
    // t0, t1, update t0 and t1 with the new intersection times.

    Vector3D mins = (min-r.o) / r.d;
    Vector3D maxs = (max-r.o) / r.d;
    
    if (mins.x > maxs.x)
        std::swap(mins.x, maxs.x);
    
    if (mins.y > maxs.y)
        std::swap(mins.y, maxs.y);
    
    if (mins.z > maxs.z)
        std::swap(mins.z, maxs.z);

    double min_t = std::max(std::max(mins.x, mins.y), mins.z);
    double max_t = std::min(std::min(maxs.x, maxs.y), maxs.z);
    
    return !(min_t > max_t || max_t < 0);
}

void BBox::draw(Color c, float alpha) const {
    
    glColor4f(c.r, c.g, c.b, alpha);
    
    // top
    glBegin(GL_LINE_STRIP);
    glVertex3d(max.x, max.y, max.z);
    glVertex3d(max.x, max.y, min.z);
    glVertex3d(min.x, max.y, min.z);
    glVertex3d(min.x, max.y, max.z);
    glVertex3d(max.x, max.y, max.z);
    glEnd();
    
    // bottom
    glBegin(GL_LINE_STRIP);
    glVertex3d(min.x, min.y, min.z);
    glVertex3d(min.x, min.y, max.z);
    glVertex3d(max.x, min.y, max.z);
    glVertex3d(max.x, min.y, min.z);
    glVertex3d(min.x, min.y, min.z);
    glEnd();
    
    // side
    glBegin(GL_LINES);
    glVertex3d(max.x, max.y, max.z);
    glVertex3d(max.x, min.y, max.z);
    glVertex3d(max.x, max.y, min.z);
    glVertex3d(max.x, min.y, min.z);
    glVertex3d(min.x, max.y, min.z);
    glVertex3d(min.x, min.y, min.z);
    glVertex3d(min.x, max.y, max.z);
    glVertex3d(min.x, min.y, max.z);
    glEnd();
    
}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
    return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
