#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with spheres.--
  Vector3D dist = pm.position - origin;
  double diff = dist.norm() - radius;

  if (diff < 0) {
    Vector3D isect = origin + radius * dist.unit();
    Vector3D crt = isect - pm.last_position;
    pm.position = pm.last_position + (1 - friction) * crt;
  }
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}
