#include "iostream"
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with planes.
  // double prev_proj = dot(pm.last_position, normal) - dot(point, normal); // this is the length projected
  // double curr_proj = dot(pm.position, normal) - dot(point, normal);

  // if ((prev_proj * curr_proj) < 0) {
  //   Vector3D dir = (pm.position - pm.last_position).unit(); // may not be orthogonal to the plane
  //   double t = -prev_proj / dot(dir, normal);
  //   Vector3D isect = pm.last_position + t * dir;
  //   Vector3D crt = (isect + (SURFACE_OFFSET * normal)) - pm.last_position;
  //   pm.position = pm.last_position + (1 - friction) * crt;
  // } else if (curr_proj < 0) {
  //   Vector3D dir = (pm.position - pm.last_position).unit(); // may not be orthogonal to the plane
  //   double t = -prev_proj / dot(dir, normal);
  //   Vector3D isect = pm.last_position + t * dir;
  //   Vector3D crt = (isect + (SURFACE_OFFSET * normal * 2)) - pm.last_position;
  //   pm.position = pm.last_position + (1 - friction) * crt;
  // }
  double curr_proj = dot(pm.position - point, normal);
  if (curr_proj < 0) {
    Vector3D dir = (pm.position - pm.last_position).unit();
    double t = dot(point - pm.last_position, normal) / dot(dir, normal);
    Vector3D isect = pm.last_position + t * dir;
    Vector3D crt = (isect + (SURFACE_OFFSET * normal)) - pm.last_position;
    pm.position = pm.last_position + (1 - friction) * crt;
  }
}

void Plane::render(GLShader &shader) {
  nanogui::Color color(0.7f, 0.7f, 0.7f, 1.0f);

  Vector3f sPoint(point.x, point.y, point.z);
  Vector3f sNormal(normal.x, normal.y, normal.z);
  Vector3f sParallel(normal.y - normal.z, normal.z - normal.x,
                     normal.x - normal.y);
  sParallel.normalize();
  Vector3f sCross = sNormal.cross(sParallel);

  MatrixXf positions(3, 4);
  MatrixXf normals(3, 4);

  positions.col(0) << sPoint + 2 * (sCross + sParallel);
  positions.col(1) << sPoint + 2 * (sCross - sParallel);
  positions.col(2) << sPoint + 2 * (-sCross + sParallel);
  positions.col(3) << sPoint + 2 * (-sCross - sParallel);

  normals.col(0) << sNormal;
  normals.col(1) << sNormal;
  normals.col(2) << sNormal;
  normals.col(3) << sNormal;

  if (shader.uniform("u_color", false) != -1) {
    shader.setUniform("u_color", color);
  }
  shader.uploadAttrib("in_position", positions);
  if (shader.attrib("in_normal", false) != -1) {
    shader.uploadAttrib("in_normal", normals);
  }

  shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);
}
