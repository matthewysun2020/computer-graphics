#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // TODO (Part 1): Build a grid of masses and springs.
  for (int j = 0; j < this->num_height_points; j++) {
    for (int i = 0; i < this->num_width_points; i++) {
      double x = i * width / (num_width_points - 1);
      double y = j * height / (num_height_points - 1);
      Vector3D position;
      if (orientation == HORIZONTAL) {
        position = Vector3D(x, 1.0, y);
      } else if (orientation == VERTICAL) {
        double sample = (double(rand()) / RAND_MAX) / 500 - (1/1000);
        position = Vector3D(x, y, sample);
      }
      point_masses.emplace_back(position, false);
    }
  }

  // pin the pinata
  for (int k = 0; k < pinned.size(); k++) {
    int x = pinned[k][0];
    int y = pinned[k][1];

    PointMass *pin = &point_masses[num_width_points * y + x];
    pin->pinned = true;
  }
  
  for (int j = 0; j < this->num_height_points; j++) {
    for (int i = 0; i < this->num_width_points; i++) {
      PointMass *target = &point_masses[num_width_points * j + i];
      if (i > 0) {
        // we will create a structural spring here
        PointMass *left = &point_masses[num_width_points * j + i - 1];
        springs.emplace_back(left, target, STRUCTURAL);
      }
      if (j > 0) {
        // we will create structural spring here
        PointMass *top = &point_masses[num_width_points * (j-1) + i];
        springs.emplace_back(top, target, STRUCTURAL);
        if (i > 0) {
          PointMass *topleft = &point_masses[num_width_points * (j-1) + i - 1];
          springs.emplace_back(topleft, target, SHEARING);
        }
        if (i < num_width_points - 1) {
          PointMass *topright = &point_masses[num_width_points * (j-1) + i + 1];
          springs.emplace_back(topright, target, SHEARING);
        }
      }
      if (i > 1) {
        PointMass *twoleft = &point_masses[num_width_points * j + i - 2];
        springs.emplace_back(twoleft, target, BENDING);
      }
      if (j > 1) {
        PointMass *twotop = &point_masses[num_width_points * (j - 2) + i];
        springs.emplace_back(twotop, target, BENDING);
      }

    }
  }

}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.
  for (int i = 0; i < point_masses.size(); i++) {
    PointMass *target = &point_masses[i];
    for (int j = 0; j < external_accelerations.size(); j++) {
      // cout <<target->forces << "before" << endl << flush;
      target->forces += mass * external_accelerations[j];
      // cout << target->forces << "after" << endl << flush;
    }
  }
  
  for (auto &s : springs) {
    Vector3D pa = s.pm_a->position;
    Vector3D pb = s.pm_b->position;
    double rest = s.rest_length;

    double sc = cp->ks;
    if (s.spring_type == BENDING && cp->enable_bending_constraints) {
      sc *= 0.2;
      s.pm_a->forces += sc * ((pa - pb).norm() - rest) * (pb - pa).unit();
      s.pm_b->forces += sc * ((pa - pb).norm() - rest) * (pa - pb).unit();
    } else if (s.spring_type == SHEARING && cp->enable_shearing_constraints) {
      s.pm_a->forces += sc * ((pa - pb).norm() - rest) * (pb - pa).unit();
      s.pm_b->forces += sc * ((pa - pb).norm() - rest) * (pa - pb).unit();
    } else if (s.spring_type == STRUCTURAL && cp->enable_structural_constraints) {
      s.pm_a->forces += sc * ((pa - pb).norm() - rest) * (pb - pa).unit();
      s.pm_b->forces += sc * ((pa - pb).norm() - rest) * (pa - pb).unit();
    }
  }

  // TODO (Part 2): Use Verlet integration to compute new point mass positions
  for (auto &target : point_masses) {
    if (!target.pinned) {
      // cout << target.forces << ":forces" << endl << flush;
      Vector3D pos = (target.position + (1. - cp->damping / 100.) * (target.position - target.last_position) + (delta_t * delta_t) * (target.forces / mass));
      target.last_position = target.position;
      target.position = pos;
    }

    target.forces = Vector3D();
  }

  // TODO (Part 4): Handle self-collisions.
  build_spatial_map();
  for (auto &target : point_masses) {
    self_collide(target, simulation_steps);
  }

  // TODO (Part 3): Handle collisions with other primitives.
  for (auto &target : point_masses) {
    for (auto &prim : *collision_objects) {
      prim->collide(target);
    }
  }

  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].

  for (auto &s : springs) {
    Vector3D pa = s.pm_a->position;
    Vector3D pb = s.pm_b->position;
    double rest = s.rest_length;
    if ((pb - pa).norm() > (1.1 * rest)) {
      double diff = (pb - pa).norm() - 1.1 * rest;
      Vector3D stretch = (pb-pa).unit() * diff;
      if (!s.pm_a->pinned && !s.pm_b->pinned) {
        s.pm_a->position += stretch / 2;
        s.pm_b->position -= stretch / 2;
      } else if (s.pm_a->pinned && !s.pm_b->pinned) {
        s.pm_b->position -= stretch;
      } else if (!s.pm_a->pinned && s.pm_b->pinned) {
        s.pm_a->position += stretch;
      }
    }
  }

}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.
  for (auto &target : point_masses) {
    double hash = hash_position(target.position);
    if (map.find(hash) == map.end()) {
      map[hash] = new vector <PointMass *>();
    }
    map[hash]->push_back(&target);
  }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.
  double hash = hash_position(pm.position);
  Vector3D crt(0.);
  int count = 0;

  for (auto &target : *(map[hash])) {
    if (target != &pm) {
      Vector3D dir = target->position - pm.position;
      if (dir.norm() < 2 * thickness) {
        count++;
        crt += dir - 2 * thickness * dir.unit();
      }
    }
  }

  if (count) {
    crt /= (count * simulation_steps);
    pm.position += crt;
  }
}

double cantor(double x, double y) {
  return .5 * (x + y) * (x + y + 1) + y;
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.

  double w = 3 * width / num_width_points;
  double h = 3 * height / num_height_points;
  double t = max(w, h);
  
  double x = floor(pos.x / w);
  double y = floor(pos.y / h);
  double z = floor(pos.z / t);

  return cantor(cantor(x, y), z);
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
