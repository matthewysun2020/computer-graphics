#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    int pt_count = points.size();
    std::vector<Vector2D> lerp_pts;
    if (pt_count == 1) {
      return points;
    }
    
    for (int i = 0; i < pt_count - 1; i++) {
      lerp_pts.push_back((1.0 - t) * points[i] + t * points[i + 1]);
    }
    return lerp_pts;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    int pt_count = points.size();
    std::vector<Vector3D> lerp_pts;
    if (pt_count == 1) {
      return points;
    }

    for (int i = 0; i < pt_count - 1; i++) {
      lerp_pts.push_back((1.0 - t) * points[i] + t * points[i + 1]);
    }
    return lerp_pts;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    if (points.size() == 1) {
      return points[0];
    }
    return evaluate1D(evaluateStep(points, t), t);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    std::vector<Vector3D> bezier_pts;
    int pt_count = controlPoints.size();
    for (int i = 0; i < pt_count; i++) {
      bezier_pts.push_back(evaluate1D(controlPoints[i], u));
    }

    return evaluate1D(bezier_pts, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    Vector3D weighted_n(0, 0, 0);
    HalfedgeCIter he = halfedge();

    do {
      // skip boundary faces
      if (!he->face()->isBoundary()) {
        // get area of face
        Vector3D p0 = position;
        Vector3D p1 = he->next()->vertex()->position;
        Vector3D p2 = he->next()->next()->vertex()->position;

        double face_area = (cross(p1 - p0, p2 - p0).norm()) / 2;
        weighted_n += face_area * he->face()->normal();
      }

      he = he->twin()->next();
    } while (he != halfedge());

    return weighted_n.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // This method should flip the given edge and return an iterator to the flipped edge.
    
    // run only if the edge is not a boundary edge
    if (!e0->isBoundary()) {
      // definitions

      // inner halfedges
      HalfedgeIter he0 = e0->halfedge();
      HalfedgeIter he1 = he0->next();
      HalfedgeIter he2 = he1->next();
      HalfedgeIter he3 = he0->twin();
      HalfedgeIter he4 = he3->next();
      HalfedgeIter he5 = he4->next();

      // outer halfedges
      HalfedgeIter he6 = he1->twin();
      HalfedgeIter he7 = he2->twin();
      HalfedgeIter he8 = he4->twin();
      HalfedgeIter he9 = he5->twin();

      // vertices
      VertexIter v0 = he0->vertex();
      VertexIter v1 = he3->vertex();
      VertexIter v2 = he2->vertex();
      VertexIter v3 = he5->vertex();

      // edges
      EdgeIter e1 = he1->edge();
      EdgeIter e2 = he2->edge();
      EdgeIter e3 = he4->edge();
      EdgeIter e4 = he5->edge();

      // faces
      FaceIter f1 = he0->face();
      FaceIter f2 = he3->face();

      // update args

      // inner halfedges
      he0->setNeighbors(he1, he3, v2, e0, f1);
      he1->setNeighbors(he2, he9, v3, e4, f1);
      he2->setNeighbors(he0, he6, v1, e1, f1);
      he3->setNeighbors(he4, he0, v3, e0, f2);
      he4->setNeighbors(he5, he7, v2, e2, f2);
      he5->setNeighbors(he3, he8, v0, e3, f2);

      // outer halfedges
      he6->setNeighbors(he6->next(), he2, v2, e1, he6->face());
      he7->setNeighbors(he7->next(), he4, v0, e2, he7->face());
      he8->setNeighbors(he8->next(), he5, v3, e3, he8->face());
      he9->setNeighbors(he9->next(), he1, v1, e4, he9->face());

      // vertices
      v0->halfedge() = he5;
      v1->halfedge() = he2;
      v2->halfedge() = he4;
      v3->halfedge() = he1;

      // edges
      e0->halfedge() = he0;
      e1->halfedge() = he2;
      e2->halfedge() = he4;
      e3->halfedge() = he5;
      e4->halfedge() = he1;

      // faces
      f1->halfedge() = he0;
      f2->halfedge() = he3;
    }

    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    // run only if the edge is not a boundary edge
    if (!e0->isBoundary()) {
      // definitions

      // inner halfedges
      HalfedgeIter he0 = e0->halfedge();
      HalfedgeIter he1 = he0->next();
      HalfedgeIter he2 = he1->next();
      HalfedgeIter he3 = he0->twin();
      HalfedgeIter he4 = he3->next();
      HalfedgeIter he5 = he4->next();

      // outer halfedges
      HalfedgeIter he6 = he1->twin();
      HalfedgeIter he7 = he2->twin();
      HalfedgeIter he8 = he4->twin();
      HalfedgeIter he9 = he5->twin();

      // vertices
      VertexIter v0 = he0->vertex();
      VertexIter v1 = he3->vertex();
      VertexIter v2 = he2->vertex();
      VertexIter v3 = he5->vertex();

      // edges
      EdgeIter e1 = he1->edge();
      EdgeIter e2 = he2->edge();
      EdgeIter e3 = he4->edge();
      EdgeIter e4 = he5->edge();

      // faces
      FaceIter f1 = he0->face();
      FaceIter f2 = he3->face();

      // new stuff

      // new inner halfedges
      HalfedgeIter he10 = newHalfedge();
      HalfedgeIter he11 = newHalfedge();
      HalfedgeIter he12 = newHalfedge();
      HalfedgeIter he13 = newHalfedge();
      HalfedgeIter he14 = newHalfedge();
      HalfedgeIter he15 = newHalfedge();

      // new vertex (return this!)
      VertexIter v = newVertex();

      // new edges
      EdgeIter e5 = newEdge();
      EdgeIter e6 = newEdge();
      EdgeIter e7 = newEdge();

      // new faces
      FaceIter f3 = newFace();
      FaceIter f4 = newFace();
      
      // update args

      // halfedges
      he0->setNeighbors(he1, he3, v, e0, f1);
      he1->setNeighbors(he2, he6, v1, e1, f1);
      he2->setNeighbors(he0, he11, v2, e5, f1);
      he3->setNeighbors(he4, he0, v1, e0, f2);
      he4->setNeighbors(he5, he15, v, e7, f2);
      he5->setNeighbors(he3, he9, v3, e4, f2);
      he6->setNeighbors(he6->next(), he1, v2, e1, he6->face());
      he7->setNeighbors(he7->next(), he12, v0, e2, he7->face());
      he8->setNeighbors(he8->next(), he14, v3, e3, he8->face());
      he9->setNeighbors(he9->next(), he5, v1, e4, he9->face());
      he10->setNeighbors(he11, he13, v0, e6, f3);
      he11->setNeighbors(he12, he2, v, e5, f3);
      he12->setNeighbors(he10, he7, v2, e2, f3);
      he13->setNeighbors(he14, he10, v, e6, f4);
      he14->setNeighbors(he15, he8, v0, e3, f4);
      he15->setNeighbors(he13, he4, v3, e7, f4);

      // new vertex
      v->position = 0.5 * (v0->position + v1->position);
      v->isNew = 1;
      v->halfedge() = he0;

      // old vertices
      v0->halfedge() = he10;
      v1->halfedge() = he1;
      v2->halfedge() = he12;
      v3->halfedge() = he5;

      // edges
      e0->halfedge() = he0;
      e1->halfedge() = he1;
      e2->halfedge() = he12;
      e3->halfedge() = he14;
      e4->halfedge() = he5;
      e5->halfedge() = he2;
      e6->halfedge() = he10;
      e7->halfedge() = he4;

      // set isNew for edges touching the new vertex 
      e0->isNew = 0;
      e6->isNew = 0;
      e5->isNew = 1;
      e7->isNew = 1;

      // faces
      f1->halfedge() = he0;
      f2->halfedge() = he3;
      f3->halfedge() = he10;
      f4->halfedge() = he13;

      return v;
    }

    return VertexIter();
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      e->isNew = 0; // none of the eld edges will be new

      // inner halfedges
      HalfedgeIter he0 = e->halfedge();
      HalfedgeIter he1 = he0->next();
      HalfedgeIter he2 = he1->next();
      HalfedgeIter he3 = he0->twin();
      HalfedgeIter he4 = he3->next();
      HalfedgeIter he5 = he4->next();

      // vertices
      Vector3D a = he0->vertex()->position;
      Vector3D b = he3->vertex()->position;
      Vector3D c = he2->vertex()->position;
      Vector3D d = he5->vertex()->position;

      // compute and store the position of the vertices
      e->newPosition = 0.375 * (a + b) + 0.125 * (c + d);
    }

    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
      v->isNew = 0; // none of the old vertices will be new

      HalfedgeIter h = v->halfedge();
      Vector3D orig_sum(0,0,0);

      do {
        orig_sum += h->twin()->vertex()->position;
        h = h->twin()->next();
      } while (h != v->halfedge());

      float n = (float) v->degree();
      float u = (n == 3.0) ? 0.1875 : (0.375 / n);

      v->newPosition = (1.0 - (n * u)) * v->position + u * orig_sum;
    }
    
    // split old edges
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      VertexIter v1 = e->halfedge()->vertex();
      VertexIter v2 = e->halfedge()->twin()->vertex();

      if (!(v1->isNew || v2->isNew)) {
        VertexIter v = mesh.splitEdge(e);
        v->newPosition = e->newPosition;
      }
    }

    // flip the right edges
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      VertexIter v1 = e->halfedge()->vertex();
      VertexIter v2 = e->halfedge()->twin()->vertex();

      if (e->isNew && (v1->isNew + v2->isNew == 1)) {
        mesh.flipEdge(e);
      }
    }

    // set all vertex positions to new positions
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
      v->position = v->newPosition;
      v->isNew = 0;
    }
  }
}
