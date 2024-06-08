#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)


    rgb_framebuffer_target[3 * (y * width + x)] = (unsigned char)(c.r * 255);
    rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char)(c.g * 255);
    rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char)(c.b * 255);
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
                                         float x1, float y1,
                                         float x2, float y2,
                                         Color color) {
    // create a unit vector in the z-direction
    Vector3D z(0, 0, 1);

    // build a triangle
    Vector3D v0(x0, y0, 0);
    Vector3D v1(x1, y1, 0);
    Vector3D v2(x2, y2, 0);

    // set enumeration order (clockwise)
    if (cross(((v1 + v2) / 2) - v0, v1 - v0).z < 0) {
        swap(v1, v2);
    }

    // build the edges
    Vector3D e0 = v0 - v1;
    Vector3D e1 = v1 - v2;
    Vector3D e2 = v2 - v0;

    // create normals that point into the triangle
    Vector3D n0 = cross(z, e0);
    Vector3D n1 = cross(z, e1);
    Vector3D n2 = cross(z, e2);

    // make box bounds for triangle
    float x_list[] = { x0, x1, x2 };
    float y_list[] = { y0, y1, y2 };
    
    // x-bounds
    int min_x = (int)floor(*std::min_element(x_list, x_list + 3));
    int max_x = (int)ceil(*std::max_element(x_list, x_list + 3));

    // y-bounds
    int min_y = (int)floor(*std::min_element(y_list, y_list + 3));
    int max_y = (int)ceil(*std::max_element(y_list, y_list + 3));

    /**
        // loop through pixels within the box bounds
        for (int x = min_x + 0.5; x < max_x; x++) {
          for (int y = min_y + 0.5; y < max_y; y++) {
            Vector3D p(x, y, 0);    // sample pixel

            // check if sample pixel is in the triangle
            if ((dot(p - v1, n0) >= 0) && (dot(p - v2, n1) >= 0) && (dot(p - v0, n2) >= 0)) {
              fill_pixel(x, y, color);
            }
          }
        }
    */

    // task 2: supersampling and anti-aliasing

    // get sample rate
    int sr = sqrt(sample_rate);

    // loop through pixels within the box bounds
    for (int x = min_x; x < max_x; x++) {
      for (int y = min_y; y < max_y; y++) {

        Vector3D p(x, y, 0);    // sample pixel
        int s = 0;              // linear index for supersampling

        for (int i = 0; i < sr; i++) {
          p.x = (float)x + ((float)i + 0.5) / (float)sr;

          for (int j = 0; j < sr; j++) {
            p.y = (float)y + ((float)j + 0.5) / (float)sr;

            // check if supersampled point is in the triangle
            if ((dot(p - v1, n0) >= 0) && (dot(p - v2, n1) >= 0) && (dot(p - v0, n2) >= 0)) {
              int sx = (int)floor(x);
              int sy = (int)floor(y);

              if (sx < 0 || sx >= width) continue;
              if (sy < 0 || sy >= height) continue;

              sample_buffer[sample_rate * (sy * width + sx) + s] = color;
            }
            s++;
          }
        }
      }
    }

  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                            float x1, float y1, Color c1,
                                                            float x2, float y2, Color c2)
  {
    // create a unit vector in the z-direction
    Vector3D z(0, 0, 1);

    // build a triangle
    Vector3D v0(x0, y0, 0);
    Vector3D v1(x1, y1, 0);
    Vector3D v2(x2, y2, 0);

    // set enumeration order (clockwise)
    if (cross(((v1 + v2) / 2) - v0, v1 - v0).z < 0) {
      swap(v1, v2);
    }

    // build the edges
    Vector3D e0 = v0 - v1;
    Vector3D e1 = v1 - v2;
    Vector3D e2 = v2 - v0;

    // create normals that point into the triangle
    Vector3D n0 = cross(z, e0);
    Vector3D n1 = cross(z, e1);
    Vector3D n2 = cross(z, e2);

    // make box bounds for triangle
    float x_list[] = { x0, x1, x2 };
    float y_list[] = { y0, y1, y2 };

    // x-bounds
    int min_x = (int)floor(*std::min_element(x_list, x_list + 3));
    int max_x = (int)ceil(*std::max_element(x_list, x_list + 3));

    // y-bounds
    int min_y = (int)floor(*std::min_element(y_list, y_list + 3));
    int max_y = (int)ceil(*std::max_element(y_list, y_list + 3));

    // get sample rate
    int sr = sqrt(sample_rate);

    // create a matrix for the computation of barycentric coefficients
    Matrix3x3 M(x0, x1, x2, y0, y1, y2, 1, 1, 1);
    M = M.inv();

    // loop through pixels within the box bounds
    for (int x = min_x; x < max_x; x++) {
      for (int y = min_y; y < max_y; y++) {

        Vector3D p(x, y, 1);    // sample pixel
        int s = 0;              // linear index for supersampling

        for (int i = 0; i < sr; i++) {
          p.x = (float)x + ((float)i + 0.5) / (float)sr;

          for (int j = 0; j < sr; j++) {
            p.y = (float)y + ((float)j + 0.5) / (float)sr;

            // check if supersampled point is in the triangle
            if ((dot(p - v1, n0) >= 0) && (dot(p - v2, n1) >= 0) && (dot(p - v0, n2) >= 0)) {
              // get barycentric coefficients
              Vector3D weights = M * p;
              weights.z = 1 - weights.x - weights.y;

              int sx = (int)floor(x);
              int sy = (int)floor(y);

              if (sx < 0 || sx >= width) continue;
              if (sy < 0 || sy >= height) continue;

              sample_buffer[sample_rate * (sy * width + sx) + s] = weights.x * c0 + weights.y * c1 + weights.z * c2;
            }
            s++;
          }
        }
      }
    }
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
      float x1, float y1, float u1, float v1,
      float x2, float y2, float u2, float v2,
      Texture& tex)
  {
      // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
      // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
      // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
    Vector3D z(0, 0, 1);

    // build a triangle (renamed for variable name conflict)
    Vector3D p0(x0, y0, 0);
    Vector3D p1(x1, y1, 0);
    Vector3D p2(x2, y2, 0);

    // set enumeration order (clockwise)
    if (cross(((p1 + p2) / 2) - p0, p1 - p0).z < 0) {
      swap(p1, p2);
    }

    // build the edges
    Vector3D e0 = p0 - p1;
    Vector3D e1 = p1 - p2;
    Vector3D e2 = p2 - p0;

    // create normals that point into the triangle
    Vector3D n0 = cross(z, e0);
    Vector3D n1 = cross(z, e1);
    Vector3D n2 = cross(z, e2);

    // make box bounds for triangle
    float x_list[] = { x0, x1, x2 };
    float y_list[] = { y0, y1, y2 };

    // x-bounds
    int min_x = (int)floor(*std::min_element(x_list, x_list + 3));
    int max_x = (int)ceil(*std::max_element(x_list, x_list + 3));

    // y-bounds
    int min_y = (int)floor(*std::min_element(y_list, y_list + 3));
    int max_y = (int)ceil(*std::max_element(y_list, y_list + 3));

    // get sample rate
    int sr = sqrt(sample_rate);

    // create a matrix for the computation of barycentric coefficients
    Matrix3x3 M(x0, x1, x2, y0, y1, y2, 1, 1, 1);
    M = M.inv();
    Vector3D u = Vector3D(u0, u1, u2);
    Vector3D v = Vector3D(v0, v1, v2);

    // create a SampleParams struct
    SampleParams SP;
    SP.lsm = lsm;
    SP.psm = psm; 

    // loop through pixels within the box bounds
    for (int x = min_x; x < max_x; x++) {
      for (int y = min_y; y < max_y; y++) {

        Vector3D p(x, y, 1);    // sample pixel
        int s = 0;              // linear index for supersampling

        for (int i = 0; i < sr; i++) {
          p.x = (float)x + ((float)i + 0.5) / (float)sr;

          for (int j = 0; j < sr; j++) {
            p.y = (float)y + ((float)j + 0.5) / (float)sr;

            // check if supersampled point is in the triangle
            if ((dot(p - p1, n0) >= 0) && (dot(p - p2, n1) >= 0) && (dot(p - p0, n2) >= 0)) {
              Vector3D wp = M * p;                            // barycentric coordinates of p = (x, y)
              Vector3D wpdx = M * (p + Vector3D(1, 0, 0));    // barycentric coordinates of p_dx = (x + 1, y)
              Vector3D wpdy = M * (p + Vector3D(0, 1, 0));    // barycentric coordinates of p_dy = (x, y + 1)

              // constrain coordinate sums to 1
              wp.z = 1 - wp.x - wp.y;
              wpdx.z = 1 - wpdx.x - wpdx.y;
              wpdy.z = 1 - wpdy.x - wpdy.y;

              // get the U, V coordinates using the barycentric coordinates
              // Set the target coordinate within SP object
              Vector2D puv = Vector2D(dot(wp, u), dot(wp, v));
              Vector2D puvdx = Vector2D(dot(wpdx, u), dot(wpdx, v));
              Vector2D puvdy = Vector2D(dot(wpdy, u), dot(wpdy, v));

              // Calculate differential between target and actual for mapping
              SP.p_uv = puv;
              SP.p_dx_uv = puvdx - puv;
              SP.p_dy_uv = puvdy - puv;

              int sx = (int)floor(x);
              int sy = (int)floor(y);

              if (sx < 0 || sx >= width) continue;
              if (sy < 0 || sy >= height) continue;

              sample_buffer[sample_rate * (sy * width + sx) + s] = tex.sample(SP);
            }
            s++;
          }
        }
      }
    }

  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
      // TODO: Task 2: You may want to update this function for supersampling support

      this->sample_rate = rate;


      this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
      size_t width, size_t height)
  {
      // TODO: Task 2: You may want to update this function for supersampling support

      this->width = width;
      this->height = height;
      this->rgb_framebuffer_target = rgb_framebuffer;


      this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
      std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
      std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
      // TODO: Task 2: You will likely want to update this function for supersampling support


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        float r = 0;
        float g = 0;
        float b = 0;

        for (int s = 0; s < sample_rate; s++) {
          r += sample_buffer[sample_rate * (y * width + x) + s].r;
          g += sample_buffer[sample_rate * (y * width + x) + s].g;
          b += sample_buffer[sample_rate * (y * width + x) + s].b;
        }

        // take the average of each color value
        Color avg_color = Color::Black;
        avg_color = r / sample_rate;
        avg_color.g = g / sample_rate;
        avg_color.b = b / sample_rate;

        fill_pixel(x, y, avg_color);
      }
    }
  }

  Rasterizer::~Rasterizer() { }


}// CGL
