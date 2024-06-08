#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.

    // if mapping falls outside of the texmap, set pixel to white
    if (sp.p_uv.x > 1 || sp.p_uv.x < 0 || sp.p_uv.y > 1 || sp.p_uv.y < 0) {
      return Color::White;
    }

    // get level
    float level = get_level(sp);

    switch (sp.lsm) {
      case L_ZERO:
        // set level to 0
        level = 0;

        if (sp.psm == P_NEAREST) {
          return sample_nearest(sp.p_uv, level);
        } else if (sp.psm == P_LINEAR) {
          return sample_bilinear(sp.p_uv, level);
        }
      break;

      case L_NEAREST:
        // round to nearest mipmap level
        level = round(level);

        // check bounds
        level = (level < 0) ? 0 : level;                          // if level is less than 0, set to 0
        level = (level > mipmap.size()) ? mipmap.size() : level;  // if level is too large, set to max

        if (sp.psm == P_NEAREST) {
          return sample_nearest(sp.p_uv, level);
        } else if (sp.psm == P_LINEAR) {
          return sample_bilinear(sp.p_uv, level);
        }
      break;

      case L_LINEAR:
        // find nearest levels
        int level_f = floor(level);
        int level_c = ceil(level);

        // check bounds
        level_f = (level_f < 0) ? 0 : level_f;                          // if level is too small, set to 0
        level_f = (level_f > mipmap.size()) ? mipmap.size() : level_f;  // if level is too big, set to max
        level_c = (level_c < 0) ? 0 : level_c;                          // if level is too small, set to 0
        level_c = (level_c > mipmap.size()) ? mipmap.size() : level_c;  // if level is too big, set to max

        if (sp.psm == P_NEAREST) {
          return (sample_nearest(sp.p_uv, level_f) + sample_nearest(sp.p_uv, level_c)) * 0.5;
        } else if (sp.psm == P_LINEAR) {
          return (sample_bilinear(sp.p_uv, level_f) + sample_bilinear(sp.p_uv, level_c)) * 0.5;
        }
      break;
    }
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    Vector2D ptxuv = Vector2D(sp.p_dx_uv.x * (width - 1), sp.p_dx_uv.y * (height - 1));
    Vector2D ptyuv = Vector2D(sp.p_dy_uv.x * (width - 1), sp.p_dy_uv.y * (height - 1));
    float Level = std::max(ptxuv.norm(), ptyuv.norm());
    float level = log2(Level);
    return level;

    return 0;
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    int tex_x = round(uv.x * (mipmap[level].width - 1));
    int tex_y = round(uv.y * (mipmap[level].height - 1));
    Color color = mipmap[level].get_texel(tex_x, tex_y);
    return color;
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    float tex_x = uv.x * (mipmap[level].width - 1);
    float tex_y = uv.y * (mipmap[level].height - 1);
    float rem_x = tex_x - floor(tex_x);
    float rem_y = tex_y - floor(tex_y);

    Color color1 = mipmap[level].get_texel(floor(tex_x), floor(tex_y));
    Color color2 = mipmap[level].get_texel(floor(tex_x), ceil(tex_y));
    Color color3 = mipmap[level].get_texel(ceil(tex_x), ceil(tex_y));
    Color color4 = mipmap[level].get_texel(ceil(tex_x), floor(tex_y));

    // horizontal lerp
    Color horiz1 = color1 + rem_x * (color4 + (-1 * color1));
    Color horiz2 = color2 + rem_y * (color4 + (-1 * color2));

    // vertical lerp
    Color vert = horiz1 + rem_y * (horiz2 + (-1 * horiz1));

    return vert;
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
