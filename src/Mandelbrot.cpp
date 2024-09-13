#include "Mandelbrot.h"
#include <cmath>
#include <complex>
#include <immintrin.h>
#include <iostream>
#include <omp.h>
#include <raylib.h>

void Mandelbrot::draw() {

  switch (m_mode) {
  case Mode::simple:
    simple();
    break;
  case Mode::simple_opt:
    simple_opt();
    break;
  case Mode::openmp:
    openmp();
    break;
  case Mode::avx128:
    avx128();
    break;
  case Mode::avx256:
    avx256();
    break;
  default:
    simple();
  }
  drawColorArr();
}

void Mandelbrot::simple() {

  std::complex<double> c;
  std::complex<double> z;
  int idx;

  for (int i = 0; i < m_diagWidth; i++) {
    for (int j = 0; j < m_diagHeight; j++) {

      idx = 0;
      c = {m_left + i * m_dx, m_top + j * m_dy};
      z = c;

      while (idx < m_iterLimit && std::abs(z) < 4.0) {
        z = pow(z, 2) + c;
        idx++;
      }
      m_colorArr[j][i] = pixelColor(idx);
    }
  }
}

void Mandelbrot::simple_opt() {

  int idx;
  for (int i = 0; i < m_diagWidth; i++) {
    for (int j = 0; j < m_diagHeight; j++) {

      idx = 0;
      double cr = m_left + i * m_dx;
      double ci = m_top + j * m_dy;
      double zr = cr;
      double zi = ci;
      double zrs = zr * zr;
      double zis = zi * zi;
      double zsqrt = sqrt(zrs + zis);

      while (idx < m_iterLimit && zsqrt < 4.0) {

        double zrt = zrs - zis;
        double zit = 2 * zr * zi;
        zr = zrt + cr;
        zi = zit + ci;
        zrs = zr * zr;
        zis = zi * zi;

        zsqrt = sqrt(zrs + zis);
        idx++;
      }
      m_colorArr[j][i] = pixelColor(idx);
    }
  }
}

void Mandelbrot::openmp() {

#pragma omp parallel for num_threads(100)
  for (int i = 0; i < m_diagWidth; i++) {
    for (int j = 0; j < m_diagHeight; j++) {

      int idx = 0;
      double cr = m_left + i * m_dx;
      double ci = m_top + j * m_dy;
      double zr = cr;
      double zi = ci;
      double zrs = zr * zr;
      double zis = zi * zi;
      double zsqrt = sqrt(zrs + zis);

      while (idx < m_iterLimit && zsqrt < 4.0) {

        double zrt = zrs - zis;
        double zit = 2 * zr * zi;
        zr = zrt + cr;
        zi = zit + ci;
        zrs = zr * zr;
        zis = zi * zi;

        zsqrt = sqrt(zrs + zis);
        idx++;
      }
      m_colorArr[j][i] = pixelColor(idx);
    }
  }
}
void Mandelbrot::avx128() {

  int idx0;
  int idx1;
  __m128d c = _mm_setr_pd(m_left, m_top);
  __m128d a = _mm_setr_pd(m_dx, m_dy);

#pragma omp parallel for num_threads(100)
  for (auto i = 0; i < m_diagWidth; i++) {
    for (auto j = 0; j < m_diagHeight; j += 2) {

      idx0 = 0;
      idx1 = 0;

      // double cr = m_left + i * m_dx;
      // double ci = m_top + j * m_dy;
      __m128d b0 = _mm_setr_pd(i, j);
      __m128d b1 = _mm_setr_pd(i, j + 1);
      //
      __m128d c0 = _mm_fmadd_pd(a, b0, c);
      __m128d c1 = _mm_fmadd_pd(a, b1, c);
      // __m128d c0 = _mm_setr_pd(m_left + i * m_dx, m_top + j * m_dy);
      // __m128d c1 = _mm_setr_pd(m_left + i * m_dx, m_top + (j + 1) * m_dy);
      __m128d z0 = c0;
      __m128d z1 = c1;
      __m128d zs0 = _mm_mul_pd(c0, c0);
      __m128d zs1 = _mm_mul_pd(c1, c1);
      __m128d za = _mm_hadd_pd(zs0, zs1);
      __m128d zsqrt = _mm_sqrt_pd(za);

      while ((idx0 < m_iterLimit && idx1 < m_iterLimit) &&
             (zsqrt[0] < 4.0 || zsqrt[1] < 4.0)) {

        __m128d zr = _mm_hsub_pd(zs0, zs1);
        // __m128d zr1 = _mm_hsub_pd(zs1, zs1);
        __m128d zi0 = _mm_setr_pd(z0[1], z0[0]);
        __m128d zi1 = _mm_setr_pd(z1[1], z1[0]);
        zi0 = _mm_mul_pd(z0, zi0);
        // zi0 = _mm_hadd_pd(zi0, zi0);
        zi1 = _mm_mul_pd(z1, zi1);
        __m128d zi = _mm_hadd_pd(zi0, zi1);

        z0 = _mm_setr_pd(zr[0], zi[0]);
        z1 = _mm_setr_pd(zr[1], zi[1]);

        // zr[1] = zi[0];
        z0 = _mm_add_pd(z0, c0);
        z1 = _mm_add_pd(z1, c1);

        zs0 = _mm_mul_pd(z0, z0);
        zs1 = _mm_mul_pd(z1, z1);
        // zs = _mm_mul_pd(z, z);
        za = _mm_hadd_pd(zs0, zs1);
        if (zsqrt[0] < 4)
          idx0++;
        if (zsqrt[1] < 4)
          idx1++;
        zsqrt = _mm_sqrt_pd(za);
      }
      m_colorArr[j][i] = pixelColor(idx0);
      m_colorArr[j + 1][i] = pixelColor(idx1);
      // if (idx1 > 0 && i > 200)
      //   std::cout << idx1 << std::endl;
    }
  }
}

void Mandelbrot::avx256() {

  int idx0;
  int idx1;
  int idx2;
  int idx3;
  __m256d c = _mm256_setr_pd(m_left, m_top, m_left, m_top);
  __m256d a = _mm256_setr_pd(m_dx, m_dy, m_dx, m_dy);

#pragma omp parallel for num_threads(100)
  for (auto i = 0; i < m_diagWidth; i++) {
    for (auto j = 0; j < m_diagHeight; j += 4) {

      idx0 = 0;
      idx1 = 0;
      idx2 = 0;
      idx3 = 0;

      __m256d b0 = _mm256_setr_pd(i, j, i, j + 1);
      __m256d b1 = _mm256_setr_pd(i, j + 2, i, j + 3);
      //
      __m256d c0 = _mm256_fmadd_pd(a, b0, c);
      __m256d c1 = _mm256_fmadd_pd(a, b1, c);
      __m256d z0 = c0;
      __m256d z1 = c1;
      __m256d zs0 = _mm256_mul_pd(c0, c0);
      __m256d zs1 = _mm256_mul_pd(c1, c1);
      __m256d za = _mm256_hadd_pd(zs0, zs1);
      __m256d zsqrt = _mm256_sqrt_pd(za);

      while ((idx0 < m_iterLimit && idx1 < m_iterLimit && idx2 < m_iterLimit &&
              idx3 < m_iterLimit) &&
             (zsqrt[0] < 4.0 || zsqrt[1] < 4.0 || zsqrt[2] < 4.0 ||
              zsqrt[3] < 4.0)) {

        __m256d zr = _mm256_hsub_pd(zs0, zs1);
        __m256d zi0 = _mm256_setr_pd(z0[1], z0[0], z0[3], z0[2]);
        __m256d zi1 = _mm256_setr_pd(z1[1], z1[0], z1[3], z1[2]);
        zi0 = _mm256_mul_pd(z0, zi0);
        zi1 = _mm256_mul_pd(z1, zi1);
        __m256d zi = _mm256_hadd_pd(zi0, zi1);

        z0 = _mm256_setr_pd(zr[0], zi[0], zr[1], zi[1]);
        z1 = _mm256_setr_pd(zr[2], zi[2], zr[3], zi[3]);

        z0 = _mm256_add_pd(z0, c0);
        z1 = _mm256_add_pd(z1, c1);

        zs0 = _mm256_mul_pd(z0, z0);
        zs1 = _mm256_mul_pd(z1, z1);
        za = _mm256_hadd_pd(zs0, zs1);
        if (zsqrt[0] < 4)
          idx0++;
        if (zsqrt[1] < 4)
          idx1++;
        if (zsqrt[2] < 4)
          idx2++;
        if (zsqrt[3] < 4)
          idx3++;
        zsqrt = _mm256_sqrt_pd(za);
      }
      m_colorArr[j][i] = pixelColor(idx0);
      m_colorArr[j + 1][i] = pixelColor(idx1);
      m_colorArr[j + 2][i] = pixelColor(idx2);
      m_colorArr[j + 3][i] = pixelColor(idx3);
    }
  }
}

//
// void Mandelbrot::avx256() {
//
//   int idx0;
//   int idx1;
//   int idx2;
//   int idx3;
//   __m256d left = _mm256_setr_pd(m_left, m_left, m_left, m_left);
//   __m256d top = _mm256_setr_pd(m_top, m_top, m_top, m_top);
//   __m256d dx = _mm256_setr_pd(m_dx, m_dx, m_dx, m_dx);
//   __m256d dy = _mm256_setr_pd(m_dy, m_dy, m_dy, m_dy);
//
//   for (auto i = 0; i < m_diagWidth; i++) {
//     for (auto j = 0; j < m_diagHeight; j += 4) {
//
//       idx0 = 0;
//       idx1 = 0;
//       idx2 = 0;
//       idx3 = 0;
//
//       __m256d x = _mm256_setr_pd(i, i, i, i);
//       __m256d y = _mm256_setr_pd(j, j + 1, j + 2, j + 3);
//       __m256d cr = _mm256_fmadd_pd(x, dx, left);
//       __m256d ci = _mm256_fmadd_pd(y, dy, top);
//       __m256d zr = cr;
//       __m256d zi = ci;
//       __m256d zrs = _mm256_mul_pd(zr, zr);
//       __m256d zis = _mm256_mul_pd(zi, zi);
//       __m256d za = _mm256_add_pd(zrs, zis);
//       __m256d zsqrt = _mm256_sqrt_pd(za);
//
//       while ((idx0 < m_iterLimit && idx1 < m_iterLimit && idx2 < m_iterLimit
//       &&
//               idx3 < m_iterLimit) &&
//              (zsqrt[0] < 4.0 || zsqrt[1] < 4.0 || zsqrt[2] < 4.0 ||
//               zsqrt[3] < 4.0)) {
//
//         __m256d zrn = _mm256_sub_pd(zrs, zis);
//         __m256d zin = _mm256_mul_pd(zr, zi);
//         zin = _mm256_add_pd(zin, zin);
//
//         zr = _mm256_add_pd(zrn, cr);
//         zi = _mm256_add_pd(zin, ci);
//         zrs = _mm256_mul_pd(zr, zr);
//         zis = _mm256_mul_pd(zi, zi);
//         // std::cout << zis0 << ' '
//         //           << zis[0] /*<< ' ' << d[2] << ' ' << d[3] */ <<
//         std::endl;
//         // if (zsqrt[0] < 4)
//         //   idx0++;
//         // if (zsqrt[1] < 4)
//         //   idx1++;
//         // if (zsqrt[2] < 4)
//         //   idx2++;
//         // if (zsqrt[3] < 4)
//         //   idx3++;
//         // idx0++;
//         za = _mm256_add_pd(zrs, zis);
//         zsqrt = _mm256_sqrt_pd(za);
//       }
//       m_colorArr[j][i] = pixelColor(idx0);
//       m_colorArr[j + 1][i] = pixelColor(idx1);
//       m_colorArr[j + 2][i] = pixelColor(idx2);
//       m_colorArr[j + 3][i] = pixelColor(idx3);
//     }
//   }
// }

void Mandelbrot::drawColorArr() {

  Image img = {.data = &m_colorArr,
               .width = m_diagWidth,
               .height = m_diagHeight,
               .mipmaps = 1,
               .format = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8};

  Texture2D texture = LoadTextureFromImage(img);
  DrawTexture(texture, 0, 0, WHITE);
}

Color Mandelbrot::pixelColor(int i) const {

  Color color = WHITE;
  if (i == m_iterLimit)
    color = BLACK;
  else if (i > 90)
    color = RED;
  else if (i > 70)
    color = PINK;
  else if (i > 50)
    color = ORANGE;
  else if (i > 30)
    color = YELLOW;
  else if (i > 20)
    color = GREEN;
  else if (i > 10)
    color = LIME;
  else if (i > 5)
    color = SKYBLUE;
  else if (i > 4)
    color = BLUE;
  else if (i > 3)
    color = DARKBLUE;
  else if (i > 2)
    color = PURPLE;
  else if (i > 1)
    color = DARKPURPLE;
  return color;
}

void Mandelbrot::moveDiagram() {
  Vector2 delta = GetMouseDelta();

  m_left -= delta.x * m_dx;
  m_top -= delta.y * m_dy;
  m_right -= delta.x * m_dx;
  m_bottom -= delta.y * m_dy;
}

void Mandelbrot::zoom(int direction) {

  Vector2 mousePos = GetMousePosition();
  m_left += mousePos.x * 0.1 * m_dx * direction;
  m_top += mousePos.y * 0.1 * m_dy * direction;
  m_right -= (m_diagWidth - mousePos.x) * 0.1 * m_dx * direction;
  m_bottom -= (m_diagHeight - mousePos.y) * 0.1 * m_dy * direction;

  m_dx = (m_right - m_left) / (m_diagWidth - 1);
  m_dy = (m_bottom - m_top) / (m_diagHeight - 1);
}

void Mandelbrot::setMode(int key) {

  switch (key) {

  case KEY_ONE:
    m_mode = Mode::simple;
    break;
  case KEY_TWO:
    m_mode = Mode::simple_opt;
    break;
  case KEY_THREE:
    m_mode = Mode::openmp;
    break;
  case KEY_FOUR:
    m_mode = Mode::avx128;
    break;
  case KEY_FIVE:
    m_mode = Mode::avx256;
    break;
  default:
    m_mode = Mode::simple;
    break;
  }
}
