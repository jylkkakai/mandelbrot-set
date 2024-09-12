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
  case Mode::avx:
    avx();
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

void Mandelbrot::avx() {

  int idx;
  for (int i = 0; i < m_diagWidth; i++) {
    for (int j = 0; j < m_diagHeight; j++) {

      idx = 0;
      double cr = m_left + i * m_dx;
      double ci = m_top + j * m_dy;
      __m128d c = _mm_setr_pd(cr, ci);
      __m128d z = c;
      __m128d zs = _mm_mul_pd(c, c);
      __m128d za = _mm_hadd_pd(zs, zs);
      __m128d zsqrt = _mm_sqrt_pd(za);

      while (idx < m_iterLimit && zsqrt[0] < 4.0) {

        __m128d z1 = _mm_hsub_pd(zs, zs);
        __m128d z2 = _mm_setr_pd(z[1], z[0]);
        z2 = _mm_mul_pd(z, z2);
        z2 = _mm_hadd_pd(z2, z2);
        z1[1] = z2[0];
        z = _mm_add_pd(z1, c);

        zs = _mm_mul_pd(z, z);
        za = _mm_hadd_pd(zs, zs);
        zsqrt = _mm_sqrt_pd(za);

        idx++;
      }
      m_colorArr[j][i] = pixelColor(idx);
    }
  }
}

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
    m_mode = Mode::avx;
    break;
  }
}
