#include "Mandelbrot.h"
#include <cmath>
#include <complex>
#include <iostream>

Color Mandelbrot::getPixelColor(int x, int y) const {
  return pixelColor(iterations(x, y));
}

int Mandelbrot::iterations(int x, int y) const {

  double dx = (m_right - m_left) / (m_diagWidth - 1);
  double dy = (m_bottom - m_top) / (m_diagHeight - 1);
  std::complex<double> c(m_left + x * dx, m_top + y * dy);
  std::complex<double> z;
  int m_iterLimit = 100;
  z = c;
  Color color = WHITE;
  int i = 0;
  while (i < m_iterLimit && std::abs(z) < 2.0) {
    z = pow(z, 2) + c;
    i++;
  }
  return i;
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

int Mandelbrot::diagWidth() { return m_diagWidth; }
int Mandelbrot::diagHeight() { return m_diagHeight; }
