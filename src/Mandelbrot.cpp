#include "Mandelbrot.h"
#include <cmath>
#include <complex>

Color Mandelbrot::getPixelColor(int x, int y) const {
  return pixelColor(iterations(x, y));
}

int Mandelbrot::iterations(int x, int y) const {

  double dx = (m_right - m_left) / (m_imgWidth - 1);
  double dy = (m_bottom - m_top) / (m_imgHeight - 1);
  std::complex<double> c(-2.0 + x * dx, -1.25 + y * dy);
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

int Mandelbrot::imgWidth() { return m_imgWidth; }
int Mandelbrot::imgHeight() { return m_imgHeight; }
