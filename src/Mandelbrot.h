#ifndef MANDELBROT_H
#define MANDELBROT_H
#include <raylib.h>

class Mandelbrot {
public:
  Color getPixelColor(int x, int y) const;
  int imgWidth();
  int imgHeight();

private:
  const int m_iterLimit = 100;
  const int m_imgWidth = 800;
  const int m_imgHeight = 600;
  double m_left = -2.0;
  double m_right = 1.0;
  double m_top = -1.25;
  double m_bottom = 1.25;
  double m_dx = (m_right - m_left) / (m_imgWidth - 1);
  double m_dy = (m_bottom - m_top) / (m_imgHeight - 1);

  int iterations(int x, int y) const;
  Color pixelColor(int i) const;
};

#endif // !MANDELBROT_H
