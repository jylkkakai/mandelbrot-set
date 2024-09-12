#ifndef MANDELBROT_H
#define MANDELBROT_H
#include <raylib.h>

enum class Mode {

  simple,
  simple_opt,
  openmp,
  avx
};

class Mandelbrot {
public:
  void draw();
  void setMode(int key);
  void moveDiagram();
  void zoom(int direction);

private:
  const int m_iterLimit = 200;
  static const int m_diagWidth = 800;
  static const int m_diagHeight = 600;
  Mode m_mode = Mode::simple;
  double m_left = -2.0;
  double m_right = 1.0;
  double m_top = -1.25;
  double m_bottom = 1.25;
  double m_dx = (m_right - m_left) / (m_diagWidth - 1);
  double m_dy = (m_bottom - m_top) / (m_diagHeight - 1);
  Color m_colorArr[m_diagHeight][m_diagWidth];

  Color pixelColor(int i) const;
  void simple();
  void simple_opt();
  void openmp();
  void avx();
  void drawColorArr();
};

#endif // !MANDELBROT_H
