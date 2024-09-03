#include <cmath>
#include <complex>
#include <raylib.h>

constexpr int width = 1600;
constexpr int height = 1200;

Color pixelColor(int x, int y) {

  double dx = 3.0 / (width - 1);
  double dy = 2.5 / (height - 1);
  std::complex<double> c(-2.0 + x * dx, -1.25 + y * dy);
  std::complex<double> z;
  int limit = 100;
  z = c;
  Color color = WHITE;
  int i = 0;
  while (i < limit && std::abs(z) < 2.0) {
    z = pow(z, 2) + c;
    i++;
  }
  if (i == limit)
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

int main() {

  InitWindow(width, height, "Mandlebrot");

  while (!WindowShouldClose()) {

    BeginDrawing();
    ClearBackground(BLACK);
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        DrawPixel(i, j, pixelColor(i, j));
      }
    }
    DrawFPS(20, 20);
    EndDrawing();
  }
}
