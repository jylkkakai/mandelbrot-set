#include "Mandelbrot.h"
#include <raylib.h>

constexpr int screenWidth = 800;
constexpr int screenHeight = 600;

int main() {

  InitWindow(screenWidth, screenHeight, "Mandelbrot");

  Mandelbrot mandelbrot;

  while (!WindowShouldClose()) {

    BeginDrawing();
    ClearBackground(BLACK);
    for (int i = 0; i < mandelbrot.diagWidth(); i++) {
      for (int j = 0; j < mandelbrot.diagHeight(); j++) {
        DrawPixel(i, j, mandelbrot.getPixelColor(i, j));
      }
    }
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
      mandelbrot.moveDiagram();
    }
    DrawFPS(20, 20);
    EndDrawing();
  }
}
