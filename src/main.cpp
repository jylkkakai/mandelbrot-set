#include "Mandelbrot.h"
#include <chrono>
#include <format>
#include <iostream>
#include <raylib.h>

constexpr int screenWidth = 800;
constexpr int screenHeight = 600;

int main() {

  InitWindow(screenWidth, screenHeight, "Mandelbrot");

  Mandelbrot mandelbrot;

  while (!WindowShouldClose()) {

    BeginDrawing();
    ClearBackground(WHITE);

    auto start = std::chrono::high_resolution_clock::now();
    mandelbrot.draw();
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::string dur = std::format("{:.3f}", float(duration.count()) / 1000);
    // DrawText(std::to_string(float(duration.count()) / 1000).c_str(), 700, 20,
    //          20, RED);
    DrawText((dur + " ms").c_str(), 680, 20, 20, RED);
    DrawFPS(20, 20);

    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
      mandelbrot.moveDiagram();
    }
    float move = GetMouseWheelMove();
    if (move != 0) {
      mandelbrot.zoom(int(move));
    }

    int key = GetKeyPressed();
    if (key != 0)
      mandelbrot.setMode(key);

    EndDrawing();
  }

  CloseWindow();
}
