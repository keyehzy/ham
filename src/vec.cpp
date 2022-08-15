#include <ham/vec.h>

Vec2 operator+(const Vec2& v1, const Vec2& v2) {
  return Vec2{v1.x + v2.x, v1.y + v2.y};
}
