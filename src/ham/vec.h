#pragma once

template <class T>
struct Vec2 {
  T x;
  T y;

  Vec2<T> operator+(const Vec2<T>& v) const;
};

template <class T>
Vec2<T> Vec2<T>::operator+(const Vec2<T>& v) const {
  return Vec2<T>{x + v.x, y + v.y};
}
