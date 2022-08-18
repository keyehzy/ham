#pragma once

#include <cmath>

template <class T>
struct Vec2 {
  T x;
  T y;

  Vec2<T> operator+(const Vec2<T>& v) const;
  Vec2<T> operator-(const Vec2<T>& v) const;
  Vec2<T> operator*(T value) const;
  T dot(const Vec2<T>& v) const;
  T norm() const;
};

template <class T>
inline Vec2<T> Vec2<T>::operator+(const Vec2<T>& v) const {
  return Vec2<T>{x + v.x, y + v.y};
}

template <class T>
inline Vec2<T> Vec2<T>::operator-(const Vec2<T>& v) const {
  return Vec2<T>{x - v.x, y - v.y};
}

template <class T>
inline T Vec2<T>::dot(const Vec2<T>& v) const {
  return x * v.x + y * v.y;
}

template <class T>
inline Vec2<T> Vec2<T>::operator*(T value) const {
  return Vec2<T>{value * x, value * y};
}

template <class T>
inline T Vec2<T>::norm() const {
  return std::sqrt(x*x + y*y);
}

template <class T>
Vec2<T> operator*(T value, const Vec2<T>& v) {
  return Vec2<T>{value * v.x, value * v.y};
}
