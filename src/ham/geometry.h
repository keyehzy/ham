#pragma once

#include <ham/types.h>
#include <ham/vec.h>

#include <cassert>
#include <vector>

std::vector<double> linspace(double start, double stop, int num = 50);

class RectangleGrid {
 private:
 public:
  RectangleGrid(double width, double height, Vec2<double> origin = Vec2d{0, 0})
      : m_width(width), m_height(height), m_origin(origin), m_size(50) {
    m_xs = linspace(m_origin.x, m_origin.x + m_width, m_size);
    m_ys = linspace(m_origin.y, m_origin.y + m_height, m_size);
  };

  Vec2<double> point(std::size_t i, std::size_t j) const {
    return Vec2<double>{m_xs[i], m_ys[j]};
  }

  double width() const { return m_width; }
  double height() const { return m_height; }
  std::size_t size() const { return m_size; }

 private:
  double m_width;
  double m_height;
  Vec2<double> m_origin;
  std::size_t m_size;
  std::vector<double> m_xs;
  std::vector<double> m_ys;
};
