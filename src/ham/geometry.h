#pragma once

#include <ham/types.h>
#include <ham/vec.h>

#include <cassert>
#include <vector>

std::vector<double> linspace(double start, double stop, int num = 50);

class RectangleGrid {
 public:
  RectangleGrid(double width, double height, Vec2d origin = Vec2d{0, 0})
      : m_width(width), m_height(height), m_origin(origin), m_size(50) {
    m_xs = linspace(m_origin.x, m_origin.x + m_width, m_size);
    m_ys = linspace(m_origin.y, m_origin.y + m_height, m_size);
  };

  Vec2d operator()(std::size_t i, std::size_t j) const {
    return Vec2d{m_xs[i], m_ys[j]};
  }

  double width() const { return m_width; }
  double height() const { return m_height; }
  Vec2d origin() const { return m_origin; }
  std::size_t size() const { return m_size; }

 private:
  double m_width;
  double m_height;
  Vec2d m_origin;
  std::size_t m_size;
  std::vector<double> m_xs;
  std::vector<double> m_ys;
};

class GeometricPath {
 public:
  GeometricPath(const std::vector<Vec2d>& points)
      : m_points(points), m_size(50) {
    // We assume that the points are sorted in the order of preference of the
    // use i.e. the user meant that the path goes P1 -> P2 -> ..., etc.

    for (std::size_t i = 0; i < m_points.size() - 1; i++) {
      Vec2d origin = m_points[i];
      Vec2d displacement = m_points[i + 1] - m_points[i];

      std::vector<double> xs =
          linspace(origin.x, origin.x + displacement.x, m_size);
      std::vector<double> ys =
          linspace(origin.y, origin.y + displacement.y, m_size);

      m_xs.insert(m_xs.end(), xs.begin(), xs.end());
      m_ys.insert(m_ys.end(), ys.begin(), ys.end());
    }
  };

  const std::vector<Vec2d>& points() const { return m_points; }

  Vec2d point(std::size_t i) const { return m_points[i]; }

  std::size_t size() const { return m_size * (m_points.size() - 1); }

  Vec2d operator()(std::size_t i) const { return Vec2d{m_xs[i], m_ys[i]}; }

 private:
  std::vector<Vec2d> m_points;
  std::size_t m_size;
  std::vector<double> m_xs;
  std::vector<double> m_ys;
};
