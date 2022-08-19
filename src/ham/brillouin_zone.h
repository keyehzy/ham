#pragma once

#include <ham/geometry.h>

constexpr double pi = std::numbers::pi;

template <class LatticeKind>
class BrillouinZone final : public RectangleGrid {
 public:
  explicit BrillouinZone(const LatticeKind& lattice, int num = 50)
      : m_width((static_cast<double>(lattice.nx()) /
                 static_cast<double>(lattice.unitcell_size())) *
                lattice.lattice_vector_1().norm()),
        m_height(static_cast<double>(lattice.ny()) *
                 lattice.lattice_vector_2().norm()) {
    RectangleGrid::initialize(-pi / m_width, pi / m_width, -pi / m_height,
                              pi / m_height, num);
  };

 private:
  double m_width;
  double m_height;
};
