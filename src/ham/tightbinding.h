#pragma once

#include <ham/lattice.h>
#include <ham/matrix.h>

#include <functional>

using namespace std::complex_literals;

struct TightBindingParameters {
  double t = 1.0;
};

template <class LatticeKind>
class Tightbinding {
 public:
  explicit Tightbinding(int nx, int ny, TightBindingParameters parameters)
      : m_lattice(LatticeKind(nx, ny)), m_parameters(parameters){};

  const LatticeKind& lattice() const { return m_lattice; }

  Matrix<double> realspace_hamiltonian() const;

  using FactorFn = std::function<Complex(Vec2d, Vec2d)>;

  Matrix<Complex> momentum_hamiltonian_base(Vec2d k, FactorFn factor) const;

  Matrix<Complex> momentum_hamiltonian(Vec2d k) const;

  Matrix<Complex> momentum_hamiltonian_x_derivative(Vec2d k) const;

  Matrix<Complex> momentum_hamiltonian_y_derivative(Vec2d k) const;

  int size() const { return m_lattice.orbital_count(); }

 private:
  LatticeKind m_lattice;
  TightBindingParameters m_parameters;
};

template <class LatticeKind>
Matrix<double> Tightbinding<LatticeKind>::realspace_hamiltonian() const {
  Matrix<double> h(this->size(), this->size());

  // 1) For each site, go to its neighbors and for each neighbor set a value at
  // the matrix entry (site_index, neighbor_index).
  // 2) For multiple orbitals we need to pad the entries.

  for (int site_index = 0; site_index < m_lattice.site_count(); site_index++) {
    for (int neighbor = 0; neighbor < m_lattice.nearest_neighbors_size();
         neighbor++) {
      int neighbor_index = m_lattice.site(site_index).neighbors[neighbor].index;
      for (int orbital_index = 0; orbital_index < m_lattice.orbitals();
           orbital_index++) {
        h(site_index * m_lattice.orbitals() + orbital_index,
          neighbor_index * m_lattice.orbitals() + orbital_index) =
            m_parameters.t;
      }
    }
  }
  return h;
}

template <class LatticeKind>
Matrix<Complex> Tightbinding<LatticeKind>::momentum_hamiltonian_base(
    Vec2d k, FactorFn factor) const {
  // We can only perform this operation when we have closed periodic boundary
  // conditions, otherwise we cannot use Bloch theorem to justify using
  // tightbinding.
  assert(m_lattice.boundary() == Boundary::Closed);

  Matrix<Complex> h(this->size(), this->size());

  // Go through each site and check its neighbors, both the indices from which
  // the electron hop from and to are known as well as the vector used in the
  // fourier transform. It is important that we sum the matrix entry in case we
  // have a small unitcell where there are two hopping coming from sites with
  // same index (for instance wrapping around from the edges). This way we don't
  // override the previous hopping with the next one.

  // Here we need to make an extra step.
  // 3) For any nearest-neighbors hamiltonian, pass a function "factor" which
  // models the hamiltonian.

  for (int site_index = 0; site_index < m_lattice.site_count(); site_index++) {
    for (int neighbor = 0; neighbor < m_lattice.nearest_neighbors_size();
         neighbor++) {
      Edge neighbor_edge = m_lattice.site(site_index).neighbors[neighbor];
      Vec2d delta = m_lattice.delta(neighbor_edge.direction);
      for (int orbital_index = 0; orbital_index < m_lattice.orbitals();
           orbital_index++) {
        h(site_index * m_lattice.orbitals() + orbital_index,
          neighbor_edge.index * m_lattice.orbitals() + orbital_index) +=
            m_parameters.t * factor(k, delta);
      }
    }
  }
  return h;
}

template <class LatticeKind>
Matrix<Complex> Tightbinding<LatticeKind>::momentum_hamiltonian(Vec2d k) const {
  return this->momentum_hamiltonian_base(k, [](Vec2d k, Vec2d delta) {
    double phase = k.dot(delta);
    return std::exp(-1.0i * phase);
  });
}

template <class LatticeKind>
Matrix<Complex> Tightbinding<LatticeKind>::momentum_hamiltonian_x_derivative(
    Vec2d k) const {
  return this->momentum_hamiltonian_base(k, [](Vec2d k, Vec2d delta) {
    double phase = k.dot(delta);
    return (-1.0i * delta.x) * std::exp(-1.0i * phase);
  });
}

template <class LatticeKind>
Matrix<Complex> Tightbinding<LatticeKind>::momentum_hamiltonian_y_derivative(
    Vec2d k) const {
  return this->momentum_hamiltonian_base(k, [](Vec2d k, Vec2d delta) {
    double phase = k.dot(delta);
    return (-1.0i * delta.y) * std::exp(-1.0i * phase);
  });
}