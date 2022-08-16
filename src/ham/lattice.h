#pragma once

#include <ham/matrix.h>
#include <ham/types.h>
#include <ham/vec.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <vector>

// Edge represent a "hop" between two sites in the graph. We store the index
// where the site lands and the direction used. direction is a index for looking
// up in a table of vectors after.
struct Edge {
  int index;
  int direction;
};

struct Site {
  std::vector<Edge> neighbors;
  Vec2<double> position;
};

class GrapheneLattice {
 public:
  // Kind of periodic boundary conditions:
  // Closed means periodic boundary condition in both X and Y directions.
  // Open no periodic boundary condition.
  // Open_x and Open_y are the options where just one of then are closed/open.
  enum class Boundary {
    Closed,
    Open,
    Open_x,
    Open_y,
  };

  static constexpr int unitcell_size = 4;
  static constexpr int nearest_neighbors_size = 3;

  explicit GrapheneLattice(int nx, int ny);

  int nx() const { return m_nx; }
  int ny() const { return m_ny; }
  int orbitals() const { return m_orbitals; }
  Boundary boundary() const { return m_boundary; }
  const std::vector<Site>& sites() const { return m_sites; }

  Site site(int i) const {
    assert(i >= 0 && i < static_cast<int>(m_sites.size()));
    return m_sites[i];
  }

  const std::array<Vec2<double>, 2 * nearest_neighbors_size>& deltas() const {
    return m_deltas;
  }

  Vec2<double> delta(int i) const {
    assert(i >= 0 && i < 2 * nearest_neighbors_size);
    return m_deltas[i];
  }

  int site_count() const { return m_nx * m_ny; }
  int orbital_count() const { return m_nx * m_ny * m_orbitals; }

  Matrix<int> adjacency_matrix() const;

 private:
  bool is_inside_graph(int x, int y) const;
  void compute_graph();

  int m_nx;
  int m_ny;
  int m_orbitals;
  Boundary m_boundary;
  std::vector<Site> m_sites;
  std::array<Vec2<double>, 2 * nearest_neighbors_size> m_deltas;
};

struct TightBindingParameters {
  double t = 1.0;
};

class GrapheneTightbinding {
 public:
  explicit GrapheneTightbinding(int nx, int ny,
                                TightBindingParameters parameters)
      : m_lattice(GrapheneLattice(nx, ny)), m_parameters(parameters){};

  Matrix<double> realspace_hamiltonian() const;

  using FactorFn = std::function<Complex(Vec2<double>, Vec2<double>)>;

  Matrix<Complex> momentum_hamiltonian_base(Vec2<double> k,
                                            FactorFn factor) const;

  Matrix<Complex> momentum_hamiltonian(Vec2<double> k) const;

  Matrix<Complex> momentum_hamiltonian_x_derivative(Vec2<double> k) const;

  Matrix<Complex> momentum_hamiltonian_y_derivative(Vec2<double> k) const;

  int size() const { return m_lattice.orbital_count(); }

 private:
  GrapheneLattice m_lattice;
  TightBindingParameters m_parameters;
};
