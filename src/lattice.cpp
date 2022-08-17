#include <ham/lattice.h>
#include <ham/types.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

static int modulo(int a, int b) { return (a % b + b) % b; }

static constexpr Complex COMPI = Complex{0.0, 1.0};

bool GrapheneLattice::is_inside_graph(int x, int y) const {
  switch (m_boundary) {
    case Boundary::Open:
      return (x >= 0 && x < m_nx) && (y >= 0 && y < m_ny);
    case Boundary::Closed:
      return 1;
    case Boundary::Open_x:
      return (x >= 0 && x < m_nx);
    case Boundary::Open_y:
      return (y >= 0 && x < m_ny);
    default:
      return 0;
  }
}

GrapheneLattice::GrapheneLattice(int nx, int ny)
    : m_nx(nx * unitcell_size),
      m_ny(ny),
      m_orbitals(1),
      m_boundary(Boundary::Closed) {
  // Nearest-neighbors vectors for honeycomb lattice with 4-sites in the
  // unitcell
  m_deltas = {
      Vec2d{0.5, 0.5 * std::sqrt(3.0)},
      Vec2d{0.5, -0.5 * std::sqrt(3.0)},
      Vec2d{-1.0, 0.0},
      Vec2d{1.0, 0.0},
      Vec2d{-0.5, 0.5 * std::sqrt(3.0)},
      Vec2d{-0.5, -0.5 * std::sqrt(3.0)},
  };

  // Total number of sites
  m_sites.resize(m_nx * m_ny);

  // Finally compute the graph
  this->compute_graph();
}

Matrix<int> GrapheneLattice::adjacency_matrix() const {
  Matrix<int> adj(this->orbital_count(), this->orbital_count());

  for (int site_index = 0; site_index < this->site_count(); site_index++) {
    for (int neighbor = 0; neighbor < nearest_neighbors_size; neighbor++) {
      int neighbor_index = m_sites[site_index].neighbors[neighbor].index;
      for (int orbital_index = 0; orbital_index < m_orbitals; orbital_index++) {
        adj(site_index * m_orbitals + orbital_index,
            neighbor_index * m_orbitals + orbital_index) = 1;
      }
    }
  }

  return adj;
}

void GrapheneLattice::compute_graph() {
  // this computes the graph of honeycomb (i.e. graphene-like) lattice. we
  // divide the 4-site unitcell in A,B,C,D. each "site kind" defines three
  // different edges. we treat edges in two ways: graph_walk is the
  // coordinatized walk on the graph, it means that for a site of kind K its
  // neighboring edges are defined in a grid. deltas_for_walk is the actual
  // physical walk where this edge is defined. graph_walk is useful for
  // defining relations between sites on the lattice while deltas_for_walk for
  // actual physical calculations afterwards, e.g. fourier transforms.

  // first we calculate the nearest-neighbors. i think we can calculate this
  // values, don't need to hardcode it
  constexpr int graph_walk[unitcell_size][nearest_neighbors_size][2] = {
      {{1, 0}, {1, -1}, {-1, 0}},
      {{1, 0}, {-1, 1}, {-1, 0}},
      {{1, 1}, {1, 0}, {-1, 0}},
      {{1, 0}, {-1, 0}, {-1, -1}}};

  // here we need to make a correction because sites A and C have three
  // hoppings that are opposite from the sites B and D.
  constexpr int deltas_for_walk[unitcell_size][nearest_neighbors_size] = {
      {0, 1, 2}, {3, 4, 5}, {0, 1, 2}, {3, 4, 5}};

  // position_walk is an auxilary array used to calculate where the next
  // site will be placed on.
  const Vec2d position_walk[unitcell_size] = {m_deltas[0], m_deltas[3],
                                              m_deltas[1], m_deltas[3]};

  for (int j = 0; j < m_ny; j++) {
    // Start from the left-most site
    Vec2d current_position{0.0, static_cast<double>(j) * std::sqrt(3.0)};

    for (int i = 0; i < m_nx; i++) {
      int index = j * m_nx + i;
      int kind = index % unitcell_size;

      m_sites[index].position = current_position;
      current_position = current_position + position_walk[kind];

      for (int neighbor = 0; neighbor < nearest_neighbors_size; neighbor++) {
        int dx = graph_walk[kind][neighbor][0];
        int dy = graph_walk[kind][neighbor][1];

        if (is_inside_graph(i + dx, j + dy)) {
          int neighbor_index =
              modulo(j + dy, m_ny) * m_nx + modulo(i + dx, m_nx);
          m_sites[index].neighbors.emplace_back(
              neighbor_index, deltas_for_walk[kind][neighbor]);
        }
      }
    }
  }
}

Matrix<double> GrapheneLattice::position_operator_x() const {
  Matrix<double> X(this->orbital_count(), this->orbital_count());

  for (int site_index = 0; site_index < this->site_count(); site_index++) {
    for (int orbital_index = 0; orbital_index < this->orbitals();
         orbital_index++) {
      X(site_index * this->orbitals() + orbital_index,
        site_index * this->orbitals() + orbital_index) =
          site(site_index).position.x;
    }
  }

  return X;
}

Matrix<double> GrapheneLattice::position_operator_y() const {
  Matrix<double> Y(this->orbital_count(), this->orbital_count());

  for (int site_index = 0; site_index < this->site_count(); site_index++) {
    for (int orbital_index = 0; orbital_index < this->orbitals();
         orbital_index++) {
      Y(site_index * this->orbitals() + orbital_index,
        site_index * this->orbitals() + orbital_index) =
          site(site_index).position.y;
    }
  }

  return Y;
}

Matrix<double> GrapheneTightbinding::realspace_hamiltonian() const {
  Matrix<double> h(this->size(), this->size());

  // 1) For each site, go to its neighbors and for each neighbor set a value at
  // the matrix entry (site_index, neighbor_index).
  // 2) For multiple orbitals we need to pad the entries.

  for (int site_index = 0; site_index < m_lattice.site_count(); site_index++) {
    for (int neighbor = 0; neighbor < m_lattice.nearest_neighbors_size;
         neighbor++) {
      int neighbor_index =
          m_lattice.site(site_index).neighbors[neighbor].index;
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

Matrix<Complex> GrapheneTightbinding::momentum_hamiltonian_base(
    Vec2d k, FactorFn factor) const {
  // We can only perform this operation when we have closed periodic boundary
  // conditions, otherwise we cannot use Bloch theorem a justify using
  // tightbinding.
  assert(m_lattice.boundary() == GrapheneLattice::Boundary::Closed);

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
    for (int neighbor = 0; neighbor < m_lattice.nearest_neighbors_size;
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

Matrix<Complex> GrapheneTightbinding::momentum_hamiltonian(Vec2d k) const {
  return this->momentum_hamiltonian_base(k, [](Vec2d k, Vec2d delta) {
    double phase = k.dot(delta);
    return std::exp(-COMPI * phase);
  });
}

Matrix<Complex> GrapheneTightbinding::momentum_hamiltonian_x_derivative(
    Vec2d k) const {
  return this->momentum_hamiltonian_base(k, [](Vec2d k, Vec2d delta) {
    double phase = k.dot(delta);
    return (-COMPI * delta.x) * std::exp(-COMPI * phase);
  });
}

Matrix<Complex> GrapheneTightbinding::momentum_hamiltonian_y_derivative(
    Vec2d k) const {
  return this->momentum_hamiltonian_base(k, [](Vec2d k, Vec2d delta) {
    double phase = k.dot(delta);
    return (-COMPI * delta.y) * std::exp(-COMPI * phase);
  });
}
