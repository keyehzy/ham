#pragma once

#include <ham/matrix.h>
#include <ham/types.h>
#include <ham/vec.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <numbers>

// Edge represent a "hop" between two sites in the graph. We store the index
// where the site lands and the direction used. direction is a index for looking
// up in a table of vectors after.
struct Edge {
  int index;
  int direction;
};

struct Site {
  std::vector<Edge> neighbors;
  Vec2d position;
};

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

class Lattice {
 public:
  explicit Lattice(int nx, int ny, int unitcell_size,
                   int nearest_neighbors_size, int orbitals = 1,
                   Boundary boundary = Boundary::Closed)
      : m_nx(nx),
        m_ny(ny),
        m_unitcell_size(unitcell_size),
        m_nearest_neighbors_size(nearest_neighbors_size),
        m_orbitals(orbitals),
        m_boundary(boundary),
        m_sites(nx * ny) {}

  int nx() const { return m_nx; }
  int ny() const { return m_ny; }
  int orbitals() const { return m_orbitals; }
  Boundary boundary() const { return m_boundary; }
  const std::vector<Site>& sites() const { return m_sites; }

  Site site(std::size_t i) const { return m_sites[i]; }

  int site_count() const { return m_nx * m_ny; }

  int orbital_count() const { return m_nx * m_ny * m_orbitals; }

  int unitcell_size() const { return m_unitcell_size; }

  int nearest_neighbors_size() const { return m_nearest_neighbors_size; }

  Matrix<int> adjacency_matrix() const;

  Matrix<double> position_operator_x() const;

  Matrix<double> position_operator_y() const;

  virtual Vec2d delta(std::size_t i) const = 0;

  virtual Vec2d lattice_vector_1() const = 0;

  virtual Vec2d lattice_vector_2() const = 0;

 protected:
  bool is_inside_graph(int x, int y) const;
  void compute_graph();

  virtual int graph_walk(int cell_offset, int nearest_neighbor_index,
                         int dir) const = 0;

  virtual int deltas_for_walk(int cell_offset,
                              int nearest_neighbor_index) const = 0;

  virtual int position_walk(int cell_offset) const = 0;

  int m_nx;
  int m_ny;
  int m_unitcell_size;
  int m_nearest_neighbors_size;
  int m_orbitals;
  Boundary m_boundary;
  std::vector<Site> m_sites;
};

class GrapheneLattice final : public Lattice {
 private:
  // clang-format off

  // This computes the graph of honeycomb (i.e. graphene-like) lattice. we
  // divide the 4-site unitcell in A,B,C,D. Each "site kind" defines three
  // different edges. We treat edges in two ways: graph_walk is the
  // coordinatized walk on the graph, it means that for a site of kind K its
  // neighboring edges are defined in a grid. deltas_for_walk is the actual
  // physical walk where this edge is defined. graph_walk is useful for defining
  // relations between sites on the lattice while deltas_for_walk for actual
  // physical calculations afterwards, e.g. fourier transforms.

  static constexpr int s_unitcell_size = 4;
  static constexpr int s_nearest_neighbors_size = 3;

  // Nearest-neighbors distance offset from the pespective of the graph itself
  // i.e. relative to the coordinates of the graph. The first site A has three
  // neighbors: The first one is in the same row one to the right, thus the
  // relative coordinate is (1, 0). The second one is in the row below one to
  // the right, thus the relative coordinate is (1, -1) and so on.
  static constexpr int s_graph_walk[s_unitcell_size][s_nearest_neighbors_size][2] =
      {{{1, 0}, {1, -1}, {-1, 0}},
       {{1, 0}, {-1, 1}, {-1, 0}},
       {{1, 1}, {1, 0}, {-1, 0}},
       {{1, 0}, {-1, 0}, {-1, -1}}};

  // Nearest-neighbors real space distance offset, indexing from the deltas.
  // Note that here we need to make a correction because sites A and C have
  // three hoppings which are the same and are opposite from the sites B and D.
  static constexpr int
      s_deltas_for_walk[s_unitcell_size][s_nearest_neighbors_size] = {
          {0, 1, 2}, {3, 4, 5}, {0, 1, 2}, {3, 4, 5}};

  // Auxilary array used to calculate where the next site will be placed on. We
  // always start from the left-most site and go to the right. In graphene we
  // make the armchair shape: up-right, right, down-right, right. This is
  // represented as indices to the deltas.
  static constexpr int s_position_walk[s_unitcell_size] = {0, 3, 1, 3};

  // Nearest-neighbors vectors for honeycomb lattice with 4-sites in the
  // unitcell. The first three are the usual delta_1, delta_2, delta_3 and
  // following it we have -delta_3, -delta_2, -delta_2 i.e the negative in the
  // opposite order.
  static constexpr double sqrt3 = std::numbers::sqrt3;
  static constexpr std::array<Vec2d, 2 * s_nearest_neighbors_size> s_deltas = {
        Vec2d{0.5, 0.5 * sqrt3},
        Vec2d{0.5, -0.5 * sqrt3},
        Vec2d{-1.0, 0.0},
        Vec2d{1.0, 0.0},
        Vec2d{-0.5, 0.5 * sqrt3},
        Vec2d{-0.5, -0.5 * sqrt3},
    };

  // Vectors corresponding to the translational invariance of the lattice. By
  // translating each site by (i * lattice_vector_1 + j * lattice_vector_2)
  // should land on the same corresponding site.
  static constexpr Vec2d s_lattice_vector_1 = Vec2d{3.0, 0.0};
  static constexpr Vec2d s_lattice_vector_2 = Vec2d{0.0, sqrt3};

  // clang-format on

 public:
  explicit GrapheneLattice(int nx, int ny)
      : Lattice(nx * s_unitcell_size, ny, s_unitcell_size,
                s_nearest_neighbors_size) {
    this->compute_graph();
  }

  Vec2d delta(std::size_t i) const override { return s_deltas[i]; }

  virtual Vec2d lattice_vector_1() const override { return s_lattice_vector_1; }

  virtual Vec2d lattice_vector_2() const override { return s_lattice_vector_2; }

 private:
  int graph_walk(int cell_offset, int nearest_neighbor_index,
                 int dir) const override {
    return s_graph_walk[cell_offset][nearest_neighbor_index][dir];
  }

  int deltas_for_walk(int cell_offset,
                      int nearest_neighbor_index) const override {
    return s_deltas_for_walk[cell_offset][nearest_neighbor_index];
  }

  int position_walk(int cell_offset) const override {
    return s_position_walk[cell_offset];
  }
};

class SquareLattice final : public Lattice {
 private:
  static constexpr int s_unitcell_size = 1;
  static constexpr int s_nearest_neighbors_size = 4;

  // clang-format off
  static constexpr int s_graph_walk[s_unitcell_size][s_nearest_neighbors_size][2] =
      {{{1, 0}, {0, 1}, {-1, 0}, {0, -1}}};

  static constexpr int
      s_deltas_for_walk[s_unitcell_size][s_nearest_neighbors_size] = {
          {0, 1, 2, 3}};
  // clang-format on

  static constexpr int s_position_walk[s_unitcell_size] = {0};
  static constexpr std::array<Vec2d, s_nearest_neighbors_size> s_deltas = {
      Vec2d{1.0, 0.0},
      Vec2d{0.0, 1.0},
      Vec2d{-1.0, 0.0},
      Vec2d{0.0, -1.0},
  };

  static constexpr Vec2d s_lattice_vector_1 = Vec2d{1.0, 0.0};
  static constexpr Vec2d s_lattice_vector_2 = Vec2d{0.0, 1.0};

 public:
  explicit SquareLattice(int nx, int ny)
      : Lattice(nx * s_unitcell_size, ny, s_unitcell_size,
                s_nearest_neighbors_size) {
    this->compute_graph();
  }

  Vec2d delta(std::size_t i) const override { return s_deltas[i]; }

  virtual Vec2d lattice_vector_1() const override { return s_lattice_vector_1; }

  virtual Vec2d lattice_vector_2() const override { return s_lattice_vector_2; }

 private:
  int graph_walk(int cell_offset, int nearest_neighbor_index,
                 int dir) const override {
    return s_graph_walk[cell_offset][nearest_neighbor_index][dir];
  }

  int deltas_for_walk(int cell_offset,
                      int nearest_neighbor_index) const override {
    return s_deltas_for_walk[cell_offset][nearest_neighbor_index];
  }

  int position_walk(int cell_offset) const override {
    return s_position_walk[cell_offset];
  }
};