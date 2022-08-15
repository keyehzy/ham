#include <ham/matrix.h>
#include <ham/vec.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

struct Edge {
  int index;
  int direction;
};

struct Site {
  std::vector<Edge> neighbors;
  Vec2<double> position;
};

enum Boundary {
  Closed,
  Open,
  Open_x,
  Open_y,
};

static int modulo(int a, int b) { return (a % b + b) % b; }

class GrapheneLattice {
  static constexpr int unitcell_size = 4;
  static constexpr int nearest_neighbors_size = 3;

 public:
  explicit GrapheneLattice(int nx, int ny);

  Matrix<int> adjacency_matrix();
  int site_count() const { return m_nx * m_ny; }
  int orbital_count() const { return m_nx * m_ny * m_orbitals; }

 private:
  bool is_inside_graph(int x, int y);
  void compute_graph();

  int m_nx;
  int m_ny;
  int m_orbitals;
  Boundary m_boundary;
  std::vector<Site> m_sites;
  std::array<Vec2<double>, 2 * nearest_neighbors_size> m_deltas;
};

bool GrapheneLattice::is_inside_graph(int x, int y) {
  switch (m_boundary) {
    case Open:
      return (x >= 0 && x < m_nx) && (y >= 0 && y < m_ny);
    case Closed:
      return 1;
    case Open_x:
      return (x >= 0 && x < m_nx);
    case Open_y:
      return (y >= 0 && x < m_ny);
    default:
      return 0;
  }
}

GrapheneLattice::GrapheneLattice(int nx, int ny)
    : m_nx(nx * unitcell_size), m_ny(ny), m_orbitals(1), m_boundary(Closed) {
  // Nearest-neighbors vectors for honeycomb lattice with 4-sites in the
  // unitcell
  m_deltas = {
      Vec2<double>{0.5, 0.5 * std::sqrt(3.0)},
      Vec2<double>{0.5, -0.5 * std::sqrt(3.0)},
      Vec2<double>{-1.0, 0.0},
      Vec2<double>{1.0, 0.0},
      Vec2<double>{-0.5, 0.5 * std::sqrt(3.0)},
      Vec2<double>{-0.5, -0.5 * std::sqrt(3.0)},
  };

  // Total number of sites
  m_sites.resize(m_nx * m_ny);

  // Compute graph
  this->compute_graph();
}

Matrix<int> GrapheneLattice::adjacency_matrix() {
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
  const Vec2<double> position_walk[unitcell_size] = {m_deltas[0], m_deltas[3],
                                                     m_deltas[1], m_deltas[3]};

  for (int j = 0; j < m_ny; j++) {
    // Start from the left-most site
    Vec2<double> current_position{0.0, static_cast<double>(j) * std::sqrt(3.0)};

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

int main(int argc, char** argv) {
  GrapheneLattice l(2, 2);
  Matrix<int> adj = l.adjacency_matrix();

  for (int i = 0; i < l.orbital_count(); i++) {
    for (int j = 0; j < l.orbital_count(); j++) {
      std::cout << adj(i, j) << " ";
    }
    std::cout << '\n';
  }

  return 0;
}
