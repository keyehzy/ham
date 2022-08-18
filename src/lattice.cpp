#include <ham/lattice.h>
#include <ham/types.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

static int modulo(int a, int b) { return (a % b + b) % b; }

bool Lattice::is_inside_graph(int x, int y) const {
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

Matrix<int> Lattice::adjacency_matrix() const {
  Matrix<int> adj(this->orbital_count(), this->orbital_count());

  for (int site_index = 0; site_index < this->site_count(); site_index++) {
    for (int neighbor = 0; neighbor < nearest_neighbors_size(); neighbor++) {
      int neighbor_index = m_sites[site_index].neighbors[neighbor].index;
      for (int orbital_index = 0; orbital_index < m_orbitals; orbital_index++) {
        adj(site_index * m_orbitals + orbital_index,
            neighbor_index * m_orbitals + orbital_index) = 1;
      }
    }
  }

  return adj;
}

void Lattice::compute_graph() {
  for (int j = 0; j < m_ny; j++) {
    Vec2d current_position = static_cast<double>(j) * lattice_vector_2();

    for (int i = 0; i < m_nx; i++) {
      int index = j * m_nx + i;
      int kind = index % unitcell_size();

      m_sites[index].position = current_position;
      current_position = current_position + delta(position_walk(kind));

      for (int neighbor = 0; neighbor < nearest_neighbors_size(); neighbor++) {
        int dx = graph_walk(kind, neighbor, 0);
        int dy = graph_walk(kind, neighbor, 1);

        if (is_inside_graph(i + dx, j + dy)) {
          int neighbor_index =
              modulo(j + dy, m_ny) * m_nx + modulo(i + dx, m_nx);
          m_sites[index].neighbors.emplace_back(
              neighbor_index, deltas_for_walk(kind, neighbor));
        }
      }
    }
  }
}

Matrix<double> Lattice::position_operator_x() const {
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

Matrix<double> Lattice::position_operator_y() const {
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
