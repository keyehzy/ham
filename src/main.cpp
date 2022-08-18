#include <ham/eigensystem.h>
#include <ham/geometry.h>
#include <ham/lattice.h>
#include <ham/tightbinding.h>
#include <ham/types.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <numbers>

template <class LatticeKind>
void print_bandstructure(const Tightbinding<LatticeKind>& tb) {
  double pi = std::numbers::pi;
  double xfactor = (static_cast<double>(tb.lattice().nx()) /
                    static_cast<double>(tb.lattice().unitcell_size())) *
                   tb.lattice().lattice_vector_1().norm();
  double yfactor = static_cast<double>(tb.lattice().ny()) *
                   tb.lattice().lattice_vector_2().norm();
  RectangleGrid grid(-pi / xfactor, pi / xfactor, -pi / yfactor, pi / yfactor);

  for (std::size_t i = 0; i < grid.size(); i++) {
    for (std::size_t j = 0; j < grid.size(); j++) {
      Vec2d k = grid(i, j);
      Matrix<Complex> h = tb.momentum_hamiltonian(k);
      Eigensystem s(h);

      std::cout << k.x << " " << k.y << " ";
      std::copy(s.eigenvalues().begin(), s.eigenvalues().end(),
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << '\n';
    }
  }
}

int main(int argc, char** argv) {
  TightBindingParameters p;
  Tightbinding<GrapheneLattice> tb(2, 2, p);
  print_bandstructure(tb);

  // GrapheneLattice gl(1, 1);
  // Matrix<double> X = gl.position_operator_x();
  // Matrix<double> Y = gl.position_operator_y();

  // for (int i = 0; i < gl.orbital_count(); i++) {
  //   for (int j = 0; j < gl.orbital_count(); j++) {
  //     std::cout << X(i, j) << " ";
  //   }
  //   std::cout << '\n';
  // }

  // TightBindingParameters p;
  // GrapheneTightbinding tb(1, 1, p);

  // Matrix<std::complex<double>> h =
  //     tb.closed_momentum_hamiltonian(Vec2<double>{1.5, 2.0});

  // for (std::size_t i = 0; i < h.rows(); i++) {
  //   for (std::size_t j = 0; j < h.cols(); j++) {
  //     std::cout << h(i, j) << " ";
  //   }
  //   std::cout << '\n';
  // }

  return 0;
}
