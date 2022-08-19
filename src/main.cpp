#include <ham/eigensystem.h>
#include <ham/lattice.h>
#include <ham/tightbinding.h>
#include <ham/types.h>
#include <ham/dos.h>
#include <ham/brillouin_zone.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <numbers>

template <class LatticeKind>
void print_bandstructure(const Tightbinding<LatticeKind>& tb) {
  BrillouinZone<LatticeKind> BZ(tb.lattice());

  for (std::size_t i = 0; i < BZ.size(); i++) {
    for (std::size_t j = 0; j < BZ.size(); j++) {
      Vec2d k = BZ(i, j);
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
  Tightbinding<GrapheneLattice> tb(1, 1, p);
  // print_bandstructure(tb);

  DOS<GrapheneLattice> dos(tb, -3.5, 3.5);

  for (std::size_t i = 0; i < dos.size(); i++) {
    std::cout << dos.omega(i) << " " << dos(i) << std::endl;
  }
  
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
