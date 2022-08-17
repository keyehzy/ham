#include <ham/eigensystem.h>
#include <ham/lattice.h>

#include <iostream>
#include "ham/geometry.h"

void print_bandstructure(const GrapheneTightbinding& tb) {
  int px = 25, py = 25;
  for (int i = 0; i < px; i++) {
    double x =
      -M_PI +
      2.0 * M_PI * (static_cast<double>(i) / static_cast<double>(px-1));
    for (int j = 0; j < py; j++) {
      double y =
        -M_PI + 2.0 * M_PI *
        (static_cast<double>(j) / static_cast<double>(py-1));

      Vec2<double> k{x, y};
      Matrix<std::complex<double>> h = tb.momentum_hamiltonian(k);
      Eigensystem s(h);

      std::cout << k.x << " " << k.y << " ";
      for (std::size_t i = 0; i < s.size(); i++) {
        std::cout << s.eigenvalue(i) << " ";
      }
      std::cout << '\n';
    }
  }

}

int main(int argc, char** argv) {
  TightBindingParameters p;
  GrapheneTightbinding tb(1, 1, p);
  print_bandstructure(tb);



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
