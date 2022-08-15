#include <ham/lattice.h>
#include <lapacke.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

class Eigensystem {
 public:
  Eigensystem(const Matrix<std::complex<double>>& m) {
    this->compute_eigenvalues(m);
  }

  const std::vector<double>& eigenvalues() const { return m_eigenvalues; }
  double eigenvalue(int i) const { return m_eigenvalues[i]; }
  std::size_t size() const { return m_eigenvalues.size(); }

 private:
  void compute_eigenvalues(const Matrix<std::complex<double>>& m) {
    int matrix_layout = LAPACK_ROW_MAJOR;
    char jobz = 'V';
    char uplo = 'U';
    lapack_int n = m.rows();

    const std::complex<double>* a = &(m.data()[0]);
    m_eigenvalues.resize(n);
    int lda = n;
    LAPACKE_zheev(matrix_layout, jobz, uplo, n, (_Complex double*)a, lda, &m_eigenvalues[0]);
    //LAPACKE_zheev(matrix_layout, jobz, uplo, n, a, lda, w.data());
  }

  std::vector<double> m_eigenvalues;
  // eigenvectors @@@
};

int main(int argc, char** argv) {

  Matrix<std::complex<double>> m(3, 3);
  m(0,0) = 1.0;
  m(0,1) = 2.0;
  m(0,2) = 3.0;
  m(1,0) = 2.0;
  m(1,1) = 4.0;
  m(1,2) = 5.0;
  m(2,0) = 3.0;
  m(2,1) = 5.0;
  m(2,2) = 6.0;

  Eigensystem s(m);
  
  for (std::size_t i = 0; i < s.size(); i++) {
    std::cout << s.eigenvalue(i) << std::endl;
  }

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
