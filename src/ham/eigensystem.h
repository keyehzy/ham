#pragma once

#include <ham/matrix.h>

#include <complex>
#include <lapacke.h>
#include <vector>
#include <cassert>

class Eigensystem {
 public:
  Eigensystem(const Matrix<std::complex<double>>& m) {
    assert(m.rows() == m.cols());

    int matrix_layout = LAPACK_ROW_MAJOR;
    char jobz = 'N'; // @@@ change to 'V' for eigenvectors
    char uplo = 'U';
    lapack_int n = m.rows();

    m_eigenvalues.resize(n);
    int lda = n;
    LAPACKE_zheev(matrix_layout, jobz, uplo, n,
                  (_Complex double*)&(m.data()[0]), lda, &m_eigenvalues[0]);
  }

  const std::vector<double>& eigenvalues() const { return m_eigenvalues; }
  double eigenvalue(int i) const { return m_eigenvalues[i]; }
  std::size_t size() const { return m_eigenvalues.size(); }

 private:
  std::vector<double> m_eigenvalues;
  // eigenvectors @@@
};
