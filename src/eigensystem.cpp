#include <ham/eigensystem.h>
#include <ham/types.h>

#include "lapacke.h"

Eigensystem::Eigensystem(const Matrix<Complex>& m) {
  assert(m.rows() == m.cols());

  int matrix_layout = LAPACK_ROW_MAJOR;
  char jobz = 'V';
  char uplo = 'U';
  lapack_int n = m.rows();
  m_eigenvalues.resize(n);
  m_eigenvectors.from(m);
  int lda = n;

  LAPACKE_zheev(matrix_layout, jobz, uplo, n,
                (_Complex double*)&(m_eigenvectors.data()[0]), lda,
                &m_eigenvalues[0]);
}
