#pragma once

#include <ham/matrix.h>
#include <ham/types.h>
#include <lapacke.h>

#include <cassert>
#include <complex>
#include <vector>

class Eigensystem {
 public:
  Eigensystem(const Matrix<Complex>& m);

  double eigenvalue(std::size_t i) const { return m_eigenvalues[i]; }

  Vector<Complex> eigenvector(std::size_t i) const {
    return m_eigenvectors.col(i).as_vec();
  }

  std::size_t size() const { return m_eigenvalues.size(); }

 private:
  std::vector<double> m_eigenvalues;
  Matrix<Complex> m_eigenvectors;
};
