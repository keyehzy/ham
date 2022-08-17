#include <ham/matrix.h>
#include <ham/types.h>

#include <complex>

template <>
Matrix<double> Matrix<double>::dot(const Matrix<double>& m) const {
  assert(m_cols == m.rows());

  std::size_t rows = m_rows;
  std::size_t cols = m.cols();
  Matrix<double> res(rows, cols);

  CBLAS_LAYOUT layout = CblasRowMajor;
  CBLAS_TRANSPOSE transA = CblasNoTrans;
  CBLAS_TRANSPOSE transB = CblasNoTrans;
  blasint M = this->rows();
  blasint N = m.cols();
  blasint K = this->cols();
  double alpha = static_cast<double>(1);
  const double* A = &(this->data()[0]);
  blasint lda = K;
  const double* B = &(m.data()[0]);
  blasint ldb = N;
  double beta = static_cast<double>(0);
  const double* C = &(res.data()[0]);
  blasint ldc = N;

  cblas_dgemm(layout, transA, transB, M, N, K, alpha, (double*)A, lda,
              (double*)B, ldb, beta, (double*)C, ldc);

  return res;
}

template <>
Matrix<Complex> Matrix<Complex>::dot(const Matrix<Complex>& m) const {
  assert(m_cols == m.rows());

  std::size_t rows = m_rows;
  std::size_t cols = m.cols();
  Matrix<Complex> res(rows, cols);

  CBLAS_LAYOUT layout = CblasRowMajor;
  CBLAS_TRANSPOSE transA = CblasNoTrans;
  CBLAS_TRANSPOSE transB = CblasNoTrans;
  blasint M = this->rows();
  blasint N = m.cols();
  blasint K = this->cols();
  Complex alpha = static_cast<Complex>(1);
  const Complex* A = &(this->data()[0]);
  blasint lda = K;
  const Complex* B = &(m.data()[0]);
  blasint ldb = N;
  Complex beta = static_cast<Complex>(0);
  const Complex* C = &(res.data()[0]);
  blasint ldc = N;

  cblas_zgemm(layout, transA, transB, M, N, K, &alpha, (Complex*)A, lda,
              (Complex*)B, ldb, &beta, (Complex*)C, ldc);

  return res;
}

template <>
std::vector<double> Matrix<double>::dot(const std::vector<double>& v) const {
  assert(m_rows == v.size());

  // std::cout << "uses blas" << std::endl;
  std::vector<double> res(m_cols);

  CBLAS_LAYOUT layout = CblasRowMajor;
  CBLAS_TRANSPOSE trans = CblasNoTrans;
  blasint M = this->rows();
  blasint N = this->cols();
  double alpha = static_cast<double>(1);
  const double* A = &(this->data()[0]);
  blasint lda = N;
  const double* X = &v[0];
  blasint incx = 1;
  double beta = static_cast<double>(0);
  const double* Y = &res[0];
  blasint incy = 1;

  cblas_dgemv(layout, trans, M, N, alpha, (double*)A, lda, (double*)X, incx,
              beta, (double*)Y, incy);

  return res;
}

template <>
std::vector<Complex> Matrix<Complex>::dot(const std::vector<Complex>& v) const {
  assert(m_rows == v.size());

  // std::cout << "uses blas" << std::endl;
  std::vector<Complex> res(m_cols);

  CBLAS_LAYOUT layout = CblasRowMajor;
  CBLAS_TRANSPOSE trans = CblasNoTrans;
  blasint M = this->rows();
  blasint N = this->cols();
  Complex alpha = static_cast<Complex>(1);
  const Complex* A = &(this->data()[0]);
  blasint lda = N;
  const Complex* X = &v[0];
  blasint incx = 1;
  Complex beta = static_cast<Complex>(0);
  const Complex* Y = &res[0];
  blasint incy = 1;

  cblas_zgemv(layout, trans, M, N, &alpha, (Complex*)A, lda, (Complex*)X, incx,
              &beta, (Complex*)Y, incy);

  return res;
}