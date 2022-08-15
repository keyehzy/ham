#include <vector>

template <class T>
class Matrix {
 public:
  Matrix(std::size_t r, std::size_t c) : m_rows(r), m_cols(c), m_data(r * c){};
  Matrix(const Matrix&);
  ~Matrix(){};

  std::size_t size() const { return m_rows * m_cols; }
  std::size_t rows() const { return m_rows; }
  std::size_t cols() const { return m_cols; }
  std::vector<T>& data() { return m_data; }
  const std::vector<T>& data() const { return m_data; }

  T& operator()(std::size_t i, std::size_t j);
  T operator()(std::size_t i, std::size_t j) const;

 private:
  std::size_t m_rows;
  std::size_t m_cols;
  std::vector<T> m_data;
};

template <class T>
inline T& Matrix<T>::operator()(std::size_t i, std::size_t j) {
  return m_data[m_cols * i + j];
}

template <class T>
inline T Matrix<T>::operator()(std::size_t i, std::size_t j) const {
  return m_data[m_cols * i + j];
}
