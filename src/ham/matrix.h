#pragma once

#include <ham/vector.h>

#include <cassert>
#include <cstddef>
#include <valarray>
#include <vector>

#include "cblas.h"

template <class T>
class Slice {
 public:
  Slice(const std::vector<T>& data, const std::slice& slice)
      : m_data(data), m_slice(slice), m_current(0){};

  Slice begin() const {
    Slice it = *this;
    it.m_current = 0;
    return it;
  }

  Slice end() const {
    Slice it = *this;
    it.current = m_slice.size();
    return it;
  }

  Slice& operator++() {
    m_current++;
    return *this;
  }

  Slice operator++(int) {
    Slice it = *this;
    m_current++;
    return it;
  }

  Vector<T> as_vec() const {
    Vector<T> v(m_slice.size());
    for (std::size_t i = 0; i < m_slice.size(); i++) {
      v[i] = this->ref(i);
    }
    return v;
  }

  const T& operator[](std::size_t i) const { return this->ref(i); }

  const T& operator*() const { return this->ref(m_current); }

  const T& ref(std::size_t index) const {
    return m_data[m_slice.start() + index * m_slice.stride()];
  }

 private:
  const std::vector<T>& m_data;
  const std::slice m_slice;
  std::size_t m_current;
};

template <class T>
class Matrix {
 public:
  Matrix(){};
  Matrix(std::size_t r, std::size_t c) : m_rows(r), m_cols(c), m_data(r * c){};
  Matrix(const Matrix&);
  ~Matrix(){};

  std::size_t size() const { return m_rows * m_cols; }
  std::size_t rows() const { return m_rows; }
  std::size_t cols() const { return m_cols; }

  const std::vector<T>& data() const { return m_data; }

  Slice<T> row(std::size_t i) const;
  Slice<T> col(std::size_t i) const;
  Slice<T> diag() const;

  T& operator()(std::size_t i, std::size_t j);
  T operator()(std::size_t i, std::size_t j) const;

  Matrix<T> dot(const Matrix<T>& m) const;
  Vector<T> dot(const Vector<T>& v) const;

  void from(const Matrix& m) {
    m_cols = m.rows();
    m_cols = m.cols();
    m_data = m.data();
  }

 private:
  std::size_t m_rows;
  std::size_t m_cols;
  std::vector<T> m_data;
};

template <class T>
inline Slice<T> Matrix<T>::row(std::size_t i) const {
  return Slice<T>(m_data, std::slice(m_cols * i, m_cols, 1));
}

template <class T>
inline Slice<T> Matrix<T>::col(std::size_t i) const {
  return Slice<T>(m_data, std::slice(i, m_rows, m_cols));
}

template <class T>
inline Slice<T> Matrix<T>::diag() const {
  return Slice<T>(m_data, std::slice(0, m_rows, m_cols + 1));
}

template <class T>
inline T& Matrix<T>::operator()(std::size_t i, std::size_t j) {
  return m_data[m_cols * i + j];
}

template <class T>
inline T Matrix<T>::operator()(std::size_t i, std::size_t j) const {
  return m_data[m_cols * i + j];
}
