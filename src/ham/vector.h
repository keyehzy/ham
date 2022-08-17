#pragma once

#include <ham/types.h>

#include <numeric>
#include <vector>

template <typename T>
class Vector {
 public:
  Vector(std::size_t s) : m_size(s), m_data(s){};
  Vector(const Vector&);
  ~Vector(){};

  std::size_t size() const { return m_size; }
  const std::vector<T>& data() const { return m_data; }

  T& operator[](std::size_t i) { return m_data[i]; };
  T operator[](std::size_t i) const { return m_data[i]; }

  T dot(const Vector<T>& v) const;
  T total() const;

 private:
  std::size_t m_size;
  std::vector<T> m_data;
};

template <typename T>
T Vector<T>::dot(const Vector<T>& v) const {
  assert(m_size == v.size());
  return std::inner_product(m_data.begin(), m_data.end(), v.data().begin(),
                            static_cast<T>(0));
}

template <typename T>
T Vector<T>::total() const {
  return std::accumulate(m_data.begin(), m_data.end(), static_cast<T>(0));
}

Complex vdot(const Vector<Complex>&, const Vector<Complex>&);
