#include <ham/vector.h>

#include <cassert>

// FIXME: we can use BLAS for this
Complex vdot(const Vector<Complex>& v1, const Vector<Complex>& v2) {
  assert(v1.size() == v2.size());

  Complex cx = 0;

  for (std::size_t i = 0; i < v1.size(); i++) {
    cx += std::conj(v1[i]) * v2[i];
  }

  return cx;
}
