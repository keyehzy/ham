#include <ham/geometry.h>

#include <algorithm>
#include <cassert>
#include <numeric>

std::vector<double> linspace(double start, double stop, int num) {
  int div = num - 1;
  assert(div > 0);

  std::vector<double> l(num);
  std::iota(l.begin(), l.end(), 0.0);

  double step = (stop - start) / div;

  std::transform(l.begin(), l.end(), l.begin(),
                 [&](double val) { return val * step + start; });

  return l;
}
