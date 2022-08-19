#pragma once

#include <ham/eigensystem.h>
#include <ham/tightbinding.h>
#include <ham/types.h>
#include <ham/brillouin_zone.h>

using namespace std::complex_literals;

template <class LatticeKind>
class DOS {
  // small value used to compute the Green's function
  static constexpr Complex s_eta = 1.0e-2i;

 public:
  explicit DOS(const Tightbinding<LatticeKind>& tightbinding, double start,
               double stop, std::size_t num = 1000)
      : m_tightbinding(tightbinding),
        m_omega(linspace(start, stop, num)),
        m_dos(num),
        m_size(num) {

    BrillouinZone<LatticeKind> BZ(m_tightbinding.lattice());

    double norma = static_cast<double>(BZ.size() * BZ.size());

    for (std::size_t i = 0; i < BZ.size(); i++) {
      for (std::size_t j = 0; j < BZ.size(); j++) {
        Vec2d k = BZ(i, j);
        Matrix<Complex> h = tightbinding.momentum_hamiltonian(k);
        Eigensystem s(h);

        std::size_t states = m_tightbinding.lattice().orbital_count();
        for (std::size_t s1 = 0; s1 < states; s1++) {
          double val = s.eigenvalue(s1);
          Vector<Complex> vec = s.eigenvector(s1);
          for (std::size_t k = 0; k < m_omega.size(); k++) {
            for (std::size_t s2 = 0; s2 < states; s2++) {
              double v2 = std::pow(std::abs(vec[s2]), 2.0);
              m_dos[k] +=
                  -v2 * greens_function(val, m_omega[k]).imag() / pi / norma;
            }
          }
        }
      }
    }
  }

  double operator()(std::size_t i) const { return m_dos[i]; }

  double omega(std::size_t i) const { return m_omega[i]; }

  std::size_t size() const { return m_size; }

 private:
  Complex greens_function(double E, double omega) const {
    return 1.0 / (omega - E + s_eta);
  }

  Tightbinding<LatticeKind> m_tightbinding;
  std::vector<double> m_omega;
  std::vector<double> m_dos;
  std::size_t m_size;
};
