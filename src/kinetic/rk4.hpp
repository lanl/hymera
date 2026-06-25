#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "util/common.hpp"

template <class System>
KOKKOS_INLINE_FUNCTION void rk4_step(
    const System &f,
    const double t, typename System::value_type &y, double h, Kokkos::Array<typename System::value_type, 5> &k) {
    const size_t n = y.size();
    // 1) k[0] = f(t, y)
    f(t, y, k[1]);
    // 2) k[1] = f(t + h/2, y + h/2 * k[0])
    for (size_t i = 0; i < n; ++i)
        k[0][i] = y[i] + 0.5 * h * k[1][i];
    f(t + 0.5 * h, k[0], k[2]);
    // 3) k[2] = f(t + h/2, y + h/2 * k[1])
    for (size_t i = 0; i < n; ++i)
        k[0][i] = y[i] + 0.5 * h * k[2][i];
    f(t + 0.5 * h, k[0], k[3]);
    // 4) k[3] = f(t + h, y + h * k[2])
    for (size_t i = 0; i < n; ++i)
        k[0][i] = y[i] + h * k[3][i];
    f(t + h, k[0], k[4]);
    // Compute 4th-order update
    for (size_t i = 0; i < n; ++i) {
        y[i] = y[i] + (h / 6.0) *(k[1][i] + 2.0 * k[2][i] + 2.0 * k[3][i] + k[4][i]);
    }
}

template <class System, class Verificator, bool VERBOSE = false>
KOKKOS_INLINE_FUNCTION typename Verificator::ResultCode_t solve_rk4_fixed(
    const System &f, const Verificator &v, typename System::value_type &y, const double t0,
    const double tf, double h,
    Kokkos::Array<typename System::value_type, 5> &work) {
    const size_t n = y.size();
    double t = t0;
    int nstep = static_cast<int>(ceil((tf - t0) / h));
    for (int stepCount = 0; stepCount < nstep; ++stepCount) {
        auto st = v.verify(y);
        if (st != Verificator::Success) return st;

        if (stepCount == nstep - 1)
            h = tf - t;
        rk4_step(f, t, y, h, work);

        st = v.verify(y);
        if (st != Verificator::Success) return st;

        t += h;

        if constexpr (VERBOSE) {
            printf("%20.14le ", t);
            for (auto x : y)
                printf("%20.14le ", x);
            printf("A\n");
        }
    }
    return Verificator::Success;
}
