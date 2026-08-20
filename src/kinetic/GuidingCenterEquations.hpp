//========================================================================================
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
// for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================

#pragma once
#include "util/common.hpp"
#include "hFlux/common.hpp"
#include "kinetic/FieldEvaluator.hpp"

/// @brief Guiding center equations of motion. Functor for the right hand side.
/// Can be used with the ODE solver.
/// Can be constructed with the external magnetic and electric fields and status
/// of the particle reference. By constructing a functor: GuidingCenterEquations
/// guiding_center(field, status, \f$\frac{c}{a\omega_{ce,0}}\f$, ,
/// \f$\frac{c\tau_c}{a}\f$,  \f${\alpha_0}\f$);
template <class Field, bool EF = true, bool SlabModel = true> struct GuidingCenterEquations {
  /// External magnetic and electric fields reference.
  const Field field;
  /// \brief ODE coefficients.
  const Real c_aw0;
  const Real ct_a;
  const Real alpha0;

  GuidingCenterEquations(const Field field, const Real c_aw0, const Real ct_a, const Real alpha0):
    field(field), c_aw0(c_aw0), ct_a(ct_a), alpha0(alpha0) {}

  typedef Dim5 value_type;

  /** @brief Guiding center equations of motion, the right hand size of the ODE.
   * @param[in] X The state vector of the guiding center particle.
   * @param[out] dXdt The right hand side of the ODE.
   * @param[in] t The time.
   * @note The state vector X is in the form of (p, xi, R, phi, Z).
   * @note The guiding center equations of motion are given by
   *   \f{eqnarray*}{ *     \frac {dp}{dt} &=& \frac c{a\omega_{ce,0}} \frac{p (1 +
   * \xi^2)}{2B_\parallel^\ast}  (\vec E\times \vec b) \cdot \nabla \ln B\\
   *                    &-& \xi \frac{B^\ast \cdot E}{B_\parallel^\ast} \\
   *                    &-& \alpha_0 B^2 p \gamma (1-\xi^2) \left(1 + \frac
   * c{a\omega_{ce,0}}  \frac{p \xi}B  \vec b \cdot \left(\nabla \times \vec
   * b\right)\right)
   *     \frac {d\xi}{dt} &=& -\frac {c\tau_c}a (1 - \xi^2) \frac p\gamma \frac
   * 1{2B^\ast_\parallel} B^\ast \cdot \nabla \ln B - \frac 12
   * \frac{c}{a\omega_0 B_\parallel} \xi (1 - \xi^2) \nabla \ln B \cdot (\vec
   * E\times \vec b)
   *                       - \frac{(1 - xi^2)}{p B_\parallel} B^\ast \cdot E
   *                    &+& \alpha_0 B^2 \xi \frac{(1-\xi^2)}{\gamma} \left(1 +
   * \frac c{a\omega_0 B} \frac 12 p \xi  \vec b \cdot \nabla \times \vec
   * b\right)
   *                    &-& \alpha_0 B^2 \frac 12 \gamma (1-\xi^2)^2 \left(\frac
   * c{a\omega_0 B} \frac pB \vec b \cdot \nabla \times \vec b\right)
   *    \f}
   * @note The guiding center equations of motion are derived from the drift
   * kinetic equation by neglecting the parallel motion and the energy
   * diffusion.
   */

  KOKKOS_INLINE_FUNCTION
  void operator()(const Real &t, const value_type &X,
                        value_type &dXdt) const {

    const Real p =   X[0];
    const Real xi =  X[1];
    const Real R =   X[2];
    const Real phi = X[3];
    const Real Z =   X[4];
    const Real gamma = gamma_(p);

    if constexpr (SlabModel) {
      dXdt[0] = -xi * field.E_0 - alpha0 * p * gamma * (1.0 - xi * xi);
      dXdt[1] = -(1.0 - xi * xi) / p * field.E_0 + alpha0 * xi * (1.0 - xi * xi) / gamma;
      dXdt[2] = 0.0;
      dXdt[3] = 0.0;
      dXdt[4] = 0.0;
      return;
    }

    Kinetic::EvalGCE ev;
    field.eval(ev, R, Z, t);

    const Real Bsq = ev.Bsq;
    const Real Bmag = ev.Bmag;

    Dim3 b = {};
    for (int i = 0; i < 3; ++i) {
      b[i] = ev.B[i] / Bmag;
    }

    Dim3 gradlnB = {dot_product(ev.B, ev.dBdR) / Bsq, 0.0,
                    dot_product(ev.B, ev.dBdZ) / Bsq};

    Dim3 b_x_gradlnB;
    cross_product(b, gradlnB, b_x_gradlnB);

    Dim3 Bstar;
    for (int i = 0; i < 3; ++i) {
      Bstar[i] = ev.B[i] - c_aw0 * p * xi * (b_x_gradlnB[i] + ev.J[i] / Bmag);
      // Modify electric field by p_parallel * dbdt
      ev.E[i] -= p * xi * ev.dbdt[i];
    }
    if constexpr (EF == false) ev.E = {};

    Real Bpar = dot_product(b, Bstar);
    Real B_d_gradlnB = dot_product(ev.B, gradlnB);
    Real Bstar_d_gradlnB = dot_product(Bstar, gradlnB);

    /** \brief Using Curl product rule
     *   \f{eqnarray*}{
     *   \nabla \times \vec b &=&  \nabla \times (\frac 1B \vec B) \\
     *                        &=&  \frac 1B \nabla \times \vec B  + \nabla \frac
     * 1B \times \vec B \\
     *                        &=&  \frac 1B \nabla \times \vec B  - \frac 1B
     * \nabla \ln B \times \vec B \\
     *                        &=&  \frac 1B \nabla \times \vec B  - \nabla \ln B
     * \times \vec b \\
     *                        &=&  \frac 1B \nabla \times \vec B  + \vec b
     * \times
     * \nabla \ln B.
     *   \f}
     *   Compuing \f[\vec b \cdot (\nabla \times b)\f], noting that \f[\vec b
     * \cdot (\vec b \times \nabla \ln B) \equiv 0 \f] is
     *
     *   Multiplying by \f[\frac{p[kg\,m/s]}{eB[T]} =
     * \frac{c}{a\omega_{ce,0}}\frac{p}{B}\f]
     *
     *   For consistency with old terms, set this term to 0
     *
     */
    Real b_rad_term = 0.0; // dot_product(b, curlB) * c_aw0 * p / Bsq;

    Dim3 E_x_b;
    cross_product(ev.E, b, E_x_b);

    Real gradlnB_d_Exb = dot_product(gradlnB, E_x_b);
    Real Bstar_d_E = dot_product(Bstar, ev.E);

    Real one_m_xisq = 1.0 - xi * xi;
    /**  \brief Momentum evolution
     *   \f{eqnarray*}{
     *     \frac {dp}{dt} &=& \frac c{a\omega_{ce,0}} \frac{p (1 +
     * \xi^2)}{2B_\parallel^\ast}  (\vec E\times \vec b) \cdot \nabla \ln B
     *                       - \xi \frac{B^\ast \cdot E}{B_\parallel^\ast} \\
     *                    &-& \alpha_0 B^2 p \gamma (1-\xi^2) \left(1 + \xi
     * \text{rad_term}\right)
     *   \f}
     *
     *   where \f[\text{rad_term} = \frac{c}{a\omega_{ce,0}}\frac{p}{B} \vec b
     * \cdot (\nabla \times b) \f]
     *
     */

    dXdt[0] = c_aw0 * p * one_m_xisq * 0.5 / Bpar * gradlnB_d_Exb -
              xi * Bstar_d_E / Bpar
              // Radiation terms
              - alpha0 * Bsq * gamma * p * one_m_xisq * (1.0 + xi * b_rad_term);

    /**  \brief Pitch evolution
     *   \f{eqnarray*}{
     *     \frac {d\xi}{dt} &=& -\frac {c\tau_c}a \frac {p(1 - \xi^2)} {2\gamma
     * B^\ast_\parallel} B^\ast \cdot \nabla \ln B \\
     *                      &-&  \frac{c}{a\omega_{ce,0}} \frac{\xi (1 -
     * \xi^2)}{2B_\parallel^\ast} \nabla \ln B \cdot (\vec E\times \vec b) \\
     *                      &-& \frac{(1 - \xi^2)}{p B_\parallel} B^\ast \cdot E
     * \\
     *                    &+& \alpha_0 B^2 \xi \frac{(1-\xi^2)}{\gamma} \left(1
     * +
     * \xi \text{rad_term}\right)\\
     *                    &-& \alpha_0 B^2 \frac {\gamma (1-\xi^2)^2}2
     * \text{rad_term}
     *   \f}
     *
     */
    dXdt[1] = -ct_a * p * one_m_xisq * 0.5 / gamma / Bpar * Bstar_d_gradlnB -
              c_aw0 * xi * one_m_xisq * 0.5 / Bpar * gradlnB_d_Exb -
              one_m_xisq / p / Bpar * Bstar_d_E
              // Radiation terms
              +
              alpha0 * Bsq * xi * one_m_xisq / gamma * (1.0 + xi * b_rad_term) -
              alpha0 * Bsq * 0.5 * gamma * one_m_xisq * one_m_xisq * b_rad_term;

    /**  \brief Position evolution
     *   \f{eqnarray*}{
     *      \frac {dX}{dt} = \frac{c \tau_c}{a} \frac {\xi p}{\gamma}
     * \frac{B^\ast}{B^\ast_\parallel}
     *   - \frac{c\tau_c}a \frac{c}{a\omega_{ce,0}} \frac{p^2 (1-\xi^2)}{2\gamma
     * B_\parallel^\ast} \vec b \times \nabla \ln B
     *      + \frac{c}{a\omega_{ce,0}}\frac{E\times b}{B_\parallel^\ast}
     *   \f}
     *
     */
    for (int i = 0; i < 3; ++i)
      dXdt[i + 2] = ct_a * xi * p / gamma * Bstar[i] / Bpar -
                    0.5 * ct_a * c_aw0 / Bpar * one_m_xisq * p * p / gamma *
                        b_x_gradlnB[i] +
                    c_aw0 / Bpar * E_x_b[i];
    dXdt[3] = 0.0;

  };


  KOKKOS_INLINE_FUNCTION
  void computeConservedQuantities(Real &p_phi, Real &mu,
                                  const Real p, const Real xi,
                                  const Real R, const Real Z,
                                  const Real &t, const Real Psi) const {
    // Evaluate B and E at the particle position (EvalBE carries both).
    Kinetic::EvalBE ev;
    field.eval(ev, R, Z, t);
    if constexpr (EF == false) ev.E = {};

    const Real Bmag = ev.Bmag;

    // First adiabatic invariant mu = p_perp^2 / |B|.
    mu = p * p * (1.0 - xi * xi) / Bmag;

    // Canonical toroidal momentum:
    //   p_phi = c_aw0 * p_parallel * B_phi/|B| * R - Psi
    // ev.B[1] is B_phi. The '- Psi' is the (normalized) charge*poloidal-flux term.
    p_phi = c_aw0 * xi * p * ev.B[1] / Bmag * R - Psi;

    // When an accelerating E-field is present, p_phi is not conserved unless we
    // account for the toroidal impulse. ev.E[1] is E_phi.
    if constexpr (EF) p_phi += t * c_aw0 * ev.E[1] * R;
  };
};
