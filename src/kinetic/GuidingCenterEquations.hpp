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
  ERROR_CODE operator()(const Real &t, const value_type &X,
                        value_type &dXdt) const {

    const Real p =   X[0];
    const Real xi =  X[1];
    const Real R =   X[2];
    const Real phi = X[3];
    const Real Z =   X[4];
    const Real gamma = std::sqrt(p * p + 1.0);

    if constexpr (SlabModel) {
      dXdt[0] = -xi * field.E_0 - alpha0 * p * gamma * (1.0 - xi * xi);
      dXdt[1] = -(1.0 - xi * xi) / p * field.E_0 + alpha0 * xi * (1.0 - xi * xi) / gamma;
      dXdt[2] = 0.0;
      dXdt[3] = 0.0;
      dXdt[4] = 0.0;
      return SUCCESS;
    }


    Dim3 vB = {}, dBdR = {}, dBdZ = {}, curlB = {}, E = {};

    ERROR_CODE status = field(X, t, vB, curlB, dBdR, dBdZ, E);
    if constexpr (EF == false) E = {};
    if (status != ERROR_CODE::SUCCESS)
      return status;


    Real Bsq = dot_product(vB, vB);
    Real B = std::sqrt(Bsq);
    Dim3 b = {};
    for (int i = 0; i < 3; ++i)
      b[i] = vB[i] / B;

    Dim3 gradlnB = {dot_product(vB, dBdR) / Bsq, 0.0,
                    dot_product(vB, dBdZ) / Bsq};

    Dim3 b_x_gradlnB;
    cross_product(b, gradlnB, b_x_gradlnB);

    Dim3 Bstar;
    for (int i = 0; i < 3; ++i)
      Bstar[i] = vB[i] - c_aw0 * p * xi * (b_x_gradlnB[i] + curlB[i] / B);

    Real Bpar = dot_product(b, Bstar);
    Real B_d_gradlnB = dot_product(vB, gradlnB);
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
    cross_product(E, b, E_x_b);

    Real gradlnB_d_Exb = dot_product(gradlnB, E_x_b);
    Real E_d_curlB = dot_product(E, curlB);
    Real B_d_E = dot_product(vB, E);
    Real Bstar_d_E = dot_product(Bstar, E);

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
    for (int i = 0; i < 3; i += 2)
      dXdt[i + 2] = ct_a * xi * p / gamma * Bstar[i] / Bpar -
                    0.5 * ct_a * c_aw0 / Bpar * one_m_xisq * p * p / gamma *
                        b_x_gradlnB[i] +
                    c_aw0 / Bpar * E_x_b[i];
    dXdt[3] = 0.0;

    // DEBUG("X = (%le %le %le %le %le), dX = (%le %le %le %le %le), t = %le",
    // X[0], X[1], X[2], X[3], X[4] ,dXdt[0], dXdt[1], dXdt[2], dXdt[3],
    // dXdt[4], t);
    return SUCCESS;
  };

  template<typename ViewType>
  KOKKOS_INLINE_FUNCTION
  ERROR_CODE operator()(const Real &t, const int idx, const ViewType& v,
                        value_type &dXdt) const {

    const Real p =   v(idx, 0);
    const Real xi =  v(idx, 1);
    const Real R =   v(idx, 2);
    const Real phi = v(idx, 3);
    const Real Z =   v(idx, 4);
    const Real gamma = std::sqrt(p * p + 1.0);

    if constexpr (SlabModel) {
      dXdt[0] = -xi * field.E_0 - alpha0 * p * gamma * (1.0 - xi * xi);
      dXdt[1] = -(1.0 - xi * xi) / p * field.E_0 + alpha0 * xi * (1.0 - xi * xi) / gamma;
      dXdt[2] = 0.0;
      dXdt[3] = 0.0;
      dXdt[4] = 0.0;
      return SUCCESS;
    }

    Dim3 vB = {}, dBdR = {}, dBdZ = {}, curlB = {}, E = {};

    auto status = field.checkWall(idx, v);
    if (status != ERROR_CODE::SUCCESS)
      return status;

    status = field(idx, v, t, vB, curlB, dBdR, dBdZ, E);
    if (status != ERROR_CODE::SUCCESS)
      return status;

    Real Bsq = dot_product(vB, vB);
    Real B = std::sqrt(Bsq);
    Dim3 b;
    for (int i = 0; i < 3; ++i)
      b[i] = vB[i] / B;

    //if constexpr (Flags::SlabGeometry) {
//      Real b_d_E = dot_product(b, E);
//      dXdt[0] = xi * b_d_E - alpha0 * p * gamma * (1.0 - xi * xi);
//      dXdt[1] =
//          -(1.0 - xi * xi) / p * b_d_E + alpha0 * xi * (1.0 - xi * xi) / gamma;
//      dXdt[2] = 0.0;
//      dXdt[3] = 0.0;
//      dXdt[4] = 0.0;
//      return SUCCESS;

    Dim3 gradlnB = {dot_product(vB, dBdR) / Bsq, 0.0,
                    dot_product(vB, dBdZ) / Bsq};

    Dim3 b_x_gradlnB;
    cross_product(b, gradlnB, b_x_gradlnB);

    Dim3 Bstar;
    for (int i = 0; i < 3; ++i)
      Bstar[i] = vB[i] - c_aw0 * p * xi * (b_x_gradlnB[i] + curlB[i] / B);

    Real Bpar = dot_product(b, Bstar);
    Real B_d_gradlnB = dot_product(vB, gradlnB);
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
    cross_product(E, b, E_x_b);

    Real gradlnB_d_Exb = dot_product(gradlnB, E_x_b);
    Real E_d_curlB = dot_product(E, curlB);
    Real B_d_E = dot_product(vB, E);
    Real Bstar_d_E = dot_product(Bstar, E);

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
    for (int i = 0; i < 3; i += 2)
      dXdt[i + 2] = ct_a * xi * p / gamma * Bstar[i] / Bpar -
                    0.5 * ct_a * c_aw0 / Bpar * one_m_xisq * p * p / gamma *
                        b_x_gradlnB[i] +
                    c_aw0 / Bpar * E_x_b[i];
    dXdt[3] = 0.0;

    // DEBUG("X = (%le %le %le %le %le), dX = (%le %le %le %le %le), t = %le",
    // X[0], X[1], X[2], X[3], X[4] ,dXdt[0], dXdt[1], dXdt[2], dXdt[3],
    // dXdt[4], t);
    return SUCCESS;
  };

  KOKKOS_INLINE_FUNCTION
  void computeConservedQuantities(const Dim5 &X, Real &p_phi, Real &mu,
                                  const Real &t) const {
    // Compute magnitude B
    Dim3 B = {}, curlB = {}, dBdR = {}, dBdZ = {}, E = {};
    field(X, t, B, curlB, dBdR, dBdZ, E);
    // if constexpr (EF == false) E = {};

    const Real Bmag = sqrt(dot_product(B, B));

    // Get variables
    const Real &p = X[0];
    const Real &xi = X[1];
    const Real &R = X[2];
    const Real &Z = X[4];

    // compute mu = p_\perp / |B|
    mu = p * p * (1 - xi * xi) / Bmag;

    Real Psi = field.Psi(X);
    // p_phi = (c over a omega) p_perp B_phi R + (-) psi
    // - stands for charge p is gamma m_e c
    // p_phi += t * c_a omega_0 * E_phi * R  if E is present
    p_phi = c_aw0 * xi * p * B[1] / Bmag * R - Psi;
    p_phi += t * c_aw0 * E[1] * field.R_a;
  };

  template<typename ViewType>
  KOKKOS_INLINE_FUNCTION
  void computeConservedQuantities(const int idx, const ViewType& v, Real &p_phi, Real &mu,
                                  const Real &t) const {
    // Compute magnitude B
    Dim3 B = {}, curlB = {}, dBdR = {}, dBdZ = {}, E = {};
    field(idx, v, t, B, curlB, dBdR, dBdZ, E);
    if constexpr (EF == false) E = {};

    const Real Bmag = sqrt(dot_product(B, B));

    // Get variables
    const Real p =  v(idx, 0);
    const Real xi = v(idx, 1);
    const Real R =  v(idx, 2);
    const Real Z =  v(idx, 4);

    // compute mu = p_\perp / |B|
    mu = p * p * (1 - xi * xi) / Bmag;

    Real Psi = field.Psi(idx, v, t);
    // p_phi = (c over a omega) p_perp B_phi R + (-) psi
    // - stands for charge p is gamma m_e c
    // p_phi += t * c_a omega_0 * E_phi * R  if E is present
    p_phi = c_aw0 * xi * p * B[1] / Bmag * R - Psi;
    // p_phi += t * c_aw0 * field.E_0 * field.R_a;
  };
};
