#pragma once
#include <vector>

namespace iter_pinn {
// psi is poloidal flux / (2*pi), in Wb/rad (numerically V s); g = R B_phi in T m.
// Native convention: B_R = -psi_Z/R, B_Z = psi_R/R, B_phi = g/R.
struct Field { double psi, g; };
struct Grid {
    int nr, nphi, nz;
    double rmin, rmax, zmin, zmax; // normalized by length_scale_m
    double length_scale_m = 2.0;
};
struct GridValues {
    // R fastest, phi next, Z slowest. psi has (nr+1)*nphi*(nz+1) entries;
    // g has nr*(nphi+1)*nz entries, including the repeated periodic end face.
    std::vector<double> psi, g;
};
// Throws for nonfinite coordinates or points outside R=[0.5,14], Z=[-9,9] m.
Field evaluate(double r_m, double z_m);
// Takes Hymera-normalized coordinates; returns dimensional fields for its loader.
// psi sign compensates Hymera's B_p = -curl(psi/R e_phi). g sign defaults to
// preserving the PINN's physical toroidal field; -1 explicitly reverses that field.
Field evaluate_hymera(double r, double z, double length_scale_m = 2.0,
                     double psi_offset = 0.0, int g_sign = 1);
GridValues evaluate_grid(const Grid& grid, double psi_offset = 0.0, int g_sign = 1);
} // namespace iter_pinn
