#include "equilibrium.hpp"
#include "model_data.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace iter_pinn {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

// Same eight AGM iterations as the executed JAX Green function. The parameter
// m is k^2, rather than the elliptic modulus k. No clipping of singular values.
double green(double r, double z, double a, double b) {
    const double m = 4*r*a/((r+a)*(r+a)+(z-b)*(z-b));
    if (!(m > 0 && m < 1)) throw std::runtime_error("Singular coil quadrature point");
    double aa=1, bb=std::sqrt(1-m), total=m/2, factor=1;
    for (int k=0;k<8;++k) {
        const double c=(aa-bb)/2, next=(aa+bb)/2;
        total += factor*c*c;
        bb=std::sqrt(aa*bb); aa=next; factor*=2;
    }
    const double K=pi/(2*aa), E=K*(1-total);
    return 2e-7*std::sqrt(r*a/m)*((2-m)*K-2*E);
}

double coil_flux(double r, double z, int coil) {
    const double* rect=data::rectangles+4*coil;
    const double r0=rect[0], r1=rect[1], z0=rect[2], z1=rect[3];
    const double ro=std::clamp(r,r0,r1), zo=std::clamp(z,z0,z1);
    double result=0;
    if (std::hypot(r-ro,z-zo) >= .5*std::max(r1-r0,z1-z0)) {
        for(int j=0;j<10;++j) for(int i=0;i<10;++i)
            result+=data::w10[j]*data::w10[i]*green(r,z,
                         r0+(r1-r0)*data::x10[i],z0+(z1-z0)*data::x10[j]);
    } else {
        const double corners[4][2]={{r0,z0},{r1,z0},{r1,z1},{r0,z1}};
        for(int k=0;k<4;++k) {
            const double ar=corners[k][0]-ro, az=corners[k][1]-zo;
            const double br=corners[(k+1)%4][0]-ro, bz=corners[(k+1)%4][1]-zo;
            const double jac=std::abs(ar*bz-az*br);
            if(jac==0) continue;
            for(int i=0;i<32;++i) for(int j=0;j<32;++j) {
                const double u=data::x32[i], v=data::x32[j];
                result+=jac*data::w32[i]*data::w32[j]*u*green(r,z,
                    ro+u*((1-v)*ar+v*br),zo+u*((1-v)*az+v*bz))/((r1-r0)*(z1-z0));
            }
        }
    }
    return result*data::currents[coil];
}

double self_flux(double r,double z) {
    const int widths[]={2,32,32,32,32,1};
    std::array<double,32> x{}, y{};
    x[0]=(r-6.1)/3.1; x[1]=(z-.05)/5.35;
    int offset=0;
    for(int layer=0;layer<5;++layer) {
        const int nin=widths[layer], nout=widths[layer+1];
        for(int i=0;i<nout;++i) {
            double value=data::network[offset+nin*nout+i];
            // Archived weights were flattened in Fortran order: output fastest.
            for(int j=0;j<nin;++j) value+=data::network[offset+j*nout+i]*x[j];
            y[i]=layer==4 ? value : std::tanh(value);
        }
        offset+=(nin+1)*nout; x=y;
    }
    const double dr=(r-6.2)/3.5, dz=(z-.5)/5.5;
    return 20*(1.5*std::exp(-dr*dr-dz*dz)+.1*x[0]);
}

bool plasma_support(double r,double z) {
    // Frozen axis-connected source selector, exactly as saved in the model.
    const int ir=std::clamp(static_cast<int>(std::floor((r-3)/(6.2/80))),0,79);
    const int iz=std::clamp(static_cast<int>(std::floor((z+5.3)/(10.7/128))),0,127);
    if(!data::support[iz*80+ir]) return false;
    bool inside=false;
    constexpr int n=sizeof(data::wall)/sizeof(double)/2;
    for(int i=0,j=n-1;i<n;j=i++) {
        const double xi=data::wall[2*i], yi=data::wall[2*i+1];
        const double xj=data::wall[2*j], yj=data::wall[2*j+1];
        if ((yi>z)!=(yj>z) && r<(xj-xi)*(z-yi)/(yj-yi)+xi) inside=!inside;
    }
    return inside;
}
void check_adapter(double scale,double offset,int sign) {
    if(!std::isfinite(scale)||scale<=0||!std::isfinite(offset)||(sign!=1 && sign!=-1))
        throw std::invalid_argument("Invalid length scale, flux offset, or g sign");
}
} // namespace

Field evaluate(double r,double z) {
    if(!std::isfinite(r)||!std::isfinite(z)||r<.5||r>14||z< -9||z>9)
        throw std::out_of_range("Point is outside the trained ITER domain (metres)");
    double psi=self_flux(r,z);
    for(int i=0;i<11;++i) psi+=coil_flux(r,z,i);
    double g=32.86;
    const double pa=data::critical_flux[0], pb=data::critical_flux[1];
    if(psi>pb && plasma_support(r,z)) {
        const double s=(pa-psi)/(pa-pb), width=.15;
        const double x=std::clamp((1-s)/width,0.,1.);
        const double primitive=s>=1 ? 0 : (s<=1-width ? s-1+width/2 : -width*(x*x*x-.5*x*x*x*x));
        g-=data::alpha*(pa-pb)*primitive;
    }
    return {psi,g};
}
Field evaluate_hymera(double r,double z,double scale,double offset,int sign) {
    check_adapter(scale,offset,sign);
    const Field value=evaluate(r*scale,z*scale);
    return {-value.psi+offset,sign*value.g};
}
GridValues evaluate_grid(const Grid& grid,double offset,int sign) {
    check_adapter(grid.length_scale_m,offset,sign);
    if(grid.nr<=0||grid.nphi<=0||grid.nz<=0||!std::isfinite(grid.rmin)||
       !std::isfinite(grid.rmax)||!std::isfinite(grid.zmin)||!std::isfinite(grid.zmax)||
       grid.rmax<=grid.rmin||grid.zmax<=grid.zmin)
        throw std::invalid_argument("Invalid grid dimensions or bounds");
    // Preflight both corners before allocation; evaluate() enforces full domain.
    evaluate_hymera(grid.rmin,grid.zmin,grid.length_scale_m,offset,sign);
    evaluate_hymera(grid.rmax,grid.zmax,grid.length_scale_m,offset,sign);
    const std::size_t nr=grid.nr, np=grid.nphi, nz=grid.nz;
    if(nr>1000000 || np>1000000 || nz>1000000 ||
       (nr+1)*np>std::numeric_limits<std::size_t>::max()/(nz+1) ||
       nr*(np+1)>std::numeric_limits<std::size_t>::max()/nz)
        throw std::length_error("Grid size overflow");
    GridValues result;
    result.psi.resize((nr+1)*np*(nz+1)); result.g.resize(nr*(np+1)*nz);
    const double dr=(grid.rmax-grid.rmin)/grid.nr, dz=(grid.zmax-grid.zmin)/grid.nz;
    for(std::size_t k=0;k<=nz;++k) for(std::size_t i=0;i<=nr;++i) {
        const double value=evaluate_hymera(grid.rmin+i*dr,grid.zmin+k*dz,
                           grid.length_scale_m,offset,sign).psi;
        for(std::size_t j=0;j<np;++j) result.psi[i+j*(nr+1)+k*np*(nr+1)]=value;
    }
    for(std::size_t k=0;k<nz;++k) for(std::size_t i=0;i<nr;++i) {
        const double value=evaluate_hymera(grid.rmin+(i+.5)*dr,grid.zmin+(k+.5)*dz,
                           grid.length_scale_m,offset,sign).g;
        for(std::size_t j=0;j<=np;++j) result.g[i+j*nr+k*(np+1)*nr]=value;
    }
    return result;
}
} // namespace iter_pinn
