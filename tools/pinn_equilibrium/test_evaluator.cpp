#include "equilibrium.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool condition) {if(!condition) throw std::runtime_error("Evaluator test failed");}
void close(double a,double b) {require(std::abs(a-b)<1e-10);}
int main() {
    using namespace iter_pinn;
    const Field native=evaluate(6.2,.6), mapped=evaluate_hymera(3.1,.3);
    close(mapped.psi,-native.psi);close(mapped.g,native.g);
    const Field adjusted=evaluate_hymera(3.1,.3,2,100,-1);
    close(adjusted.psi,100-native.psi);close(adjusted.g,-native.g);
    for(auto p: {std::pair<double,double>{0,0},{14.1,0},{6,9.1},
                 {std::numeric_limits<double>::quiet_NaN(),0}}) {
        bool threw=false;try{evaluate(p.first,p.second);}catch(const std::out_of_range&){threw=true;}
        require(threw);
    }
    bool threw=false;try{evaluate_hymera(3,0,0);}catch(const std::invalid_argument&){threw=true;}
    require(threw);
    const Grid grid{3,2,4,2.5,3.5,-.5,.5,2};
    const auto result=evaluate_grid(grid,7,-1);
    require(result.psi.size()==40 && result.g.size()==36);
    for(int k=0;k<=4;++k)for(int j=0;j<2;++j)for(int i=0;i<=3;++i)
        close(result.psi[i+j*4+k*8],evaluate_hymera(2.5+i/3.,-.5+k*.25,2,7,-1).psi);
    for(int k=0;k<4;++k)for(int j=0;j<=2;++j)for(int i=0;i<3;++i)
        close(result.g[i+j*3+k*9],evaluate_hymera(2.5+(i+.5)/3.,-.5+(k+.5)*.25,2,7,-1).g);
    for(auto p: {std::pair<double,double>{.5,-9},{14,9},{.5,9},{14,-9}}) {
        const auto f=evaluate(p.first,p.second);require(std::isfinite(f.psi)&&std::isfinite(f.g));
    }
    std::cout << "Unit, adapter, bounds and staggered-grid checks passed\n";
}
