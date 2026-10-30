#include "equilibrium.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void write_vector(const std::filesystem::path& path,const std::vector<double>& values) {
    std::ofstream out(path);
    if(!out) throw std::runtime_error("Cannot create output file");
    out << values.size() << '\n' << std::setprecision(17);
    for(double x:values) out << x << '\n';
    if(!out) throw std::runtime_error("Failed writing output");
}
int integer(const char* s) {
    std::size_t used=0;const std::string text(s);int n=std::stoi(text,&used);
    if(used!=text.size()) throw std::invalid_argument("Invalid integer");
    return n;
}
}
int main(int argc,char** argv) {
    try {
        if(argc==4 && std::string(argv[1])=="points") {
            std::ifstream in(argv[2]);
            if(!in) throw std::runtime_error("Cannot open input coordinates");
            if(std::filesystem::exists(argv[3])) throw std::runtime_error("Output already exists");
            std::ofstream out(argv[3]);
            if(!out) throw std::runtime_error("Cannot create output");
            out << "# R_m Z_m psi_Vs g_Tm psi_Hymera_Vs\n" << std::setprecision(17);
            double r,z;
            while(in>>r) {
                if(!(in>>z)) throw std::runtime_error("Incomplete coordinate pair");
                auto value=iter_pinn::evaluate(r,z);
                out << r << ' ' << z << ' ' << value.psi << ' ' << value.g << ' ' << -value.psi << '\n';
            }
            if(!in.eof()) throw std::runtime_error("Invalid coordinate input");
            if(!out) throw std::runtime_error("Failed writing output");
            return 0;
        }
        if((argc==6||argc==11) && std::string(argv[1])=="grid") {
            const std::filesystem::path dir(argv[2]);
            if(std::filesystem::exists(dir)) throw std::runtime_error("Output directory already exists");
            iter_pinn::Grid grid{integer(argv[3]),integer(argv[4]),integer(argv[5]),1.525,4.975,-2.975,2.975,2};
            if(argc==11) {
                grid.rmin=std::stod(argv[6]); grid.rmax=std::stod(argv[7]);
                grid.zmin=std::stod(argv[8]); grid.zmax=std::stod(argv[9]); grid.length_scale_m=std::stod(argv[10]);
            }
            const auto result=iter_pinn::evaluate_grid(grid);
            std::filesystem::create_directories(dir);
            // Match the exact existing count-header, R-fastest ASCII format.
            char suffix[100];
            std::snprintf(suffix,sizeof(suffix),"_grid%03dx%02dx%03d.txt",grid.nr,grid.nphi,grid.nz);
            write_vector(dir/("vecpsi"+std::string(suffix)),result.psi);
            write_vector(dir/("vecg"+std::string(suffix)),result.g);
            std::ofstream metadata(dir/"grid.json");
            metadata << std::setprecision(17) << "{\"nr\":" << grid.nr << ",\"nphi\":" << grid.nphi
                     << ",\"nz\":" << grid.nz << ",\"normalized_bounds\":[" << grid.rmin << ',' << grid.rmax
                     << ',' << grid.zmin << ',' << grid.zmax << "],\"length_scale_m\":" << grid.length_scale_m
                     << ",\"psi_sign_relative_to_PINN\":-1,\"g_sign_relative_to_PINN\":1,\"psi_offset\":0}\n";
            std::cout << "Wrote " << result.psi.size() << " psi and " << result.g.size() << " g values\n";
            return 0;
        }
        std::cerr << "Usage: iter_pinn_evaluate points coordinates.txt values.txt\n"
                     "       iter_pinn_evaluate grid NEW_DIRECTORY NR NPHI NZ [RMIN RMAX ZMIN ZMAX L0_M]\n"
                     "Grid bounds are normalized by L0_M; point coordinates are metres.\n";
        return 2;
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
