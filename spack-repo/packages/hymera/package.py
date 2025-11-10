# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *


class Hymera(CMakePackage):
    """Hymera is a hybrid MHD-kinetic code for relativistic particles"""

    # FIXME: Add a proper url for your package's homepage here.
    #url = "https://www.example.com/example-1.2.3.tar.gz"
    git = "https://github.com/lanl/hymera.git"

    maintainers("tukss", "obeznosov-LANL")

    license("BSD-3-Clause", checked_by="tukss")
    
    version("main", branch="main")

    depends_on("parthenon@25.05:")
    depends_on("hdf5+cxx")
    depends_on("kokkos")
    depends_on("petsc")
    depends_on("hflux@main")

    depends_on("c", type="build")
    depends_on("cxx", type="build")

    def cmake_args(self):
        args = [self.define("ENABLE_INTERNAL_PARTHENON", False)]
        return args
