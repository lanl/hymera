# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *


class Hymera(CMakePackage):
    """FIXME: Put a proper description of your package here."""

    # FIXME: Add a proper url for your package's homepage here.
    #url = "https://www.example.com/example-1.2.3.tar.gz"

    maintainers("tukss", "obeznosov-LANL")

    # FIXME: Add the SPDX identifier of the project's license below.
    # See https://spdx.org/licenses/ for a list. Upon manually verifying
    # the license, set checked_by to your Github username.
    #license("UNKNOWN", checked_by="github_user1")

    version("main", branch="hybrid_main")

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
