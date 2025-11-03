# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *


class Hflux(CMakePackage):
    """hFlux is lightweight toolkit for tokamak simulation code diagnostics."""

    homepage = "https://obeznosov-lanl.github.io/hFlux/"
    git = "https://github.com/obeznosov-LANL/hFlux.git"
    url = "https://github.com/obeznosov-LANL/hFlux/archive/refs/tags/v0.1.0-rc0.tar.gz"

    maintainers("tukss", "obeznosov-LANL")

    license("BSD-3-Clause", checked_by="tukss")

    version("main", branch="main")
    version("0.1.0-rc0", sha256="9c71e6e8b2a6ed18b3e7106fbb027981c6e398d2ade8381bf1b47021ad3d1ffb")

    depends_on("hdf5+cxx")
    depends_on("kokkos")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
