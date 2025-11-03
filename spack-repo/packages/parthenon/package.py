from spack.package import *
from spack_repo.builtin.packages.parthenon.package import Parthenon

class Parthenon(Parthenon):
    version("25.05", sha256="bb1b4fa0c7018cb63b162065d9fceee4229c8efb5f0adcd2232f92d50d775fe8")
    depends_on("c", type="build")
