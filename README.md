# Hymera
Hymera: Hybrid MHD-kinetic code for relativistic particles

## Building the code
### Spack
Spack is a package manager tailored to the needs of high-performance computing users. It is not necessary to build hymera but it is a convenient way to install all the dependencies. This is a simple example of how to build an environment for hymera. More details can be found in the [Spack documentation](https://spack.readthedocs.io/).

```sh
# Get the latest release of Spack.
git clone -b releases/latest --depth 2 https://github.com/spack/spack
source spack/share/spack/setup-env.sh

# Create a new environment for hymera and activate it. The name is arbitrary.
spack env create hymera
spack env activate hymera

# Change to the directory containing this readme file

# Add the custom package repository containing hflux, hymera, and an updated Parthenon.
spack repo add spack-repo

# Add the hymera package to the list of packages to be installed.
spack add hymera

# Add cmake to the installed packages so we can use it directly later.
spack add cmake

# Tell spack to build hymera from the local directory instead of from the repository.
# This is needed until hymera moves to its own repository.
spack develop -p "$PWD" hymera@main

# Load all the modules for MPI, compilers etc. you would like to use and add them as externals.
# This is just a suggestion. The exact packages depend on your local configuration.
spack compiler find
spack external find openssl openssh curl
spack external find openmpi

# Concretize the packages to be installed. Check if the list is reasonable. Possibly add more external packages.
spack concretize -f

# Install the packages if you're fine with the list.
# This doesn't install the code itself but all the dependencies so we can build it manually with cmake.
spack install --only dependencies

# Now build hymera itself
cmake -Bbuild
cmake --build build -j $(nproc)
```

To reactivate the same environment in a later session:
```sh
source spack/share/spack/setup-env.sh
spack env activate hymera
```



# Release
O5009 hFlux was approved for Open-Source Assertion

