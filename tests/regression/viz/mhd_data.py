"""Readers for the MHD solver's on-disk data formats.

Three formats are covered:

1. The ASCII grid-data files in ``inputs/mhd/`` (``veceta``, ``vecpsi``, ``vecg``),
   read by ``ReadInitialData`` (src/mhd/ts_functions.c) and indexed in C as
   ``data[er + ephi*N0 + ez*N1*N0]`` with a per-file stagger; see FIELD_LAYOUTS.
2. ``AxisSymmetricGeometry.dat``, a tab-separated Nr x Nz integer matrix read by
   the C++ side at src/kinetic/kinetic.cpp.
3. PETSc binary Vec files, either a single Vec (``inputs/X_16.petscbin``) or the
   11-record sequence written by ``stag_vec_io`` (src/mhd/mhd.c).

The stag_vec_io format has no header, no magic, no version, and no dimension
check, and ``PetscObjectSetName`` is ignored by the binary viewer -- so records
are matched purely by loop order. That order is reproduced in STAG_RECORDS and
must stay in sync with the ``locs[]`` array in ``stag_vec_io``.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

# PETSc's VEC_FILE_CLASSID, written big-endian at the head of every Vec record.
VEC_FILE_CLASSID = 1211214

# DMStagStencilLocation values for the eight 3-D strata, in the order
# stag_vec_io iterates them. The trailing int is the DOF count per stratum for
# the production layout dof0=4, dof1=dof2=dof3=1 (src/mhd/mhd.c).
STAG_RECORDS = [
    ("BACK_DOWN_LEFT", 4),  # vertex: V_r, V_phi, V_z, EP
    ("BACK_DOWN", 1),       # edge
    ("BACK_LEFT", 1),       # edge
    ("DOWN_LEFT", 1),       # edge
    ("BACK", 1),            # face
    ("DOWN", 1),            # face
    ("LEFT", 1),            # face
    ("ELEMENT", 1),         # cell centre: n_i
]

# Component names for the 4-DOF vertex stratum.
VERTEX_COMPONENTS = ["V_r", "V_phi", "V_z", "EP"]

# Material region tags in veceta / AxisSymmetricGeometry.dat, with the
# resistivity each selects in alphaec2 (src/mhd/mass_matrix_coefficients.c).
MATERIAL_TAGS = {
    1: ("plasma", "etaplasma"),
    2: ("separatrix-wall", "etasepwal"),
    0: ("blanket wall", "etawall"),
    -1: ("vacuum vessel", "etaVV"),
    -2: ("exterior", "etaout"),
}

# The predicate fabs(dataC - 1.5) < 0.7, which appears at 112 sites in
# ts_functions.c, selects exactly these tags.
PLASMA_REGION_TAGS = (1, 2)

# Per-file stagger of the ASCII grid data. Each entry gives the array shape as
# a function of (Nr, Nphi, Nz), matching the C indexing arithmetic.
FIELD_LAYOUTS = {
    # dataC[er + ephi*Nr + ez*Nphi*Nr]
    "veceta": lambda nr, nphi, nz: (nr, nphi, nz),
    # datapsi[er + ephi*(Nr+1) + ez*Nphi*(Nr+1)]  -- r-z edge, so +1 in r and z
    "vecpsi": lambda nr, nphi, nz: (nr + 1, nphi, nz + 1),
    # datag[er + ephi*Nr + ez*Nr*(Nphi+1)]  -- phi-face, so +1 in phi
    "vecg": lambda nr, nphi, nz: (nr, nphi + 1, nz),
}


def expected_count(kind: str, nr: int, nphi: int, nz: int) -> int:
    """Element count an ASCII grid-data file must contain for a given grid."""
    shape = FIELD_LAYOUTS[kind](nr, nphi, nz)
    return shape[0] * shape[1] * shape[2]


def read_ascii_grid(path, kind: str, nr: int, nphi: int, nz: int) -> np.ndarray:
    """Read one inputs/mhd/*.txt file into an (er, ephi, ez) array.

    The file's first token is a self-declared count. ``ReadInitialData`` trusts
    it without checking; we verify it against both the real token count and the
    grid formula, since a mismatch means the grid and the data disagree.
    """
    tokens = Path(path).read_text().split()
    declared = int(tokens[0])
    values = np.array(tokens[1:], dtype=np.float64)

    if declared != values.size:
        raise ValueError(
            f"{path}: header declares {declared} values but file holds {values.size}"
        )

    shape = FIELD_LAYOUTS[kind](nr, nphi, nz)
    want = shape[0] * shape[1] * shape[2]
    if values.size != want:
        raise ValueError(
            f"{path}: holds {values.size} values, but a {nr}x{nphi}x{nz} grid "
            f"needs {want} for '{kind}' (shape {shape})"
        )

    # C index is er + ephi*stride_r + ez*stride_r*stride_phi, i.e. er fastest
    # and ez slowest, so read in Fortran order into (er, ephi, ez).
    return values.reshape(shape, order="F")


def read_geometry_dat(path) -> np.ndarray:
    """Read AxisSymmetricGeometry.dat into an (er, ez) array of material tags."""
    rows = [line.split() for line in Path(path).read_text().splitlines() if line.strip()]
    return np.array(rows, dtype=np.float64)


def read_petsc_vecs(path) -> list[np.ndarray]:
    """Read a PETSc binary file as a sequence of Vec records.

    Each record is: int32 classid, int32 length, then length float64 values,
    all big-endian. Returns one 1-D array per record.
    """
    raw = Path(path).read_bytes()
    vecs: list[np.ndarray] = []
    off = 0
    while off < len(raw):
        if off + 8 > len(raw):
            raise ValueError(f"{path}: truncated record header at byte {off}")
        classid, length = struct.unpack_from(">ii", raw, off)
        off += 8
        if classid != VEC_FILE_CLASSID:
            raise ValueError(
                f"{path}: expected Vec classid {VEC_FILE_CLASSID} at byte "
                f"{off - 8}, found {classid}. Not a PETSc Vec file, or the "
                f"record sequence is out of step."
            )
        nbytes = length * 8
        if off + nbytes > len(raw):
            raise ValueError(f"{path}: record of {length} values is truncated")
        vecs.append(
            np.frombuffer(raw, dtype=">f8", count=length, offset=off).astype(np.float64)
        )
        off += nbytes
    return vecs


def plasma_mask(tags: np.ndarray) -> np.ndarray:
    """The predicate fabs(dataC - 1.5) < 0.7 from ts_functions.c, vectorized."""
    return np.abs(tags - 1.5) < 0.7
