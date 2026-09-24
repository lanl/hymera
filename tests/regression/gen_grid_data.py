#!/usr/bin/env python
"""Generate inputs/mhd grid-data files at an arbitrary resolution.

The solver reads three ASCII files whose names encode the grid
(src/mhd/mhd.c:100-107), so changing NR/Nphi/NZ requires generating them. Only
`veceta` is needed for the manufactured initial conditions (ictype 1..8), since
those never read psi or g -- but veceta is read at 139 sites in
mass_matrix_coefficients.c and is therefore mandatory at every resolution.

Two modes:

  --mode uniform   every cell tagged plasma (1). This is what the manufactured
                   convergence tests want: a single material, so the mass-matrix
                   coefficients are smooth and the measured convergence order
                   reflects the discretization rather than material jumps.

  --mode downsample
                   nearest-neighbour resample of the production 100x200 tag map
                   (equivalently inputs/AxisSymmetricGeometry.dat, which is the
                   same data transposed). Preserves the tokamak layering, so the
                   code paths for all five materials stay exercised.
                   NOTE: this is NOT a physics-equivalent coarsening. The
                   hardcoded cell indices in betaephi_isolcell
                   (mass_matrix_coefficients.c:3426,3451,3469) address absolute
                   (er,ez) positions and will land on different cells, or none.

psi/g are emitted only with --with-equilibrium, by bilinear interpolation, and
are likewise NOT a discrete Grad-Shafranov solution on the new mesh -- so do not
use them to certify physics. They exist only so that an ictype 9 run can be
started at a reduced size for smoke-testing.

Usage:
    gen_grid_data.py --nr 25 --nz 50 -o tests/regression/fixtures/mms_025x02x050
    gen_grid_data.py --nr 50 --nz 100 --mode downsample --with-equilibrium -o DIR
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "viz"))
from mhd_data import expected_count, read_ascii_grid  # noqa: E402

PROD_NR, PROD_NPHI, PROD_NZ = 100, 2, 200


def find_repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "src" / "mhd").is_dir():
            return parent
    raise SystemExit("cannot locate repo root (no src/mhd above this script)")


def write_field(path: Path, arr: np.ndarray, fmt: str = "%.10g") -> None:
    """Write in ReadInitialData's format: count on line 1, then one value per line.

    C reads with index er + ephi*s_r + ez*s_r*s_phi -- er fastest, ez slowest --
    so flatten in Fortran order.
    """
    flat = arr.flatten(order="F")
    with path.open("w") as fh:
        fh.write(f"{flat.size}\n")
        for v in flat:
            fh.write(f"{fmt % v}\n")
    print(f"wrote {path.name}: {flat.size} values, shape {arr.shape}")


def resample_nearest(src: np.ndarray, nr: int, nz: int) -> np.ndarray:
    """Nearest-neighbour resample of an (er, ez) array. Exact for tag fields."""
    ir = np.clip((np.arange(nr) + 0.5) * src.shape[0] / nr, 0, src.shape[0] - 1)
    iz = np.clip((np.arange(nz) + 0.5) * src.shape[1] / nz, 0, src.shape[1] - 1)
    return src[np.round(ir).astype(int)[:, None], np.round(iz).astype(int)[None, :]]


def resample_bilinear(src: np.ndarray, nr: int, nz: int) -> np.ndarray:
    """Bilinear resample of an (er, ez) array, for smooth fields only."""
    r = np.linspace(0, src.shape[0] - 1, nr)
    z = np.linspace(0, src.shape[1] - 1, nz)
    r0 = np.floor(r).astype(int)
    z0 = np.floor(z).astype(int)
    r1 = np.minimum(r0 + 1, src.shape[0] - 1)
    z1 = np.minimum(z0 + 1, src.shape[1] - 1)
    fr = (r - r0)[:, None]
    fz = (z - z0)[None, :]
    return (
        src[np.ix_(r0, z0)] * (1 - fr) * (1 - fz)
        + src[np.ix_(r1, z0)] * fr * (1 - fz)
        + src[np.ix_(r0, z1)] * (1 - fr) * fz
        + src[np.ix_(r1, z1)] * fr * fz
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nr", type=int, required=True)
    ap.add_argument("--nphi", type=int, default=2)
    ap.add_argument("--nz", type=int, required=True)
    ap.add_argument("--mode", choices=("uniform", "downsample"), default="uniform")
    ap.add_argument("--tag", type=int, default=1,
                    help="material tag for --mode uniform (1 = plasma)")
    ap.add_argument("--with-equilibrium", action="store_true",
                    help="also emit vecpsi/vecg (NOT a valid GS solution)")
    ap.add_argument("-o", "--outdir", required=True)
    args = ap.parse_args()

    root = find_repo_root()
    outdir = Path(args.outdir)
    if not outdir.is_absolute():
        outdir = root / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    nr, nphi, nz = args.nr, args.nphi, args.nz
    suffix = f"grid{nr:03d}x{nphi:02d}x{nz:03d}.txt"

    # --- veceta: material tags, cell-centred, Nr x Nphi x Nz ---------------
    if args.mode == "uniform":
        plane = np.full((nr, nz), float(args.tag))
        print(f"veceta: uniform tag {args.tag}")
    else:
        prod = root / "inputs" / "mhd" / (
            f"veceta_grid{PROD_NR:03d}x{PROD_NPHI:02d}x{PROD_NZ:03d}.txt"
        )
        src = read_ascii_grid(prod, "veceta", PROD_NR, PROD_NPHI, PROD_NZ)[:, 0, :]
        plane = resample_nearest(src, nr, nz)
        tags, counts = np.unique(plane, return_counts=True)
        print(f"veceta: downsampled from {src.shape} -> {plane.shape}, "
              f"tags {dict(zip(tags.astype(int).tolist(), counts.tolist()))}")

    # The problem is axisymmetric: replicate the single plane across phi.
    eta = np.repeat(plane[:, None, :], nphi, axis=1)
    write_field(outdir / f"veceta_{suffix}", eta, "%.1f")

    if args.with_equilibrium:
        prod_dir = root / "inputs" / "mhd"
        psuffix = f"grid{PROD_NR:03d}x{PROD_NPHI:02d}x{PROD_NZ:03d}.txt"

        # vecpsi: (Nr+1) x Nphi x (Nz+1), r-z edge
        src = read_ascii_grid(
            prod_dir / f"vecpsi_{psuffix}", "vecpsi", PROD_NR, PROD_NPHI, PROD_NZ
        )[:, 0, :]
        pl = resample_bilinear(src, nr + 1, nz + 1)
        write_field(outdir / f"vecpsi_{suffix}",
                    np.repeat(pl[:, None, :], nphi, axis=1))

        # vecg: Nr x (Nphi+1) x Nz, phi-face. The production file's last phi
        # plane is constant padding; reproduce that structure.
        src = read_ascii_grid(
            prod_dir / f"vecg_{psuffix}", "vecg", PROD_NR, PROD_NPHI, PROD_NZ
        )
        pl = resample_bilinear(src[:, 0, :], nr, nz)
        g = np.empty((nr, nphi + 1, nz))
        for p in range(nphi):
            g[:, p, :] = pl
        g[:, nphi, :] = src[:, PROD_NPHI, :].flat[0]  # the padding value
        write_field(outdir / f"vecg_{suffix}", g)
        print("WARNING: vecpsi/vecg are interpolated, NOT a discrete "
              "Grad-Shafranov solution on this mesh. Do not use them to "
              "certify physics.")

    # Verify what we wrote round-trips through the reader the solver mimics.
    for kind in (["veceta"] + (["vecpsi", "vecg"] if args.with_equilibrium else [])):
        path = outdir / f"{kind}_{suffix}"
        arr = read_ascii_grid(path, kind, nr, nphi, nz)
        assert arr.size == expected_count(kind, nr, nphi, nz)
        print(f"verified {path.name}: shape {arr.shape} matches the C layout")

    print(f"\nrun with:  MHD_Config/input_folder={outdir}  "
          f"Numerical/NR={nr} Numerical/Nphi={nphi} Numerical/NZ={nz}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
