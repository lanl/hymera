#!/usr/bin/env python
"""Plot the Grad-Shafranov equilibrium input fields and verify their layouts.

inputs/mhd/ holds three ASCII files read by ReadInitialData (src/mhd/mhd.c:100-107):

    veceta -> dataC    cell-centred material tags   Nr      x Nphi     x Nz
    vecpsi -> datapsi  poloidal flux, r-z edge      (Nr+1)  x Nphi     x (Nz+1)
    vecg   -> datag    R*B_phi, phi-face            Nr      x (Nphi+1) x Nz

Only the element count is self-declared in each file; the stagger is implied by
the C indexing arithmetic at ts_functions.c:20326-20348 and nothing validates
it. This script reads each with the assumed layout and plots it: a wrong stride
produces a visibly scrambled or sheared image, so the picture is the check.

psi and g are read ONLY on the ictype 9/15 path -- the production EFIT initial
condition. There is no generator for them in the repo, so they are
irreplaceable primary data.

Usage:
    plot_equilibrium.py [--nr 100] [--nphi 2] [--nz 200] [-o OUTDIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from mhd_data import expected_count, plasma_mask, read_ascii_grid  # noqa: E402


def find_repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "src" / "mhd").is_dir():
            return parent
    raise SystemExit("cannot locate repo root (no src/mhd above this script)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nr", type=int, default=100)
    ap.add_argument("--nphi", type=int, default=2)
    ap.add_argument("--nz", type=int, default=200)
    # Physical extent from inputs/mhd.input via kinetic.cpp; used for axis labels.
    ap.add_argument("--rmin", type=float, default=3.05)
    ap.add_argument("--rmax", type=float, default=9.95)
    ap.add_argument("--zmin", type=float, default=-5.95)
    ap.add_argument("--zmax", type=float, default=5.95)
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--input-folder", default=None)
    args = ap.parse_args()

    root = find_repo_root()
    infolder = Path(args.input_folder) if args.input_folder else root / "inputs" / "mhd"
    outdir = Path(args.outdir) if args.outdir else Path(__file__).parent / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    nr, nphi, nz = args.nr, args.nphi, args.nz
    suffix = f"grid{nr:03d}x{nphi:02d}x{nz:03d}.txt"

    fields = {}
    for kind in ("veceta", "vecpsi", "vecg"):
        path = infolder / f"{kind}_{suffix}"
        want = expected_count(kind, nr, nphi, nz)
        arr = read_ascii_grid(path, kind, nr, nphi, nz)
        fields[kind] = arr
        print(
            f"{path.name}: {arr.size} values, shape {arr.shape} "
            f"(formula gives {want}) OK"
        )

    # Axisymmetry: every phi-plane should be identical, since the problem is 2-D.
    print("\nphi-plane structure:")
    for kind, arr in fields.items():
        nplanes = arr.shape[1]
        equal = [
            bool(np.array_equal(arr[:, 0, :], arr[:, p, :])) for p in range(nplanes)
        ]
        spans = [
            f"plane{p}:[{arr[:, p, :].min():.4g},{arr[:, p, :].max():.4g}]"
            for p in range(nplanes)
        ]
        print(f"  {kind}: {nplanes} planes, equal-to-plane0={equal}")
        print(f"      {'  '.join(spans)}")

    psi = fields["vecpsi"][:, 0, :]
    g = fields["vecg"][:, 0, :]
    tags = fields["veceta"][:, 0, :]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6), constrained_layout=True)

    def rz_extent(shape):
        # Cell/edge counts differ per stagger; extent is the physical box either way.
        return [args.rmin, args.rmax, args.zmin, args.zmax]

    # Panel 1: poloidal flux psi, with contours (flux surfaces) and the plasma
    # boundary from the tag field overlaid. If the stagger were wrong the
    # contours would not close.
    ax = axes[0]
    im = ax.imshow(
        psi.T, origin="lower", aspect="auto", cmap="viridis",
        extent=rz_extent(psi.shape),
    )
    fig.colorbar(im, ax=ax, label="psi  (poloidal flux)")
    rr = np.linspace(args.rmin, args.rmax, psi.shape[0])
    zz = np.linspace(args.zmin, args.zmax, psi.shape[1])
    ax.contour(rr, zz, psi.T, levels=14, colors="white", linewidths=0.5, alpha=0.7)
    rc = np.linspace(args.rmin, args.rmax, tags.shape[0])
    zc = np.linspace(args.zmin, args.zmax, tags.shape[1])
    ax.contour(
        rc, zc, plasma_mask(tags).T.astype(float), levels=[0.5],
        colors="#c44601", linewidths=1.8,
    )
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_title(f"psi  {psi.shape}  (r-z edge)\norange = plasma/separatrix boundary")

    # Panel 2: g = R*B_phi, the toroidal field function.
    ax = axes[1]
    im = ax.imshow(
        g.T, origin="lower", aspect="auto", cmap="magma",
        extent=rz_extent(g.shape),
    )
    fig.colorbar(im, ax=ax, label="g = R*B_phi")
    ax.contour(
        rc, zc, plasma_mask(tags).T.astype(float), levels=[0.5],
        colors="#00b3ff", linewidths=1.8,
    )
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_title(f"g  {g.shape}  (phi-face)")

    # Panel 3: midplane profiles, the sharpest check that the data is smooth
    # and therefore correctly strided.
    ax = axes[2]
    kz_psi, kz_g = psi.shape[1] // 2, g.shape[1] // 2
    ax.plot(
        np.linspace(args.rmin, args.rmax, psi.shape[0]), psi[:, kz_psi],
        color="#4421af", lw=1.8, label=f"psi at ez={kz_psi}",
    )
    ax.set_xlabel("R (m)")
    ax.set_ylabel("psi", color="#4421af")
    ax.tick_params(axis="y", labelcolor="#4421af")
    ax2 = ax.twinx()
    ax2.plot(
        np.linspace(args.rmin, args.rmax, g.shape[0]), g[:, kz_g],
        color="#c44601", lw=1.8, label=f"g at ez={kz_g}",
    )
    ax2.set_ylabel("g", color="#c44601")
    ax2.tick_params(axis="y", labelcolor="#c44601")
    ax.set_title("midplane profiles\n(smooth => stagger decoded correctly)")
    ax.grid(alpha=0.3)

    out = outdir / "equilibrium.png"
    fig.savefig(out, dpi=140)
    print(f"\nwrote {out}")

    # Report which phi-planes carry real data vs padding. datag has Nphi+1
    # planes for the face stagger; a constant plane is padding, not physics.
    print("\nplane content (constant plane => padding, not physics):")
    for kind, arr in fields.items():
        for p in range(arr.shape[1]):
            plane = arr[:, p, :]
            uniq = np.unique(plane)
            kind_of = "CONSTANT (padding)" if uniq.size == 1 else f"{uniq.size} distinct"
            print(f"  {kind} plane {p}: {kind_of}"
                  + (f" = {uniq[0]:.6g}" if uniq.size == 1 else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
