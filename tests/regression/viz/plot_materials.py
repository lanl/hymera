#!/usr/bin/env python
"""Plot the material region tag field and verify the two redundant copies agree.

The solver reads the same material map twice, in two formats, through two
languages:

  * inputs/mhd/veceta_grid<NR>x<Nphi>x<NZ>.txt  -> user->dataC, read by
    ReadInitialData from src/mhd/mhd.c and consumed at 139 sites in
    src/mhd/mass_matrix_coefficients.c
  * inputs/AxisSymmetricGeometry.dat            -> indicator, read by a bare
    ifstream loop in src/kinetic/kinetic.cpp

Nothing checks that they match. This script checks, and draws the map with the
plasma/separatrix predicate shaded and the three hardcoded "isolated cell"
indices marked.

Usage:
    plot_materials.py [--nr 100] [--nphi 2] [--nz 200] [-o OUTDIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from mhd_data import (  # noqa: E402
    MATERIAL_TAGS,
    PLASMA_REGION_TAGS,
    plasma_mask,
    read_ascii_grid,
    read_geometry_dat,
)

# Hardcoded in alphaecphi_isolcell (src/mhd/mass_matrix_coefficients.c:3361) at
# lines 3426, 3451 and 3469. betaephi_isolcell (line 1326) is the edge-weighting
# function that calls it; the literals are not in betaephi_isolcell itself.
# Reached on the production IC path via FormIFunction_InitializeEP_halo
# (src/mhd/ts_functions.c:10383). Grid-resolution dependent: at any other NR/NZ
# these indices address different cells, or none.
ISOLCELL_INDICES = [(28, 183), (13, 176), (47, 176)]

# Draw tags in physical order from core to exterior so the colourbar reads
# like a radial cut.
TAG_ORDER = [1, 2, 0, -1, -2]
TAG_COLORS = {
    1: "#c44601",   # plasma
    2: "#f57600",   # separatrix-wall
    0: "#8babf1",   # blanket wall
    -1: "#0073e6",  # vacuum vessel
    -2: "#e6e6e6",  # exterior
}


def find_repo_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "src" / "mhd").is_dir():
            return parent
    raise SystemExit("cannot locate repo root (no src/mhd above this script)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nr", type=int, default=100)
    ap.add_argument("--nphi", type=int, default=2)
    ap.add_argument("--nz", type=int, default=200)
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--input-folder", default=None)
    args = ap.parse_args()

    root = find_repo_root()
    infolder = Path(args.input_folder) if args.input_folder else root / "inputs" / "mhd"
    outdir = Path(args.outdir) if args.outdir else Path(__file__).parent / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    eta_path = infolder / f"veceta_grid{args.nr:03d}x{args.nphi:02d}x{args.nz:03d}.txt"
    geom_path = root / "inputs" / "AxisSymmetricGeometry.dat"

    dataC = read_ascii_grid(eta_path, "veceta", args.nr, args.nphi, args.nz)
    print(f"read {eta_path.name}: shape {dataC.shape}")

    # Check every phi-plane is identical; the problem is axisymmetric, so the
    # 3-D array should carry no phi dependence.
    phi_planes_equal = all(
        np.array_equal(dataC[:, 0, :], dataC[:, p, :]) for p in range(args.nphi)
    )
    print(f"all phi-planes identical: {phi_planes_equal}")

    tags = dataC[:, 0, :]  # (er, ez)

    # Cross-check against the C++ side's copy.
    geom_ok = None
    if geom_path.is_file():
        geom = read_geometry_dat(geom_path)
        print(f"read {geom_path.name}: shape {geom.shape}")
        if geom.shape == tags.shape:
            mismatches = int(np.count_nonzero(geom != tags))
            geom_ok = mismatches == 0
            print(f"veceta vs AxisSymmetricGeometry.dat mismatches: {mismatches}")
        else:
            print(
                f"shape mismatch: veceta gives {tags.shape}, geometry file "
                f"gives {geom.shape} -- cannot compare"
            )
    else:
        print(f"note: {geom_path} not found, skipping cross-check")

    present = sorted(set(np.unique(tags).astype(int)), key=lambda t: TAG_ORDER.index(t)
                     if t in TAG_ORDER else 99)
    counts = {int(t): int(np.count_nonzero(tags == t)) for t in present}
    print("tag counts:", {MATERIAL_TAGS.get(t, (str(t),))[0]: c for t, c in counts.items()})

    mask = plasma_mask(tags)
    print(
        f"plasma predicate |dataC-1.5|<0.7 selects {int(mask.sum())} cells "
        f"= tags {PLASMA_REGION_TAGS}"
    )

    # --- figure -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)

    # Panel 1: material map, one discrete colour per tag.
    cmap = matplotlib.colors.ListedColormap([TAG_COLORS[t] for t in TAG_ORDER])
    lookup = {t: i for i, t in enumerate(TAG_ORDER)}
    indexed = np.vectorize(lambda v: lookup.get(int(v), len(TAG_ORDER) - 1))(tags)

    ax = axes[0]
    ax.imshow(
        indexed.T, origin="lower", aspect="auto", cmap=cmap,
        vmin=-0.5, vmax=len(TAG_ORDER) - 0.5, interpolation="nearest",
    )
    ax.set_xlabel("er (radial index)")
    ax.set_ylabel("ez (vertical index)")
    ax.set_title(f"Material regions  {args.nr}x{args.nz}")
    ax.legend(
        handles=[
            mpatches.Patch(
                color=TAG_COLORS[t],
                label=f"{t:+d}  {MATERIAL_TAGS[t][0]} ({MATERIAL_TAGS[t][1]})",
            )
            for t in TAG_ORDER if t in counts
        ],
        loc="upper right", fontsize=8, framealpha=0.95,
    )

    # Panel 2: the predicate that 112 code sites test, plus the isolcell marks.
    ax = axes[1]
    ax.imshow(
        mask.T, origin="lower", aspect="auto", cmap="Greys_r",
        interpolation="nearest", vmin=0, vmax=1,
    )
    ax.set_xlabel("er (radial index)")
    ax.set_ylabel("ez (vertical index)")
    ax.set_title("|dataC - 1.5| < 0.7   (white = plasma or separatrix-wall)")

    for er, ez in ISOLCELL_INDICES:
        inside = 0 <= er < tags.shape[0] and 0 <= ez < tags.shape[1]
        tag = int(tags[er, ez]) if inside else None
        ax.plot(er, ez, "o", ms=11, mfc="none", mec="#c44601", mew=2.0)
        ax.annotate(
            f"({er},{ez}) tag={tag}" if inside else f"({er},{ez}) OUT OF RANGE",
            (er, ez), textcoords="offset points", xytext=(9, 5),
            color="#c44601", fontsize=8, fontweight="bold",
        )
    ax.plot([], [], "o", ms=9, mfc="none", mec="#c44601", mew=2.0,
            label="alphaecphi_isolcell hardcoded cells")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.95)

    out = outdir / "materials.png"
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")

    # --- what are those three cells? -------------------------------------
    print("\nisolcell neighbourhood (4-connected, plane ephi=0):")
    for er, ez in ISOLCELL_INDICES:
        if not (0 < er < tags.shape[0] - 1 and 0 < ez < tags.shape[1] - 1):
            print(f"  (er={er}, ez={ez}): out of range for this grid")
            continue
        nb = {
            "r-": int(tags[er - 1, ez]), "r+": int(tags[er + 1, ez]),
            "z-": int(tags[er, ez - 1]), "z+": int(tags[er, ez + 1]),
        }
        print(f"  (er={er}, ez={ez}) tag={int(tags[er, ez])} neighbours={nb}")

    # Are they isolated in the tag field? Label connected components of the
    # wall region and report which component each lands in.
    wall = tags == 0
    labels = np.zeros(wall.shape, dtype=np.int32)
    ncomp = 0
    for start in zip(*np.nonzero(wall)):
        if labels[start]:
            continue
        ncomp += 1
        stack = [start]
        labels[start] = ncomp
        while stack:
            a, b = stack.pop()
            for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                x, y = a + da, b + db
                if (
                    0 <= x < wall.shape[0] and 0 <= y < wall.shape[1]
                    and wall[x, y] and not labels[x, y]
                ):
                    labels[x, y] = ncomp
                    stack.append((x, y))
    sizes = {c: int(np.count_nonzero(labels == c)) for c in range(1, ncomp + 1)}
    print(f"\nwall region (tag==0): {ncomp} connected component(s), sizes {sizes}")
    for er, ez in ISOLCELL_INDICES:
        if 0 <= er < tags.shape[0] and 0 <= ez < tags.shape[1]:
            c = int(labels[er, ez])
            print(
                f"  (er={er}, ez={ez}) is in component {c} of size "
                f"{sizes.get(c, 0)} -> {'ISOLATED' if sizes.get(c) == 1 else 'NOT isolated'}"
            )

    if geom_ok is False:
        print("\nFAIL: the two copies of the material map disagree")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
