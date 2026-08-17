#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
"""
dycelfem_plot.py - plot the .dat files written by dycelfem_analyze.py.

    ./dycelfem_plot.py ANALYSISDIR [-o ANALYSISDIR] [--format png|pdf|eps]

Produces:
    WoundSpecies.png       per-species mean +/- sd inside the wound over time
    TissueSpecies.png      the same over the whole tissue
    WoundClosure.png       epithelial gap vs time
    Granulation.png        fibroblast/keratinocyte/ECM number and area in the bed
    CellCounts.png         cell counts by type
    SpeciesRadial.png      species vs distance from wound centre, over time

Needs matplotlib. If you would rather use gnuplot, the .dat column layout is
unchanged from the archived scripts - see scripts/gnuplot/.
"""

import argparse
import os
import sys

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("ERROR: matplotlib is required for plotting.\n"
             "       pip3 install matplotlib\n"
             "       (the .dat files are plain text; gnuplot scripts are in scripts/gnuplot/)")

# A colour per species, stable across every figure so the panels can be compared.
SPECIES_COLOURS = {
    "WoundSignal": "#8c8c8c", "Clot": "#d62728", "Collagen": "#7f7f7f",
    "PDGF": "#e377c2", "KGF": "#1f77b4", "IL1": "#ff7f0e",
    "TGFb1": "#2ca02c", "Procollagen": "#9467bd", "MMP": "#17becf",
}
TYPE_COLOURS = {
    "ECM": "#bdbdbd", "BM": "#8c564b", "Fibroblast": "#e377c2",
    "Keratinocyte": "#1f77b4", "Clot": "#2ca02c", "Hypodermis": "#404040",
    "Macrophage": "#ff7f0e",
}


def read_dat(path):
    """Return (column_labels, rows). Labels come from the '# 1:Name 2:Name' line."""
    labels, rows = [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                body = line.lstrip("#").strip()
                if body and body[0].isdigit() and ":" in body:
                    labels = [tok.split(":", 1)[1] for tok in body.split()]
                continue
            if not line.strip():
                rows.append(None)          # blank separator (radial profile)
                continue
            rows.append([float(v) for v in line.split()])
    return labels, rows


def col(rows, i):
    return [r[i] for r in rows if r is not None]


def band(ax, x, m, sd, label, colour):
    ax.plot(x, m, lw=2, color=colour, label=label)
    lo = [a - b for a, b in zip(m, sd)]
    hi = [a + b for a, b in zip(m, sd)]
    if any(s > 0 for s in sd):
        ax.fill_between(x, lo, hi, color=colour, alpha=0.18, linewidth=0)


def plot_species(datpath, outpath, title, fmt):
    if not os.path.exists(datpath):
        return False
    labels, rows = read_dat(datpath)
    rows = [r for r in rows if r is not None]
    if not rows:
        return False
    x = col(rows, 0)
    names = [labels[i].replace("_mean", "") for i in range(1, len(labels), 2)] \
        if labels else []

    n = len(names)
    ncol = 3
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.8 * nrow), squeeze=False)

    for k, name in enumerate(names):
        ax = axes[k // ncol][k % ncol]
        m = col(rows, 1 + 2 * k)
        sd = col(rows, 2 + 2 * k)
        band(ax, x, m, sd, name, SPECIES_COLOURS.get(name, "#333333"))
        ax.set_title(name, fontsize=11)
        ax.set_xlabel("Hours")
        ax.set_ylabel("mean per cell (a.u.)")
        ax.margins(x=0)
        ax.grid(alpha=0.25, linewidth=0.5)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")

    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

    # combined single-axes version, closer to the archived gnuplot figure
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for k, name in enumerate(names):
        band(ax, x, col(rows, 1 + 2 * k), col(rows, 2 + 2 * k), name,
             SPECIES_COLOURS.get(name, "#333333"))
    ax.set_xlabel("Hours")
    ax.set_ylabel("mean concentration per cell (a.u.)")
    ax.set_title(title + " - all species")
    ax.legend(ncol=3, fontsize=9)
    ax.margins(x=0)
    ax.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath.replace("." + fmt, "_combined." + fmt), dpi=150)
    plt.close(fig)
    return True


def plot_closure(datpath, outpath):
    if not os.path.exists(datpath):
        return False
    _, rows = read_dat(datpath)
    rows = [r for r in rows if r is not None]
    if not rows:
        return False
    x, m, sd = col(rows, 0), col(rows, 1), col(rows, 2)
    fig, ax = plt.subplots(figsize=(8, 5))
    band(ax, x, m, sd, "epithelial gap", "#d62728")
    ax.set_xlabel("Hours")
    ax.set_ylabel("Epithelial gap (model length units)")
    ax.set_title("Wound closure / re-epithelialisation")
    ax.margins(x=0)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return True


def plot_granulation(datpath, outpath):
    if not os.path.exists(datpath):
        return False
    _, rows = read_dat(datpath)
    rows = [r for r in rows if r is not None]
    if not rows:
        return False
    x = col(rows, 0)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    for name, mi in (("Fibroblast", 1), ("Keratinocyte", 3), ("ECM", 5)):
        band(a1, x, col(rows, mi), col(rows, mi + 1), name, TYPE_COLOURS[name])
    a1.set_title("Cells in the wound bed")
    a1.set_ylabel("cell number")
    for name, mi in (("Fibroblast", 7), ("Keratinocyte", 9), ("ECM", 11)):
        band(a2, x, col(rows, mi), col(rows, mi + 1), name, TYPE_COLOURS[name])
    a2.set_title("Area in the wound bed")
    a2.set_ylabel("total area (model units$^2$)")
    for ax in (a1, a2):
        ax.set_xlabel("Hours")
        ax.legend()
        ax.margins(x=0)
        ax.grid(alpha=0.25, linewidth=0.5)
    fig.suptitle("Granulation tissue formation", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return True


def plot_counts(datpath, outpath):
    if not os.path.exists(datpath):
        return False
    labels, rows = read_dat(datpath)
    rows = [r for r in rows if r is not None]
    if not rows or not labels:
        return False
    x = col(rows, 0)
    whole = [(i, l) for i, l in enumerate(labels) if i > 0 and not l.endswith("_inWound")]
    wound = [(i, l) for i, l in enumerate(labels) if l.endswith("_inWound")]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    for i, l in whole:
        a1.plot(x, col(rows, i), lw=2, label=l, color=TYPE_COLOURS.get(l, None))
    a1.set_title("Whole tissue")
    for i, l in wound:
        nm = l.replace("_inWound", "")
        a2.plot(x, col(rows, i), lw=2, label=nm, color=TYPE_COLOURS.get(nm, None))
    a2.set_title("Inside the wound")
    for ax in (a1, a2):
        ax.set_xlabel("Hours")
        ax.set_ylabel("cell number")
        ax.legend(fontsize=9)
        ax.margins(x=0)
        ax.grid(alpha=0.25, linewidth=0.5)
    fig.suptitle("Cell populations", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return True


def plot_radial(datpath, outpath):
    if not os.path.exists(datpath):
        return False
    labels, rows = read_dat(datpath)
    rows = [r for r in rows if r is not None]
    if not rows:
        return False
    names = labels[2:] if len(labels) > 2 else []
    hours = sorted(set(r[0] for r in rows))
    dists = sorted(set(r[1] for r in rows))
    hi = {h: i for i, h in enumerate(hours)}
    di = {d: i for i, d in enumerate(dists)}

    show = [n for n in ("PDGF", "KGF", "TGFb1", "Collagen", "IL1", "Procollagen")
            if n in names] or names[:6]
    ncol = 3
    nrow = (len(show) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.2 * nrow), squeeze=False)
    for k, name in enumerate(show):
        j = names.index(name) + 2
        grid = [[0.0] * len(hours) for _ in dists]
        for r in rows:
            grid[di[r[1]]][hi[r[0]]] = r[j]
        ax = axes[k // ncol][k % ncol]
        im = ax.imshow(grid, aspect="auto", origin="lower", cmap="viridis",
                       extent=(min(hours), max(hours), min(dists), max(dists)))
        ax.set_title(name, fontsize=11)
        ax.set_xlabel("Hours")
        ax.set_ylabel("Distance from wound centre")
        fig.colorbar(im, ax=ax, fraction=0.046)
    for k in range(len(show), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle("Spatial spread of each species over time", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(description="Plot DyCelFEM analysis .dat files.")
    ap.add_argument("analysisdir")
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--format", default="png", choices=("png", "pdf", "eps", "svg"))
    args = ap.parse_args()

    d = args.analysisdir
    o = args.outdir or d
    os.makedirs(o, exist_ok=True)
    f = args.format
    made = []

    def note(ok, name):
        if ok:
            made.append(name)

    note(plot_species(os.path.join(d, "WoundCytokinesTempoDyns.dat"),
                      os.path.join(o, "WoundSpecies." + f),
                      "Molecular species inside the wound", f), "WoundSpecies." + f)
    note(plot_species(os.path.join(d, "AllCytokinesTempoDyns.dat"),
                      os.path.join(o, "TissueSpecies." + f),
                      "Molecular species across the whole tissue", f), "TissueSpecies." + f)
    note(plot_closure(os.path.join(d, "WoundSizeHours.dat"),
                      os.path.join(o, "WoundClosure." + f)), "WoundClosure." + f)
    note(plot_granulation(os.path.join(d, "Granulation.dat"),
                          os.path.join(o, "Granulation." + f)), "Granulation." + f)
    note(plot_counts(os.path.join(d, "CellCounts.dat"),
                     os.path.join(o, "CellCounts." + f)), "CellCounts." + f)
    note(plot_radial(os.path.join(d, "SpeciesRadialProfile.dat"),
                     os.path.join(o, "SpeciesRadial." + f)), "SpeciesRadial." + f)

    if not made:
        sys.exit("ERROR: no .dat files found in %s - run dycelfem_analyze.py first." % d)
    print("Wrote to %s/:" % o)
    for m in made:
        print("  " + m)
    if "WoundSpecies." + f in made:
        print("  WoundSpecies_combined." + f)
        print("  TissueSpecies_combined." + f)


if __name__ == "__main__":
    main()
