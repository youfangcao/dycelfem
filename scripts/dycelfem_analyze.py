#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
"""
dycelfem_analyze.py - extract time courses from a DyCelFEM run.

Replaces the Windows-only Perl scripts (parsetissue-allspecies.pl,
parsetissue-ncells.pl, parsetissue-healingrate.pl, parsetissue-proliferation.pl)
with one pass over the printout files. Standard library only.

    ./dycelfem_analyze.py RUNDIR [RUNDIR2 ...] -o analysis

Each RUNDIR is a directory holding printout1.txt, printout2.txt, ... (that is,
the run's output.dir). Give more than one to treat them as replicates.

WOUND REGION
    Cells with |X| < wound-radius and Y > wound-depth, matching the region that
    scripts/labelwound.pl carves out when it builds the initial condition
    (WOUNDRAD = 200, WOUNDDEP = 0). Override with --wound-radius/--wound-depth.

OUTPUTS (in -o, gnuplot/matplotlib friendly, '#' comment headers)
    WoundCytokinesTempoDyns.dat   per-species mean+/-sd over cells IN THE WOUND
    AllCytokinesTempoDyns.dat     per-species mean+/-sd over the WHOLE TISSUE
    WoundSpeciesTotals.dat        per-species SUM over cells in the wound
    Granulation.dat               fibroblast/keratinocyte/ECM counts and areas
                                  in the wound bed
    WoundSizeHours.dat            re-epithelialisation gap (see below)
    CellCounts.dat                cell counts per type, tissue-wide and in wound
    SpeciesRadialProfile.dat      species vs distance from wound centre, per step

    Column layout of the *CytokinesTempoDyns.dat files is
        1:Hour  2:S1_mean 3:S1_sd  4:S2_mean 5:S2_sd  ...
    which is what the archived plot_speciestimedyns*.gnu scripts expect.

MEAN AND SD
    One run   -> sd is the spread ACROSS CELLS within the region at that step.
    Several   -> the per-run regional means are combined, and sd is the spread
                 ACROSS RUNS, which is how the published figures were made.
    Which one was used is recorded in the header of every file.

WOUND SIZE
    Distance between the innermost keratinocytes (type 3) on either side of the
    wound centre - the width of the epithelial gap. It starts at ~2x the wound
    radius and falls to 0 as the epidermis closes over.
"""

import argparse
import math
import os
import re
import sys

# Species order is fixed by the D[0..8] assignment in biology.cpp. Names are read
# from an SBML model when one is given, since that is where the simulator gets
# them (Biology::getSpeciesName).
DEFAULT_SPECIES = ["WoundSignal", "Clot", "Collagen", "PDGF", "KGF",
                   "IL1", "TGFb1", "Procollagen", "MMP"]

# Cell type codes, from the adhesion_pair_array comments in biology.cpp.
CELL_TYPES = {0: "ECM", 1: "BM", 2: "Fibroblast", 3: "Keratinocyte",
              4: "Clot", 5: "Hypodermis", 6: "Macrophage"}

# CELL_LIST column indices (0-based):
#   0 ID, 1 Dead, 2 S_Dead, 3 Mig, 4 mAng, 5 Type, 6 X, 7 Y, 8 Area,
#   9 N_Species, 10.. species values
C_ID, C_DEAD, C_MIG, C_TYPE, C_X, C_Y, C_AREA, C_NSPEC = 0, 1, 3, 5, 6, 7, 8, 9
C_SPECIES0 = 10


def species_names_from_sbml(path, count_hint=None):
    """Pull <species id="..."> in document order. Avoids an XML dependency;
    the models are small and flat."""
    try:
        with open(path, "r", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return None
    names = re.findall(r'<species\s+id="([^"]+)"', text)
    if not names:
        return None
    if count_hint and len(names) > count_hint:
        names = names[:count_hint]
    return names


def find_steps(rundir):
    """Return {step: path} for printout<N>.txt in rundir."""
    steps = {}
    try:
        entries = os.listdir(rundir)
    except OSError as exc:
        sys.exit("ERROR: cannot read %s: %s" % (rundir, exc))
    for name in entries:
        m = re.match(r"^printout(\d+)\.txt$", name)
        if m:
            steps[int(m.group(1))] = os.path.join(rundir, name)
    if not steps:
        sys.exit("ERROR: no printout*.txt found in %s\n"
                 "       Point this at the run's output directory." % rundir)
    return steps


def read_cell_list(path):
    """Yield the CELL_LIST rows of one printout file as lists of floats.

    Handles the CRLF line endings the Windows-era files carry, and tolerates
    both 'full' and 'simple' output.sections.
    """
    rows = []
    inside = False
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line == "CELL_LIST":
                inside = True
                continue
            if line == "END_CELL_LIST":
                break
            if not inside or line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) <= C_NSPEC:
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue          # a trailing geometry field we do not need
    return rows


def mean_sd(values):
    n = len(values)
    if n == 0:
        return 0.0, 0.0
    m = sum(values) / n
    if n < 2:
        return m, 0.0
    var = sum((v - m) ** 2 for v in values) / (n - 1)
    return m, math.sqrt(max(var, 0.0))


class StepStats(object):
    """Everything one printout file contributes, for one run."""

    __slots__ = ("wound_mean", "wound_sum", "tissue_mean", "counts",
                 "wound_counts", "gran_n", "gran_a", "gap", "radial", "ncells",
                 "wound_mean_sd", "tissue_mean_sd")

    def __init__(self):
        self.wound_mean = []
        self.wound_sum = []
        self.tissue_mean = []
        self.counts = {}
        self.wound_counts = {}
        self.gran_n = {}
        self.gran_a = {}
        self.gap = 0.0
        self.radial = {}
        self.ncells = 0


def analyse_step(path, nspecies, wound_radius, wound_depth, radial_bin):
    rows = read_cell_list(path)
    st = StepStats()

    wound_vals = [[] for _ in range(nspecies)]
    tissue_vals = [[] for _ in range(nspecies)]
    radial = {}

    left_edge = None        # innermost keratinocyte with x < 0
    right_edge = None       # innermost keratinocyte with x > 0

    for r in rows:
        if r[C_DEAD] != 0:                      # skip dead cells, as the model does
            continue
        st.ncells += 1
        ctype = int(r[C_TYPE])
        x, y, area = r[C_X], r[C_Y], r[C_AREA]
        st.counts[ctype] = st.counts.get(ctype, 0) + 1

        avail = min(nspecies, len(r) - C_SPECIES0)
        vals = [r[C_SPECIES0 + k] for k in range(avail)]
        vals += [0.0] * (nspecies - avail)

        for k in range(nspecies):
            tissue_vals[k].append(vals[k])

        in_wound = abs(x) < wound_radius and y > wound_depth
        if in_wound:
            st.wound_counts[ctype] = st.wound_counts.get(ctype, 0) + 1
            for k in range(nspecies):
                wound_vals[k].append(vals[k])
            if ctype in (0, 2, 3):
                st.gran_n[ctype] = st.gran_n.get(ctype, 0) + 1
                st.gran_a[ctype] = st.gran_a.get(ctype, 0.0) + area

        # radial profile: distance from the wound centre (0,0)
        d = int(math.hypot(x, y) / radial_bin)
        bucket = radial.setdefault(d, [[] for _ in range(nspecies)])
        for k in range(nspecies):
            bucket[k].append(vals[k])

        # epithelial gap: innermost keratinocytes either side of the midline
        if ctype == 3 and y > wound_depth:
            if x < 0 and (left_edge is None or x > left_edge):
                left_edge = x
            if x > 0 and (right_edge is None or x < right_edge):
                right_edge = x

    st.wound_mean = [mean_sd(v)[0] for v in wound_vals]
    st.wound_sum = [sum(v) for v in wound_vals]
    st.tissue_mean = [mean_sd(v)[0] for v in tissue_vals]
    st.radial = dict((d, [sum(v) for v in b]) for d, b in radial.items())

    if left_edge is not None and right_edge is not None:
        st.gap = max(0.0, right_edge - left_edge)
    else:
        st.gap = 0.0

    # keep the per-cell spread for the single-run case
    st.wound_mean_sd = [mean_sd(v)[1] for v in wound_vals]
    st.tissue_mean_sd = [mean_sd(v)[1] for v in tissue_vals]
    return st


def write_series(path, header_cols, rows, note):
    with open(path, "w") as fh:
        fh.write("# %s\n" % os.path.basename(path))
        fh.write("# %s\n" % note)
        fh.write("# " + "  ".join(header_cols) + "\n")
        for r in rows:
            fh.write("  ".join(("%d" % v) if isinstance(v, int) else ("%.6g" % v)
                               for v in r) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Extract wound/tissue time courses from DyCelFEM printout files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n"
               "  ./dycelfem_analyze.py myrun/out -o myrun/analysis --sbml myrun/mcommon.xml\n"
               "  ./dycelfem_analyze.py rep1/out rep2/out rep3/out -o combined\n")
    ap.add_argument("rundirs", nargs="+",
                    help="directories containing printout*.txt (the run's output.dir)")
    ap.add_argument("-o", "--outdir", default="analysis", help="where to write .dat files")
    ap.add_argument("--sbml", default=None,
                    help="an SBML model (e.g. mcommon.xml) to take species names from")
    ap.add_argument("--wound-radius", type=float, default=200.0,
                    help="wound half-width in X (default 200, from labelwound.pl)")
    ap.add_argument("--wound-depth", type=float, default=0.0,
                    help="wound floor in Y (default 0, from labelwound.pl)")
    ap.add_argument("--radial-bin", type=float, default=25.0,
                    help="bin width for the radial profile (default 25)")
    ap.add_argument("--hours-per-step", type=float, default=1.0,
                    help="DeltaT in hours; 1 step = 60 min in the shipped models")
    ap.add_argument("--max-steps", type=int, default=0, help="stop after N steps (0 = all)")
    args = ap.parse_args()

    nspecies = len(DEFAULT_SPECIES)
    species = DEFAULT_SPECIES
    if args.sbml:
        got = species_names_from_sbml(args.sbml, nspecies)
        if got:
            species = got
            nspecies = len(got)
            print("Species names from %s: %s" % (args.sbml, ", ".join(species)))
        else:
            print("WARNING: no species found in %s; using built-in names." % args.sbml)

    # per-run, per-step statistics
    runs = []
    for rd in args.rundirs:
        steps = find_steps(rd)
        keys = sorted(steps)
        if args.max_steps:
            keys = keys[:args.max_steps]
        print("%s: %d steps (%d..%d)" % (rd, len(keys), keys[0], keys[-1]))
        per_step = {}
        for i, s in enumerate(keys):
            per_step[s] = analyse_step(steps[s], nspecies, args.wound_radius,
                                       args.wound_depth, args.radial_bin)
            if (i + 1) % 25 == 0 or i + 1 == len(keys):
                sys.stdout.write("\r  parsed %d/%d" % (i + 1, len(keys)))
                sys.stdout.flush()
        print()
        runs.append(per_step)

    common = sorted(set.intersection(*[set(r.keys()) for r in runs]))
    if not common:
        sys.exit("ERROR: the run directories share no common step numbers.")
    nrep = len(runs)
    across = "across runs (%d replicates)" % nrep if nrep > 1 else "across cells (single run)"
    print("Combining %d step(s); sd is %s." % (len(common), across))

    os.makedirs(args.outdir, exist_ok=True)
    H = args.hours_per_step

    def species_table(pick_mean, pick_sd):
        rows = []
        for s in common:
            row = [s * H]
            for k in range(nspecies):
                if nrep > 1:
                    vals = [pick_mean(r[s])[k] for r in runs]
                    m, sd = mean_sd(vals)
                else:
                    m = pick_mean(runs[0][s])[k]
                    sd = pick_sd(runs[0][s])[k]
                row += [m, sd]
            rows.append(row)
        return rows

    cols = ["1:Hour"]
    for i, nm in enumerate(species):
        cols += ["%d:%s_mean" % (2 + 2 * i, nm), "%d:%s_sd" % (3 + 2 * i, nm)]

    note_w = ("mean concentration per cell inside the wound "
              "(|X| < %g, Y > %g); sd %s" % (args.wound_radius, args.wound_depth, across))
    write_series(os.path.join(args.outdir, "WoundCytokinesTempoDyns.dat"), cols,
                 species_table(lambda st: st.wound_mean, lambda st: st.wound_mean_sd), note_w)

    note_t = "mean concentration per cell over the whole tissue; sd %s" % across
    write_series(os.path.join(args.outdir, "AllCytokinesTempoDyns.dat"), cols,
                 species_table(lambda st: st.tissue_mean, lambda st: st.tissue_mean_sd), note_t)

    # totals in the wound
    tcols = ["1:Hour"] + ["%d:%s_total" % (2 + i, nm) for i, nm in enumerate(species)]
    trows = []
    for s in common:
        row = [s * H]
        for k in range(nspecies):
            row.append(sum(r[s].wound_sum[k] for r in runs) / nrep)
        trows.append(row)
    write_series(os.path.join(args.outdir, "WoundSpeciesTotals.dat"), tcols, trows,
                 "total amount of each species inside the wound, averaged over runs")

    # Granulation: counts and areas of fibroblast/keratinocyte/ECM in the wound bed
    gcols = ["1:Hour", "2:FNMean", "3:FNSTD", "4:KNMean", "5:KNSTD", "6:ENMean",
             "7:ENSTD", "8:FAMean", "9:FASTD", "10:KAMean", "11:KASTD",
             "12:EAMean", "13:EASTD"]
    grows = []
    for s in common:
        row = [s * H]
        for source, key in ((lambda st, t: st.gran_n.get(t, 0), "n"),
                            (lambda st, t: st.gran_a.get(t, 0.0), "a")):
            for t in (2, 3, 0):                    # Fibroblast, Keratinocyte, ECM
                vals = [source(r[s], t) for r in runs]
                m, sd = mean_sd(vals)
                row += [m, sd]
        grows.append(row)
    write_series(os.path.join(args.outdir, "Granulation.dat"), gcols, grows,
                 "cell number (N) and total area (A) of Fibroblast/Keratinocyte/ECM "
                 "in the wound bed; column order matches the archived Granulation.dat")

    # Wound size: epithelial gap
    wrows = []
    for s in common:
        vals = [r[s].gap for r in runs]
        m, sd = mean_sd(vals)
        wrows.append([s * H, m, sd])
    write_series(os.path.join(args.outdir, "WoundSizeHours.dat"),
                 ["1:Hour", "2:Mean_Gap", "3:STD_Gap"], wrows,
                 "distance between the innermost keratinocytes across the wound "
                 "(epithelial gap); 0 means the epidermis has closed")

    # Cell counts
    types = sorted(set(t for r in runs for s in common for t in r[s].counts))
    ccols = ["1:Hour"] + \
            ["%d:%s" % (2 + i, CELL_TYPES.get(t, "type%d" % t)) for i, t in enumerate(types)] + \
            ["%d:%s_inWound" % (2 + len(types) + i, CELL_TYPES.get(t, "type%d" % t))
             for i, t in enumerate(types)]
    crows = []
    for s in common:
        row = [s * H]
        row += [sum(r[s].counts.get(t, 0) for r in runs) / nrep for t in types]
        row += [sum(r[s].wound_counts.get(t, 0) for r in runs) / nrep for t in types]
        crows.append(row)
    write_series(os.path.join(args.outdir, "CellCounts.dat"), ccols, crows,
                 "cell counts by type, whole tissue then restricted to the wound")

    # Radial profile (gnuplot pm3d friendly: blank line between steps)
    rp = os.path.join(args.outdir, "SpeciesRadialProfile.dat")
    with open(rp, "w") as fh:
        fh.write("# SpeciesRadialProfile.dat\n")
        fh.write("# total amount per %g-unit shell from the wound centre\n" % args.radial_bin)
        fh.write("# 1:Hour  2:Distance  " +
                 "  ".join("%d:%s" % (3 + i, nm) for i, nm in enumerate(species)) + "\n")
        for s in common:
            bins = set()
            for r in runs:
                bins |= set(r[s].radial.keys())
            for d in sorted(bins):
                tot = [0.0] * nspecies
                for r in runs:
                    b = r[s].radial.get(d)
                    if b:
                        for k in range(nspecies):
                            tot[k] += b[k]
                fh.write("%.6g  %.6g  %s\n" % (s * H, d * args.radial_bin,
                                               "  ".join("%.6g" % (v / nrep) for v in tot)))
            fh.write("\n")

    print("\nWrote to %s/:" % args.outdir)
    for f in sorted(os.listdir(args.outdir)):
        if f.endswith(".dat"):
            print("  %-32s %8d bytes" % (f, os.path.getsize(os.path.join(args.outdir, f))))
    print("\nPlot with:  ./dycelfem_plot.py %s" % args.outdir)


if __name__ == "__main__":
    main()
