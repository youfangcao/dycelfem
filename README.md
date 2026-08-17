# DyCelFEM

Multiscale modelling and simulation of cell dynamics: a dynamic finite-element
tissue mechanics model coupled to per-cell stochastic reaction networks read
from SBML, with diffusion between neighbouring cells.

The shipped example is the skin wound-healing model — epidermis, dermis,
hypodermis and a fibrin clot — where PDGF from the clot recruits fibroblasts,
fibroblasts produce KGF/TGF-β1/IL-1 and procollagen, and keratinocytes migrate
across the wound bed to re-epithelialise it.

---

## Quick start

```bash
make deps        # fetch + build libSBML into third_party/  (~4 min, once)
make -j8         # build bin/dycelfem
make check       # 2-step smoke test of the example
```

Then run the example end to end — simulate, analyse, plot, encode:

```bash
cd examples/wound-healing
../../scripts/dycelfem_run.sh run.cfg
```

Or just the simulator:

```bash
cd examples/wound-healing
../../bin/dycelfem run.cfg
```

It starts simulating immediately and quits when it reaches `run.steps`. Press
`s` to pause/resume, `q` to quit, `i`/`o` to zoom.

---

## Requirements

| | macOS | Linux |
|---|---|---|
| compiler | Xcode Command Line Tools (`xcode-select --install`) | g++ or clang++ |
| build | `brew install cmake pkg-config` | `cmake`, `pkg-config` |
| OpenGL/GLUT | built in (system frameworks) | `freeglut3-dev libglu1-mesa-dev` (Debian/Ubuntu)<br>`freeglut-devel mesa-libGLU-devel` (Fedora/RHEL) |
| libSBML | `make deps` | `make deps`, plus `libxml2-dev zlib1g-dev libbz2-dev` |
| analysis | `python3`, `matplotlib` (`pip3 install matplotlib`) | same |
| video | `ffmpeg` (optional) | same |

A run needs a graphical session: frames are captured with `glReadPixels` from a
real GLUT window, so this does not work over a plain SSH connection. On a
headless Linux box use `xvfb-run -s "-screen 0 1024x768x24" ./bin/dycelfem run.cfg`.

`make deps` builds **libSBML 5.20.5**. The code was originally written against
libSBML 3.x, which cannot be configured on a current toolchain; every entry
point it uses is API-identical in 5.x. The build sets `WITH_CPP_NAMESPACE=OFF`,
which matters — libSBML 4.0 introduced the optional `libsbml` C++ namespace and
this code predates it.

If you already have libSBML installed system-wide, skip `make deps`; the
Makefile finds it via `pkg-config` or in `/usr/local`, `/opt/homebrew`, `/usr`,
and prints which one it used. To point at a specific copy:

```bash
make SBML_PREFIX=/path/to/libsbml
```

Objects depend on the detected libSBML configuration, so switching sources
forces a rebuild rather than silently linking objects compiled against
different headers.

**A caveat seen in the wild:** some older libSBML 3.x installs under
`/usr/local` have a *relative* `install_name` (`../libsbml.dylib`). Such a
library links fine but cannot be loaded at run time, and `DYLD_LIBRARY_PATH`
does not help because the name contains a slash. Symptoms: the build succeeds,
then the binary aborts immediately with `Library not loaded: ../libsbml.dylib`.
Either use `make deps`, or repair the install:

```bash
sudo install_name_tool -id /usr/local/lib/libsbml.3.4.1.dylib \
     /usr/local/lib/libsbml.3.4.1.dylib
```

---

## Configuring a run

Everything about a simulation lives in one plain-text file. Generate an
annotated template with:

```bash
./bin/dycelfem --write-default-config my.cfg
```

```ini
input = tissue-with-wound.txt      # initial tissue (cell database)

[sbml]                             # one reaction network per cell type,
common       = mcommon.xml         # plus 'common', which carries the
macrophage   = mmacrophage.xml     # diffusion constants D_*, D0 and DeltaT
fibroblast   = mfibroblast.xml
keratinocyte = mkeratinocyte.xml
ecm          = mecm.xml
bm           = mbm.xml
clot         = mclot.xml

[run]
steps          = 300               # 1 step = DeltaT = 60 min
autostart      = 1                 # start without waiting for the 's' key
exit_when_done = 1                 # quit at the step limit, so scripts continue
mode           = 1                 # 0 initialise, 1 simulate, 2 relax
diffusion      = edgelen           # see "Diffusion kernel" below
seed           = 0                 # 0 = clock; >0 = fixed RNG stream

[output]
dir            = out               # created if missing
data_every     = 1                 # write printout<N>.txt every N steps
image_every    = 1                 # write frames every N steps
species_images = 1                 # one frame per species as well
sections       = full              # 'simple' = CELL_LIST+PAIR_LIST only, ~7x smaller
bmp_reference  =                   # leave empty; a correct BMP header is generated

[view]
width      = 1000
height     = 700
autofit    = 1                     # frame the tissue instead of the u/v window
margin     = 0.03                  # baked into the input file
zoom       = 1.5                   # 1 = whole tissue visible; >1 closes in
center_x   = 0                     # view centre in world coordinates
center_y   = 150                   # (omit to use the tissue centre)
show_index = 0                     # cell index labels off
```

Every key can be overridden from the environment, which wins over the file:
`run.steps` → `DYCELFEM_RUN_STEPS`. Useful for sweeps:

```bash
for s in 1 2 3; do
  DYCELFEM_RUN_SEED=$s DYCELFEM_OUTPUT_DIR=out_rep$s ./bin/dycelfem run.cfg
done
```

Relative paths resolve against the directory holding the config, so a run
directory can live anywhere.

`dycelfem <celldatabase>` still works without a config and uses built-in
defaults.

### Diffusion kernel

The inter-cell diffusion coefficient can be computed three ways. This matters:
the choice changes the effective diffusivity by ~3 orders of magnitude, and
`D0` in `mcommon.xml` was calibrated for the first form.

| `run.diffusion` | coefficient | |
|---|---|---|
| `edgelen` | `D·D0·commedgelen` | **default**; reproduces the published results |
| `edgelen_distsq` | `D·D0·commedgelen/distsq` | |
| `distsq` | `D·D0/distsq` | needs `D0` raised by ~3×10³ to give the same physical diffusivity |

Do not pair `distsq` with the shipped `D0`; that combination produces almost no
diffusion and therefore almost no chemotactic migration.

---

## Analysis

```bash
scripts/dycelfem_analyze.py RUNDIR [RUNDIR2 ...] -o analysis --sbml mcommon.xml
scripts/dycelfem_plot.py analysis
```

`RUNDIR` is the run's `output.dir`. Passing several treats them as replicates,
and the reported standard deviation becomes the spread across runs (how the
published figures were made) rather than across cells.

The wound region is `|X| < 200` and `Y > 0`, matching the region
`labelwound.pl` carves out when building the initial condition. Override with
`--wound-radius` / `--wound-depth`.

| file | contents |
|---|---|
| `WoundCytokinesTempoDyns.dat` | per-species mean ± sd **inside the wound**, per hour |
| `AllCytokinesTempoDyns.dat` | the same over the whole tissue |
| `WoundSpeciesTotals.dat` | total amount of each species in the wound |
| `Granulation.dat` | fibroblast / keratinocyte / ECM number and area in the wound bed |
| `WoundSizeHours.dat` | epithelial gap (re-epithelialisation) per hour |
| `CellCounts.dat` | cell counts by type, tissue-wide and in the wound |
| `SpeciesRadialProfile.dat` | species vs distance from the wound centre vs time |

Column layout is `1:Hour 2:S1_mean 3:S1_sd 4:S2_mean …`, unchanged from the
original Perl/gnuplot workflow, so gnuplot scripts still work — adapted ones
are in `scripts/gnuplot/`. `dycelfem_plot.py` produces PNG/PDF/EPS with
matplotlib instead.

**Wound size** is the distance between the innermost keratinocytes on either
side of the wound centre. It starts at ~2× the wound radius and falls to 0 as
the epidermis closes.

## Video

```bash
scripts/makevideo.sh examples/wound-healing
```

Encodes `Healing.mp4` plus one per species from the `.bmp` frames, at 10 fps
H.264. `LEGACY_MPEG4=1` reproduces the original `-vcodec mpeg4 -qscale 1`
encoder.

---

## Things that look wrong but are not

Checked against the published 2014 run output:

* **No macrophages.** `set_cell_type(6)` is commented out in the source and has
  been since 2014 ("ECM → Macrophage previously. Now gone."). No type-6 cells
  appear in the published run either. `mmacrophage.xml` is still loaded.
* **`WoundSignal` and `MMP` stay 0** for the whole run. No reaction in this
  model set produces them, and they are 0 throughout the published run too.
  Their videos are blank fields.
* **Printed timings are 1000× too large.** Elapsed times are formatted as
  `(t_end - t_begin) * 0.001`, which assumes MSVC's `CLOCKS_PER_SEC == 1000`;
  POSIX uses 1,000,000. `152.38 seconds!` means 0.152 s. Display only — these
  values never feed back into the simulation.

## Reproducibility

`run.seed = 0` (the default) seeds from the clock, as the original code did, so
runs are not reproducible. Setting a positive `seed` fixes the RNG stream —
verified: the stochastic reaction step counts become identical run to run, and
the first step is bit-identical.

It does **not** give full multi-step bit-reproducibility. Some address-dependent
ordering in the mechanics perturbs the geometry at the 1e-16 level, which the
diffusion assembly amplifies into ±1 differences in integer species counts
within a couple of steps. Uninitialised memory was ruled out (the divergence
survives `MallocPreScribble`). Compare ensembles and distributions rather than
individual trajectories.

Note also that `rand()` differs between Apple libc, glibc and MSVC, so the same
seed gives different streams on different platforms.

## Layout

```
src/            simulator (FEM mechanics, cell biology, rendering)
libBioModel/    SBML reaction networks and stochastic simulation
common/         run-configuration reader
examples/       wound-healing example: config, SBML models, initial tissue
scripts/        analysis, plotting, video, pipeline driver
doc/            provenance diffs against the original SVN trunk
```

`src/RK4A.cpp` is an unrelated standalone utility (a Runge–Kutta steady-state
solver for the chemical master equation) with its own `main()`. It is excluded
from the build and needs `anyoption.h`, which is not part of this package.

`libBioModel/debug.h` is a reconstruction: the original was missing from the
archive, and every reference to it is a bare `#include` or a commented-out
diagnostic, so an empty header reproduces a release build.

## Credits

DyCelFEM was developed by Jieling Zhao and Youfang Cao. `libBioModel` is by
Youfang Cao. Please cite the accompanying paper when using this software.

## In memoriam

This work was carried out in the laboratory of **Professor Jie Liang** in the
Department of Bioengineering at the University of Illinois at Chicago. Jie was
the principal investigator on the project and a co-author of the paper this
software accompanies. He died in late 2024.

The model and the ideas behind it took shape under his guidance. This release is
dedicated to his memory.
