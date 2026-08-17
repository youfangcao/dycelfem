# libBioModel: what this package ships, and why

`libBioModel` (Youfang Cao) is a general library for stochastic reaction-network
modelling. DyCelFEM uses only a small part of it: SBML model loading, the
stoichiometry matrix, the per-cell state vector, and a Gaussian RNG. The full
library additionally contains a discrete chemical-master-equation state-space
sampler that the wound-healing simulator never calls.

This distribution ships **only the reachable subset**. That keeps the package
small and, more importantly, keeps its provenance simple: everything included
was written at the University of Illinois at Chicago (file dates Nov 2011 -
May 2013), whereas some of the removed sampler code descends from earlier work
(`SparseMatrix.cpp` carried working notes dated November 2007).

## How the subset was determined

Not by reading the code, but by asking the linker. Starting from the simulator's
own objects and an empty libBioModel set, link, collect the undefined symbols,
add only the translation units that define them, and repeat until the link
succeeds. That reaches a fixed point in one round:

    GaussianRNG.cpp  Network.cpp  SamplingState.cpp  State.cpp  Stoichiometry.cpp

The retained headers are the `#include` closure of those five files plus the
headers the simulator includes directly. Four headers are kept whose `.cpp` was
removed (`DestinationSet.h`, `MatrixElement.h`, `SLinkedList.h`,
`SparseMatrix.h`): they are pulled in by the retained sources for declarations
only, and nothing calls the missing definitions - which is precisely what the
successful link proves.

## Contents

Kept (16 files):

| | |
|---|---|
| sources | `GaussianRNG.cpp` `Network.cpp` `SamplingState.cpp` `State.cpp` `Stoichiometry.cpp` |
| headers | `DestinationSet.h` `GaussianRNG.h` `MatrixElement.h` `Network.h` `SLinkedList.h` `SamplingState.h` `SparseMatrix.h` `State.h` `StateType.h` `Stoichiometry.h` `debug.h` |

Removed (23 files): `DestinationSet.cpp`, `HashTable.{cpp,h}`,
`HashTableS.{cpp,h}`, `InitialConds.{cpp,h}`, `MatrixElement.cpp`,
`Queue.{cpp,h}`, `SLinkedList.cpp`, `SLinkedListS.{cpp,h}`,
`SLinkedListT.{cpp,h}`, `Sampling.{cpp,h}`, `SparseMatrix.cpp`,
`Stack.{cpp,h}`, `StateSpace.{cpp,h}`, `winlin.h`.

`debug.h` is a reconstruction: the original was missing from the archived
source, and every surviving reference to it is a bare `#include` or a
commented-out diagnostic, so an empty header reproduces a release build.

## Verification

* Clean `make` and `make check` both pass.
* The pruned binary reproduces the full build **byte for byte**: same
  configuration with `run.seed = 4242`, step 1 `printout1.txt` is identical
  between a build over the complete 39-file libBioModel and this 16-file subset.
* Removing `winlin.h` also removed the last need for the `-include unistd.h`
  compiler workaround (it existed only because `winlin.h` mapped `RMDIR` to
  `rmdir` without including `<unistd.h>`). The build no longer needs it.

## The complete library

Nothing here is lost - this is a packaging decision, not a deletion. The full
39-file `libBioModel` remains in the original working tree alongside `source/`
and `branch_dsq9_simpleoutput/`. If you later want to distribute the sampler as
well, copy the removed files back and the Makefile will pick them up
automatically; only `SLinkedListT.cpp` needs re-excluding, since it is an
out-of-line template implementation that cannot compile standalone.
