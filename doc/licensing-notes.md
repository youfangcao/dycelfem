# Licensing notes

**Approved.** The University of Illinois Chicago Office of Technology Management
authorised public release under the BSD 3-Clause licence (Clare Harper, PhD,
Technology Manager, August 2026).

## Files in this package

| file | state |
|---|---|
| `LICENSE` | BSD 3-Clause, © 2026 The Board of Trustees of the University of Illinois |
| `CITATION.cff` | Complete; page numbers pending the BJ corrected proof |
| SPDX headers | Applied to 45 files |

## What OTM specified, verbatim

> 1. Yes, that open source license is okay.
> 2. The copyright label is: ©The Board of Trustees of the University of
>    Illinois 2026. All rights reserved. Permission to use these materials may
>    be obtained from the University of Illinois.
> 3. I think it would be good to have this technology in our records, if you
>    don't mind submitting a disclosure.

### One deviation from that label — confirmed by OTM

The final sentence of the supplied label — *"Permission to use these materials
may be obtained from the University of Illinois"* — is standard university
boilerplate for **all-rights-reserved** material. It states that permission must
be requested. Under BSD 3-Clause, permission is already granted to everyone by
the licence itself, without asking anyone. Including both would make the package
self-contradictory: a reader could not tell whether they may use the code
without writing to UIC first.

`LICENSE` therefore carries the ownership and year exactly as specified:

    Copyright (c) 2026, The Board of Trustees of the University of Illinois
    All rights reserved.

`All rights reserved.` is retained — it is part of the standard BSD 3-Clause
text, and is not in conflict, because the licence grant immediately follows and
defines which rights are licensed. Only the "permission may be obtained"
sentence is omitted.

This reading matches point 1 of the same email, where the BSD 3-Clause licence
was approved. It was put to OTM explicitly and confirmed:

> **Youfang Cao:** Just to confirm — since BSD 3-Clause itself grants permission
> to use, I've used "Copyright (c) 2026, The Board of Trustees of the University
> of Illinois. All rights reserved." as the notice and omitted the sentence
> about obtaining permission, which would contradict the licence grant.
>
> **Clare Harper, PhD (UIC OTM):** That is fine, thank you for checking.

So the omission is authorised, not assumed. Nothing further is outstanding on
the licence text.

### Year

The notice reads 2026, as specified, rather than the 2011–2015 development
period. That is OTM's call and is normal practice for a first publication year.

## Why BSD 3-Clause

The deciding constraint is libSBML, which is **LGPL 2.1**. A permissive licence
links against LGPL code without friction — the LGPL is designed for exactly
that, and only the library itself stays under LGPL.

* **Apache 2.0 was rejected.** Its patent-termination clause is widely held to
  be incompatible with GPLv2/LGPLv2.1. For research code with negligible patent
  exposure, that ambiguity buys nothing.
* **GPL was rejected.** It would force derivatives open, deterring reuse by
  industry and other labs, which works against citation and adoption.
* **BSD 3-Clause over MIT** because clause 3 bars use of the authors' or the
  University's name to endorse derived products.

If copyleft is preferred after all, **LGPL 2.1-or-later** matches libSBML and is
the natural choice — but it is a much larger commitment.

## Why no LGPL obligation attaches here

libSBML is **not redistributed**. `build-libsbml.sh` downloads it from the
upstream project at build time, and the binary links it dynamically:

    $ otool -L bin/dycelfem | grep sbml
        @rpath/libsbml.5.dylib

So this source package triggers no LGPL distribution obligation at all. If
prebuilt binaries are ever shipped, **keep the dynamic linkage** — that is what
satisfies LGPL-2.1 §6, by letting users relink against their own libSBML. Do
not switch to static linking without revisiting this.

## Ownership — resolved

Both authors were UIC employees during development (roughly 2011–2015), so under
standard university IP policy copyright vests in the university. OTM confirmed
this by supplying the Board of Trustees copyright label, and a single approval
covers both contributions.

1. **UIC OTM approval — obtained.** Clare Harper, PhD, Technology Manager,
   August 2026. Correspondence: `build/otm-disclosure-email.md`.
2. **Jieling Zhao's written agreement — obtained.** This mattered because if UIC
   releases rights back to the authors, which is common for non-commercial
   academic code, copyright is then held jointly by Cao and Zhao and licensing
   requires both. With his agreement on file, either outcome from OTM leads to a
   BSD 3-Clause release. Keep the email; it is the record.

Prof. Jie Liang, the principal investigator, died in late 2024. He is credited
as a co-author of the paper through `preferred-citation`. The code itself was
written by Cao and Zhao per the file headers, so the software `authors:` list is
unchanged; if UIC owns the work, his estate does not enter into the licensing
question at all.

### Provenance cleanup already done

`libBioModel` originally contained files whose lineage predates UIC
(`SparseMatrix.cpp` carried working notes dated November 2007). Those files were
unused by the simulator and have been removed from this distribution; see
`doc/libBioModel-pruning.md`. Everything shipped here is dated Nov 2011 – May
2013. The affiliation line of a former institute was also removed from six
analysis scripts in the original working tree.

## Software authors vs. paper authors

The paper has five authors; the **code** has two. Copyright follows authorship
of the code, so `LICENSE` and the `authors:` block of `CITATION.cff` name Cao
and Zhao only. The other three are credited through `preferred-citation`, which
points at the Biophysical Journal paper. Anyone who contributed code not visible
in the file headers should be added.

## Release checklist

Done:

- [x] Jieling Zhao's written agreement
- [x] UIC OTM approval in writing — Clare Harper, August 2026
- [x] Copyright holder set in `LICENSE`, with the "permission may be obtained"
      sentence omitted and that omission confirmed by OTM
- [x] SPDX headers on 45 files
- [x] `CITATION.cff` — ORCID, repository URL, release date, both DOIs
- [x] `.zenodo.json` so Zenodo uses our metadata rather than guessing
- [x] Published: https://github.com/youfangcao/dycelfem
- [x] Tagged `v1.0`, archived on Zenodo
      — concept DOI [10.5281/zenodo.21985815](https://doi.org/10.5281/zenodo.21985815)
      — v1.0 DOI [10.5281/zenodo.21985816](https://doi.org/10.5281/zenodo.21985816)
- [x] `.svn/` and the publisher PDF kept out of the repository (the repo is
      rooted at the package directory, which contains neither)

Outstanding:

- [ ] Submit the UIC invention disclosure, which OTM asked for:
      https://otm.uic.edu/resources/community-support/submit-your-idea-2/
- [ ] Publish the simulation data record on Zenodo and add its DOI here,
      in `README.md`, and in `CITATION.cff`
- [ ] Add volume/issue/pages to `CITATION.cff` once the BJ proof is paginated
