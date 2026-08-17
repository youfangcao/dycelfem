#!/bin/bash
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
# Two-step run of the shipped example: proves the binary, the SBML models, the
# reaction/diffusion path and the image writer all work. Takes about a minute.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EX="$ROOT/examples/wound-healing"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

sed -e 's/^steps *=.*/steps = 2/' \
    -e 's/^dir *=.*/dir = out/' "$EX/run.cfg" > "$TMP/run.cfg"
cp "$EX"/*.xml "$TMP/"
ln -s "$EX/tissue-with-wound.txt" "$TMP/tissue-with-wound.txt"

echo "==> smoke test in $TMP"
( cd "$TMP" && "$ROOT/bin/dycelfem" run.cfg ) > "$TMP/log" 2>&1
rc=$?

fail() { echo "SMOKE TEST FAILED: $1"; echo "--- last 25 log lines ---"; tail -25 "$TMP/log"; exit 1; }

[ $rc -eq 0 ] || fail "binary exited with status $rc"
grep -q "Updating Cell States by Reactions" "$TMP/log" || fail "reactions never ran (check run.mode)"
grep -q "All Diffusions Done"               "$TMP/log" || fail "diffusion never ran"
[ -f "$TMP/out/printout2.txt" ]                        || fail "no printout2.txt written"
[ -f "$TMP/out/output_2.bmp" ]                         || fail "no image written"

cells=$(awk '/^CELL_LIST$/{f=1;next}/^END_CELL_LIST$/{f=0}f' "$TMP/out/printout2.txt" | wc -l)
[ "$cells" -gt 1000 ] || fail "only $cells cells in the output"

echo "PASS: 2 steps, $cells cells, reactions + diffusion + images all produced."
