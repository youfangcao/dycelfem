#!/bin/bash
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
# -----------------------------------------------------------------------------
#  dycelfem_run.sh - simulate, analyse, plot and encode in one command.
#
#      scripts/dycelfem_run.sh myrun/run.cfg
#      scripts/dycelfem_run.sh myrun/run.cfg --skip-sim      # re-analyse only
#      scripts/dycelfem_run.sh myrun/run.cfg --no-video
#
#  Everything is driven by run.cfg; this script only sequences the stages and
#  reads output.dir out of the config so the pieces agree on where things live.
# -----------------------------------------------------------------------------
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

CFG="${1:-}"
[ -n "$CFG" ] || { echo "usage: $(basename "$0") <run.cfg> [--skip-sim] [--no-video] [--no-plots]"; exit 1; }
[ -f "$CFG" ] || { echo "ERROR: no such config: $CFG"; exit 1; }
shift

SKIP_SIM=0; NO_VIDEO=0; NO_PLOTS=0
for a in "$@"; do
	case "$a" in
		--skip-sim)  SKIP_SIM=1 ;;
		--no-video)  NO_VIDEO=1 ;;
		--no-plots)  NO_PLOTS=1 ;;
		*) echo "ERROR: unknown option $a"; exit 1 ;;
	esac
done

CFGDIR="$(cd "$(dirname "$CFG")" && pwd)"
CFGFILE="$(basename "$CFG")"

# Read a key from the config, handling both the [output] dir = out form
# and the flat output.dir = out form.
OUTDIR="$(awk -v want="output.dir" -v def="out" '
	/^[ \t]*[#;]/ { next }
	/^[ \t]*\[/   { s=tolower($0); gsub(/[][ \t]/,"",s); section=s; next }
	/=/ { split($0,kv,"="); k=tolower(kv[1]); gsub(/[ \t]/,"",k);
	      v=kv[2]; sub(/[#;].*$/,"",v); gsub(/^[ \t]+|[ \t]+$/,"",v);
	      full=(section=="" ? k : section "." k); if (full==want) out=v }
	END { print (out=="" ? def : out) }' "$CFG")"
SBML="$(awk -v want="sbml.common" -v def="mcommon.xml" '
	/^[ \t]*[#;]/ { next }
	/^[ \t]*\[/   { s=tolower($0); gsub(/[][ \t]/,"",s); section=s; next }
	/=/ { split($0,kv,"="); k=tolower(kv[1]); gsub(/[ \t]/,"",k);
	      v=kv[2]; sub(/[#;].*$/,"",v); gsub(/^[ \t]+|[ \t]+$/,"",v);
	      full=(section=="" ? k : section "." k); if (full==want) out=v }
	END { print (out=="" ? def : out) }' "$CFG")"

RUNOUT="$CFGDIR/$OUTDIR"
ANALYSIS="$CFGDIR/analysis"

BIN=""
for c in "$ROOT/bin/dycelfem" "$ROOT/bin/dycelfem_source" "$ROOT/dycelfem"; do
	[ -x "$c" ] && { BIN="$c"; break; }
done

echo "=============================================================="
echo " config     : $CFG"
echo " output dir : $RUNOUT"
echo " analysis   : $ANALYSIS"
echo " binary     : ${BIN:-<not found>}"
echo "=============================================================="

if [ "$SKIP_SIM" -eq 0 ]; then
	[ -n "$BIN" ] || { echo "ERROR: no dycelfem binary found; run 'make' first."; exit 1; }
	echo "==> [1/4] simulating"
	( cd "$CFGDIR" && "$BIN" "$CFGFILE" ) || { echo "ERROR: simulation failed"; exit 1; }
else
	echo "==> [1/4] simulation skipped"
fi

echo "==> [2/4] analysing"
SBMLARG=()
[ -f "$CFGDIR/$SBML" ] && SBMLARG=(--sbml "$CFGDIR/$SBML")
python3 "$HERE/dycelfem_analyze.py" "$RUNOUT" -o "$ANALYSIS" "${SBMLARG[@]}" \
	|| { echo "ERROR: analysis failed"; exit 1; }

if [ "$NO_PLOTS" -eq 0 ]; then
	echo "==> [3/4] plotting"
	python3 "$HERE/dycelfem_plot.py" "$ANALYSIS" || echo "WARNING: plotting failed (matplotlib missing?)"
else
	echo "==> [3/4] plots skipped"
fi

if [ "$NO_VIDEO" -eq 0 ]; then
	if command -v ffmpeg >/dev/null; then
		echo "==> [4/4] encoding videos"
		"$HERE/../build/makevideo.sh" "$CFGDIR" 2>/dev/null \
			|| "$HERE/makevideo.sh" "$CFGDIR" \
			|| echo "WARNING: video encoding failed"
	else
		echo "==> [4/4] ffmpeg not found; skipping videos"
	fi
else
	echo "==> [4/4] videos skipped"
fi

echo
echo "Done."
echo "  data   : $RUNOUT"
echo "  tables : $ANALYSIS/*.dat"
echo "  plots  : $ANALYSIS/*.png"
