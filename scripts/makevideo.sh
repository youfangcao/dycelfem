#!/bin/bash
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
# -----------------------------------------------------------------------------
#  Build MP4 videos from a DyCelFEM out/ directory of .bmp frames.
#
#  Usage:  build/makevideo.sh <run-directory> [last-frame]
#     e.g. build/makevideo.sh WT_TGFb1_x10_ProcMax_x0.3_ProcKD_x10_DTGFb1_x0.2_rep1
#
#  Frames are read from <run-directory>/out/ and the .mp4 files are written to
#  <run-directory>/, using the same names as the 2014 run (Healing.mp4,
#  Healing_Clot.mp4, ...).
#
#  This replaces the Windows scripts/makevideo.cmd. That script used
#  "-vcodec mpeg4 -qscale 1 -filter:v setpts=2.5*PTS", i.e. MPEG-4 Part 2 at an
#  effective 10 fps (ffmpeg's image2 default input rate is 25). Here the same
#  10 fps is set directly with -framerate, and the codec is H.264 + yuv420p,
#  which every current player and QuickTime handles and which is far smaller at
#  equal quality. To reproduce the original encoder exactly, set
#  LEGACY_MPEG4=1.
#
#  If a species is flat zero for the whole run its video is a blank field; that
#  is the model, not an encoding failure (see README-build.md).
# -----------------------------------------------------------------------------
set -uo pipefail

RUNDIR="${1:?usage: makevideo.sh <run-directory> [last-frame]}"
RUNDIR="${RUNDIR%/}"
OUTDIR="$RUNDIR/out"

[ -d "$OUTDIR" ] || { echo "ERROR: $OUTDIR does not exist"; exit 1; }
command -v ffmpeg >/dev/null || { echo "ERROR: ffmpeg not found (brew install ffmpeg)"; exit 1; }

PREFIXES=(output_ outputWoundSignal_ outputClot_ outputCollagen_ outputPDGF_
          outputKGF_ outputIL1_ outputTGFb1_ outputProcollagen_ outputMMP_)
NAMES=(Healing Healing_WoundSignal Healing_Clot Healing_Collagen Healing_PDGF
       Healing_KGF Healing_IL1 Healing_TGFb1 Healing_Procollagen Healing_MMP)

# Longest run of frames 1..N present for EVERY prefix, so all videos line up.
if [ $# -ge 2 ]; then
	LAST="$2"
else
	LAST=0
	while :; do
		n=$((LAST + 1))
		for p in "${PREFIXES[@]}"; do
			[ -f "$OUTDIR/${p}${n}.bmp" ] || { n=-1; break; }
		done
		[ "$n" -eq -1 ] && break
		LAST=$n
	done
fi
[ "$LAST" -ge 2 ] || { echo "ERROR: fewer than 2 complete frame sets in $OUTDIR"; exit 1; }
echo "==> encoding frames 1..$LAST from $OUTDIR"

for i in "${!PREFIXES[@]}"; do
	p="${PREFIXES[$i]}"
	out="$RUNDIR/${NAMES[$i]}.mp4"
	if [ -n "${LEGACY_MPEG4:-}" ]; then
		ffmpeg -nostdin -v error -y -start_number 1 -i "$OUTDIR/${p}%d.bmp" \
		       -frames:v "$LAST" -qscale 1 -vcodec mpeg4 \
		       -filter:v "setpts=2.5*PTS" "$out"
	else
		ffmpeg -nostdin -v error -y -framerate 10 -start_number 1 \
		       -i "$OUTDIR/${p}%d.bmp" -frames:v "$LAST" \
		       -c:v libx264 -crf 18 -preset medium -pix_fmt yuv420p \
		       -movflags +faststart "$out"
	fi
	if [ -s "$out" ]; then
		printf "  %-26s %s\n" "$(basename "$out")" "$(du -h "$out" | cut -f1)"
	else
		echo "  FAILED: $out"
	fi
done
echo "==> done"
