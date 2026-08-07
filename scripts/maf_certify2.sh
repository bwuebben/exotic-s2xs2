#!/bin/bash
# maf_certify2.sh — session-17 rerun of the independent-engine verification
# on the CORRECTED-word exports (maf_runs2/, written by maf_export2.g with
# the honest dirTbBase = r^-1 M^(-e5) r B and R3): the 8 representative
# (±1,+1)-cell presentations plus the surface-group control, run through
# MAF (Alun Williams' Monoid Automata Factory, v2.2.1) — an implementation
# of Knuth–Bendix / automatic structures fully independent of GAP and kbmag.
#
# Engine provenance note (session 17): this run uses the AUTHOR'S OFFICIAL
# v2.2.1 macpro-x64 binaries (sourceforge.net/projects/maffsa, run under
# Rosetta 2; dylib install-name fixed to @executable_path, ad-hoc resigned),
# validated by first reproducing the committed session-14 record
# maf_out.txt from regenerated maf_runs/ exports LINE-IDENTICALLY.
# The self-built -O0 binaries on the second MacBook (see maf_certify.sh's
# build note) remain a third, independently produced binary.
#
# Expected per case: "accepted language contains 1 words" and every
# generator reducing to IdWord; surface control: infinite language.
set -e
cd "$(dirname "$0")/maf_runs2"
MAF=${MAF:-$HOME/opt/maf/bin-official}
for i in 1 2 3 4 5 6 7 8; do
  "$MAF/automata" "case$i" > /dev/null 2>&1
  echo "case $i: $("$MAF/fsacount" "case$i.wa" 2>&1 | tail -1)"
  for g in _g1 _g3 _g5; do
    echo "  reduce($g): $("$MAF/reduce" "case$i" "$g" 2>&1)"
  done
done
"$MAF/automata" surface > /dev/null 2>&1
echo "surface control: $("$MAF/fsacount" surface.wa 2>&1 | tail -1)"
