#!/bin/bash
# maf_certify.sh — independent-engine verification (session 14) of the 8
# representative (±1,+1)-cell presentations, plus the surface-group control,
# using MAF (Alun Williams' Monoid Automata Factory, v2.2.1 sources,
# https://sourceforge.net/projects/maffsa/) — an implementation of
# Knuth–Bendix / automatic structures fully independent of GAP and kbmag.
#
# Inputs: the GASP/KBMAG-format rewriting systems in maf_runs/, produced by
#   gap -q -A -T maf_export.g
# Prereq: MAF binaries (automata, fsacount, reduce); set $MAF to their dir.
# NOTE (build): MAF must be compiled WITHOUT optimization (-O0 -w
# -funsigned-char -fno-strict-aliasing -fpermissive -std=gnu++98, g++); an
# -Os build with GCC 16 on arm64 miscompiles and segfaults.
#
# Expected output per case: "accepted language contains 1 words" and every
# generator reducing to IdWord; surface control: infinite language.
set -e
cd "$(dirname "$0")/maf_runs"
MAF=${MAF:-$HOME/opt/maf/bin-safe}
for i in 1 2 3 4 5 6 7 8; do
  "$MAF/automata" "case$i" > /dev/null 2>&1
  echo "case $i: $("$MAF/fsacount" "case$i.wa" 2>&1 | tail -1)"
  for g in _g1 _g3 _g5; do
    echo "  reduce($g): $("$MAF/reduce" "case$i" "$g" 2>&1)"
  done
done
"$MAF/automata" surface > /dev/null 2>&1
echo "surface control: $("$MAF/fsacount" surface.wa 2>&1 | tail -1)"
