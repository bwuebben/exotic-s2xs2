#!/bin/bash
# Resume the phase-3 certification farm after an interruption, on any machine.
# Needs only GAP on PATH (INSTALL_GAP.md) + this directory from git.
#
# State model: each (case, target) job is atomic. A job counts as done if it
# is listed in phase3_done.txt (the committed manifest — *.log is gitignored,
# so the manifest is what travels across machines) or a local log in
# phase3_logs/ contains a clean "JOB DONE" without "ERRORED". Everything else
# (never started, killed mid-run, or errored) is (re)run, 8-wide. Idempotent:
# interrupt and re-run as often as needed; the manifest and the committed
# evidence file phase3_out.txt are regenerated at the end of every run.
cd "$(dirname "$0")" || exit 1
mkdir -p phase3_logs
TARGETS="M12 L2_53 U3_4 L2_49 L2_47 L2_43 L2_41 L2_32 Sz_8 U4_2 L2_37 L3_4 A8 U3_3 LI7"
REMAINING=$({
  for t in $TARGETS; do for c in 1 2 3 4 5 6 7 8; do echo "$c $t"; done; done
  echo "3 A5"   # the two jobs phase 2's GQuotients bug skipped (see worker)
  echo "4 A5"
} | while read -r c t; do
  grep -qx "$c $t" phase3_done.txt 2>/dev/null && continue
  f="phase3_logs/case${c}_${t}.log"
  [ -f "$f" ] && grep -q "JOB DONE" "$f" && ! grep -q "ERRORED" "$f" && continue
  echo "$c $t"
done)
if [ -n "$REMAINING" ]; then
  echo "$(echo "$REMAINING" | wc -l | tr -d ' ') jobs to run:"; echo "$REMAINING"
  echo "$REMAINING" | xargs -P 8 -n 2 sh -c \
    'P3_CASE=$0 P3_TARGET=$1 gap -q -A -T phase3_worker.g > "phase3_logs/case${0}_${1}.log" 2>&1'
else
  echo "nothing to do — all jobs already done"
fi
echo "PHASE3 RESUME SWEEP DONE"
# regenerate the committed state files
for f in phase3_logs/case*.log; do
  b=$(basename "$f" .log); b=${b#case}; b=${b%_retry}
  if grep -q "JOB DONE" "$f" && ! grep -q "ERRORED" "$f"; then
    echo "${b%%_*} ${b#*_}"
  fi
done | sort -u > phase3_done.txt
n=$(wc -l < phase3_done.txt | tr -d ' ')
{ if [ "$n" -ge 122 ]; then echo "# Phase-3 verdicts (COMPLETE, $n/122 jobs)"
  else echo "# Phase-3 verdicts so far (PARTIAL, $n/122 jobs — resume with phase3_resume.sh)"; fi
  cat phase3_logs/case*.log; } > phase3_out.txt
echo "$n/122 jobs done. Hits (none expected unless a certificate was found):"
grep -h ">>>" phase3_logs/case*.log || echo "  (no nontrivial hits)"
