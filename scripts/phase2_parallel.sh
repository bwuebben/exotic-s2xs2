#!/bin/bash
# Run the 8 independent phase-2 probes as concurrent GAP processes (GAP is
# single-threaded; the cases share nothing, so this is an 8-way wall-clock
# split). One log per case in phase2_par_logs/. -T is REQUIRED (GQuotients
# raises break-loop errors otherwise).
cd "$(dirname "$0")" || exit 1
mkdir -p phase2_par_logs
for c in 1 2 3 4 5 6 7 8; do
  PHASE2_CASE=$c gap -q -A -T phase2_worker.g \
    > "phase2_par_logs/case$c.log" 2>&1 &
done
wait
echo "ALL 8 CASES DONE"
grep -h -e '===' -e '>>>' -e 'UNRESOLVED' -e 'RESOLVED' phase2_par_logs/case*.log
