#!/bin/bash
# Phase-3 certification round, farmed as per-(case,target) jobs, 8 concurrent
# GAP processes. Heavy job types (LI7, M12, big PSL(2,q)) are emitted first so
# they don't straggle at the end. One log per job in phase3_logs/.
cd "$(dirname "$0")" || exit 1
mkdir -p phase3_logs
TARGETS="LI7 M12 L2_53 U3_4 L2_49 L2_47 L2_43 L2_41 L2_32 Sz_8 U4_2 L2_37 L3_4 A8 U3_3"
{
  for t in $TARGETS; do for c in 1 2 3 4 5 6 7 8; do echo "$c $t"; done; done
  echo "3 A5"   # phase-2 GQuotients errored on these two; re-run for clean logs
  echo "4 A5"
} | xargs -P 8 -n 2 sh -c \
  'P3_CASE=$0 P3_TARGET=$1 gap -q -A -T phase3_worker.g > "phase3_logs/case${0}_${1}.log" 2>&1'
echo "ALL PHASE3 JOBS DONE"
grep -h -e '>>>' -e 'ERRORED' phase3_logs/*.log
grep -hc 'no quotient' phase3_logs/*.log | awk '{s+=$1} END {print s " clean no-quotient verdicts"}'
