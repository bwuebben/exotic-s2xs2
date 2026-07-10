# phase2_worker.g — run ONE of the 8 phase-2 probes, selected by the
# environment variable PHASE2_CASE (1..8). Case order matches the sequential
# sweep in phase2_resume.g exactly; case 1 is the resumed cell(1,1)
# e=(1,1,1,1,1) probe (targets L2_31, M11, L3_3 only). Launch all 8 at once
# with phase2_parallel.sh.
Read("phase2_common.g");
c := Int(GAPInfo.SystemEnvironment.PHASE2_CASE);;
CASES := [
  rec( mn := [ 1,1],  e := [ 1, 1, 1, 1, 1], resumed := true  ),
  rec( mn := [ 1,1],  e := [-1,-1,-1,-1,-1], resumed := false ),
  rec( mn := [ 1,1],  e := [ 1,-1, 1,-1, 1], resumed := false ),
  rec( mn := [ 1,1],  e := [-1, 1,-1, 1,-1], resumed := false ),
  rec( mn := [-1,1],  e := [ 1, 1, 1, 1, 1], resumed := false ),
  rec( mn := [-1,1],  e := [-1,-1,-1,-1,-1], resumed := false ),
  rec( mn := [-1,1],  e := [ 1,-1, 1,-1, 1], resumed := false ),
  rec( mn := [-1,1],  e := [-1, 1,-1, 1,-1], resumed := false ) ];;
cs := CASES[c];;
if cs.resumed then
  probe(Concatenation("cell(", String(cs.mn[1]), ",", String(cs.mn[2]),
        ") e=", String(cs.e), " (RESUMED at L2_31)"),
        mkG(cs.mn[1], cs.mn[2], cs.e[1], cs.e[2], cs.e[3], cs.e[4], cs.e[5]),
        REMAINING);
else
  probe(Concatenation("cell(", String(cs.mn[1]), ",", String(cs.mn[2]),
        ") e=", String(cs.e)),
        mkG(cs.mn[1], cs.mn[2], cs.e[1], cs.e[2], cs.e[3], cs.e[4], cs.e[5]),
        QUOTS);
fi;
Print("CASE ", c, " DONE ", stamp(), "\n");
QUIT_GAP(0);
