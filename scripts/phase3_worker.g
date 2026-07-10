# phase3_worker.g — one (case, target) job of the phase-3 certification round.
# Env: P3_CASE = 1..8 (same case order as phase2_worker.g),
#      P3_TARGET = target name below, or "LI7" (low-index <= 7), or "A5"
#      (re-run of the two jobs GQuotients errored out of in phase 2).
# Phase-3 target ladder = U3(3) (closes the one hole in the phase-2 ladder:
# with it, ALL nonabelian simple groups of order <= 14880 are covered) + every
# remaining nonabelian simple group of order <= 10^5, per the classification:
#   U3_3 6048, A8 20160, L3_4 20160, L2_37 25308, U4_2 25920, Sz_8 29120,
#   L2_32 32736, L2_41 34440, L2_43 39732, L2_47 51888, L2_49 58800,
#   U3_4 62400, L2_53 74412, M12 95040.
# Since H1 = 0 (perfect), a nontrivial finite quotient of order <= X exists
# iff a simple quotient of order <= X exists — so this ladder decides
# "no nontrivial finite quotient of order <= 10^5" for each case.
Read("phase2_common.g");
c := Int(GAPInfo.SystemEnvironment.P3_CASE);;
tname := GAPInfo.SystemEnvironment.P3_TARGET;;
CASES := [
  rec( mn := [ 1,1],  e := [ 1, 1, 1, 1, 1] ),
  rec( mn := [ 1,1],  e := [-1,-1,-1,-1,-1] ),
  rec( mn := [ 1,1],  e := [ 1,-1, 1,-1, 1] ),
  rec( mn := [ 1,1],  e := [-1, 1,-1, 1,-1] ),
  rec( mn := [-1,1],  e := [ 1, 1, 1, 1, 1] ),
  rec( mn := [-1,1],  e := [-1,-1,-1,-1,-1] ),
  rec( mn := [-1,1],  e := [ 1,-1, 1,-1, 1] ),
  rec( mn := [-1,1],  e := [-1, 1,-1, 1,-1] ) ];;
mkTarget := function(name)
  local r;
  if name = "A5"    then return AlternatingGroup(5);
  elif name = "U3_3"  then return PSU(3,3);
  elif name = "A8"    then return AlternatingGroup(8);
  elif name = "L3_4"  then return PSL(3,4);
  elif name = "L2_37" then return PSL(2,37);
  elif name = "U4_2"  then return PSU(4,2);
  elif name = "Sz_8"  then
    r := CALL_WITH_CATCH(SuzukiGroup, [IsPermGroup, 8]);
    if r[1] = true then return r[2]; fi;
    return Image(IsomorphismPermGroup(SuzukiGroup(8)));
  elif name = "L2_32" then return PSL(2,32);
  elif name = "L2_41" then return PSL(2,41);
  elif name = "L2_43" then return PSL(2,43);
  elif name = "L2_47" then return PSL(2,47);
  elif name = "L2_49" then return PSL(2,49);
  elif name = "U3_4"  then return PSU(3,4);
  elif name = "L2_53" then return PSL(2,53);
  elif name = "M12"   then return MathieuGroup(12);
  fi;
  Error("unknown target ", name);
end;;
cs := CASES[c];;
G := mkG(cs.mn[1], cs.mn[2], cs.e[1], cs.e[2], cs.e[3], cs.e[4], cs.e[5]);;
iso := IsomorphismSimplifiedFpGroup(G);; H := Image(iso);;
# GQuotients' pre-processing step ExcludedOrders probes generator orders with
# a NEWTC coset enumeration whose limit:=50000 is hard-coded at the call site
# (lib/grpfp.gi:5352, GAP 4.16.0); the cell-(1,1) alternating-sign
# presentations blow past it and the whole GQuotients call dies (this, not
# break-loop flakiness, also caused the phase-2 A5 "errored, skipped" lines —
# and phase 2's later targets survived only because ExcludedOrders marks an
# order as tested in the mutable StoredExcludedOrders cache BEFORE probing it,
# so one process per case self-immunized after the first casualty).
# Fix: pre-mark every order as already tested. Marking as tested (NOT as
# excluded — that direction would be unsound) makes the probe a no-op;
# exclusion only prunes the search space, so skipping it cannot change the
# set of epimorphisms found.
sxo := StoredExcludedOrders(H);;
for i in [1..Length(sxo)] do UniteSet(sxo[i][2], [1..200]); od;
Print("=== case ", c, " cell(", String(cs.mn[1]), ",", String(cs.mn[2]),
      ") e=", String(cs.e), " target ", tname, " ", stamp(), " ===\n");
if tname = "LI7" then
  li := CALL_WITH_CATCH(LowIndexSubgroupsFpGroup, [H, 7]);
  if li[1] = true then
    Print("  low-index<=7: ", Length(li[2]), " subgroup(s) ", stamp(), "\n");
    if Length(li[2]) > 1 then
      Print("  >>> NONTRIVIAL (proper low-index subgroup) <<< ", stamp(), "\n");
    fi;
  else Print("  low-index<=7 ERRORED ", stamp(), "\n"); fi;
else
  T := mkTarget(tname);;
  q := CALL_WITH_CATCH(GQuotients, [H, T]);
  if q[1] = true and Length(q[2]) > 0 then
    Print("  >>> NONTRIVIAL: onto ", tname, " (", Length(q[2]), " maps) <<< ",
          stamp(), "\n");
  elif q[1] = true then
    Print("  no quotient onto ", tname, " ", stamp(), "\n");
  else
    Print("  ", tname, " ERRORED ", stamp(), "\n");
  fi;
fi;
Print("JOB DONE case ", c, " target ", tname, " ", stamp(), "\n");
QUIT_GAP(0);
