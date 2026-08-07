# phase2_common.g — shared machinery for the phase-2 scout, factored out of
# phase2_resume.g so the 8 independent probes can run as concurrent GAP
# processes (one per case; see phase2_worker.g / phase2_parallel.sh).
F := FreeGroup("x","y","r","s","A","B","M","N");;
x:=F.1;;y:=F.2;;r:=F.3;;s:=F.4;;A:=F.5;;B:=F.6;;M:=F.7;;N:=F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
R0 := comm(x,y)*comm(r,s);;
base := [R0, A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1, B*x*B^-1*y, B*r*B^-1*r^-1,
         B*(s^-1*r^-1*y*x)*B^-1*(r^-1*s^-1*x)^-1];;
mkG := function(m, n, e3, e4, e5, eA, eB)
  return F / Concatenation(base,
    [ A*s*A^-1*(N^e3*y)^-1, B*y*B^-1*(M^e4*y*x)^-1,
      B*s*B^-1*(r^-1*M^e5*r*s)^-1,
      M*((A*r^-1)*((r*x)^-1)^n)^eA, N*(B*(s*r^-1*s^-1)^m)^eB ]);
end;;
QUOTS := [ [AlternatingGroup(5),"A5"], [PSL(2,7),"L2_7"], [AlternatingGroup(6),"A6"],
           [PSL(2,8),"L2_8"], [PSL(2,11),"L2_11"], [PSL(2,13),"L2_13"],
           [PSL(2,16),"L2_16"], [PSL(2,17),"L2_17"], [PSL(2,19),"L2_19"],
           [PSL(2,23),"L2_23"], [PSL(2,25),"L2_25"], [AlternatingGroup(7),"A7"],
           [PSL(2,27),"L2_27"], [PSL(2,29),"L2_29"], [PSL(2,31),"L2_31"],
           [MathieuGroup(11),"M11"], [PSL(3,3),"L3_3"] ];;
# Case 1 (cell (1,1), e=(1,1,1,1,1)) was already tested through L2(29) in the
# paused 2026-07-09 run; only these targets remain for it.
REMAINING := [ [PSL(2,31),"L2_31"], [MathieuGroup(11),"M11"], [PSL(3,3),"L3_3"] ];;
stamp := function() return Concatenation("[", String(Int(Runtime()/1000)), "s]"); end;;
probe := function(label, G, targets)
  local iso, H, tab, t, q, li;
  Print("=== ", label, " ", stamp(), " ===\n");
  iso := IsomorphismSimplifiedFpGroup(G);; H := Image(iso);;
  Print("  simplified gens=", Length(GeneratorsOfGroup(H)),
        " rels=", Length(RelatorsOfFpGroup(H)),
        " len=", Sum(RelatorsOfFpGroup(H), Length), " ", stamp(), "\n");
  tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(H),
           RelatorsOfFpGroup(H), [] : max := 8000000, silent := true);
  if tab <> fail then Print("  DEEP ENUM RESOLVED: |G| = ", Length(tab[1]), " ", stamp(), "\n"); return; fi;
  Print("  enum(8e6): blowup ", stamp(), "\n");
  li := CALL_WITH_CATCH(LowIndexSubgroupsFpGroup, [H, 6]);
  if li[1] = true then
    Print("  low-index<=6: ", Length(li[2]), " subgroup(s) ", stamp(), "\n");
    if Length(li[2]) > 1 then Print("  >>> NONTRIVIAL (proper low-index subgroup) <<<\n"); return; fi;
  else Print("  low-index errored ", stamp(), "\n"); fi;
  for t in targets do
    Print("  trying ", t[2], " ", stamp(), "\n");
    q := CALL_WITH_CATCH(GQuotients, [H, t[1]]);
    if q[1] = true and Length(q[2]) > 0 then
      Print("  >>> NONTRIVIAL: onto ", t[2], " (", Length(q[2]), " maps) <<< ", stamp(), "\n");
      return;
    elif q[1] = false then Print("    (", t[2], " errored, skipped)\n"); fi;
  od;
  Print("  UNRESOLVED: no certificate found ", stamp(), "\n");
end;;
