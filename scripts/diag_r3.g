# diag_r3.g — DIAGNOSTIC part 4 (R3-completed presentation) (2026-07-15). NOT a verdict run.
# Full 9-cell x 32-sign grid with the corrected dirTbBase candidates
# (lasso-conjugated M inserted; see diag_dirTbBase.g header). Enumeration
# first, KB certification of every blowup. Two representative candidates
# (u=d pre, eps=+1) and (u=d' pre, eps=+1); the exact convention is to be
# derived, but diag_kb.g showed the (0,0) verdict is insensitive to it.

LoadPackage("kbmag");

F := FreeGroup("x","y","r","s","A","B","M","N");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;;
A := F.5;; B := F.6;; M := F.7;; N := F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
R0 := comm(x,y)*comm(r,s);;
base := [ R0,
  A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
  B*x*B^-1*y,    B*r*B^-1*r^-1 ];;
delta := r^-1;;
mkAs := function(e) return A*s*A^-1 * (N^e*y)^-1; end;;
mkBy := function(e) return B*y*B^-1 * (M^e*y*x)^-1; end;;
mkBs := function(e) return B*s*B^-1 * (delta*M^e*delta^-1*s)^-1; end;;
dirTaBase := A*r^-1;;
dirTaFib  := (r*x)^-1;;
dirTbFib  := s*r^-1*s^-1;;
R3 := B*(s^-1*r^-1*y*x)*B^-1*(r^-1*s^-1*x)^-1;;

certify := function(G)
  local H, rws, kb, gens;
  H := Image(IsomorphismSimplifiedFpGroup(G));
  rws := KBMAGRewritingSystem(H);
  kb := CALL_WITH_CATCH(KnuthBendix, [rws]);
  if kb[1] <> true then return "kb-error"; fi;
  if not IsConfluent(rws) then return "not-confluent"; fi;
  gens := GeneratorsOfGroup(FreeStructureOfRewritingSystem(rws));
  if Size(rws) = 1 and ForAll(gens, z -> IsOne(ReducedForm(rws, z))) then
    return "TRIVIAL";
  fi;
  return Concatenation("confluent-Size=", String(Size(rws)));
end;;

for cand in [ ["pre u=d eps=1",  r^-1*M*r * B],
              ["pre u=d' eps=1", x*M*x^-1 * B] ] do
  Print("=== CANDIDATE ", cand[1], " ===\n");
  for mn in Cartesian([-1,0,1],[-1,0,1]) do
    nTriv := 0;; nKBTriv := 0;; other := [];;
    for eA in [1,-1] do for eB in [1,-1] do
    for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
      rels := Concatenation(base, [R3],
               [mkAs(e3), mkBy(e4), mkBs(e5),
                M*(dirTaBase*dirTaFib^mn[2])^eA,
                N*(cand[2]*dirTbFib^mn[1])^eB]);
      G := F / rels;
      tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
               RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
      if tab <> fail and Length(tab[1]) = 1 then
        nTriv := nTriv + 1;
      elif tab <> fail then
        Add(other, [ [eA,eB,e3,e4,e5], Concatenation("|G|=", String(Length(tab[1]))) ]);
      else
        v := certify(F/rels);
        if v = "TRIVIAL" then nKBTriv := nKBTriv + 1;
        else Add(other, [ [eA,eB,e3,e4,e5], v ]); fi;
      fi;
    od; od; od; od; od;
    Print("cell(", mn[1], ",", mn[2], "): enum-trivial=", nTriv,
          " kb-trivial=", nKBTriv, " other=", Length(other), "\n");
    for rr in other do
      Print("   signs ", rr[1], " -> ", rr[2], "\n");
    od;
  od;
od;
QUIT_GAP(0);
