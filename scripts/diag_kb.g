# diag_kb.g — DIAGNOSTIC part 2 (2026-07-15). NOT a verdict run.
# KB-certify the blowup conventions produced by the corrected dirTbBase
# candidates at the LP cell (0,0)  (see diag_dirTbBase.g).
# Question: hard-trivial (like the old (+-1,+1) cells) or genuinely not trivial?

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

certify := function(G)
  local H, rws, kb, conf, gens, allone;
  H := Image(IsomorphismSimplifiedFpGroup(G));
  rws := KBMAGRewritingSystem(H);
  kb := CALL_WITH_CATCH(KnuthBendix, [rws]);
  if kb[1] <> true then return "kb-error"; fi;
  conf := IsConfluent(rws);
  if not conf then return "not-confluent"; fi;
  gens := GeneratorsOfGroup(FreeStructureOfRewritingSystem(rws));
  allone := ForAll(gens, z -> IsOne(ReducedForm(rws, z)));
  if Size(rws) = 1 and allone then return "TRIVIAL"; fi;
  return Concatenation("confluent, Size=", String(Size(rws)),
                       ", gens all 1: ", String(allone));
end;;

cands := [ ["ORIGINAL B", B] ];;
for u in [ ["d",r^-1], ["d'",x] ] do
  for eps in [1,-1] do
    Add(cands, [Concatenation("pre  u=",u[1]," eps=",String(eps)),
                u[2]*M^eps*u[2]^-1 * B]);
    Add(cands, [Concatenation("post u=",u[1]," eps=",String(eps)),
                B * u[2]*M^eps*u[2]^-1]);
  od;
od;

for cand in cands do
  nTriv := 0;; nBlow := 0;; results := [];;
  for eA in [1,-1] do for eB in [1,-1] do
  for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
    rels := Concatenation(base,
             [mkAs(e3), mkBy(e4), mkBs(e5),
              M*(dirTaBase)^eA, N*(cand[2])^eB]);
    G := F / rels;
    tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
             RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
    if tab <> fail and Length(tab[1]) = 1 then
      nTriv := nTriv + 1;
    else
      nBlow := nBlow + 1;
      Add(results, [ [eA,eB,e3,e4,e5], certify(F/rels) ]);
    fi;
  od; od; od; od; od;
  Print("CAND ", cand[1], ": enum-trivial=", nTriv, "\n");
  for rr in results do
    Print("   signs ", rr[1], " -> KB: ", rr[2], "\n");
  od;
od;
QUIT_GAP(0);
