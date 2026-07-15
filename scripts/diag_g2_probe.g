# diag_g2_probe.g — DIAGNOSTIC (2026-07-15). Why does kb_certify G2 (diagram 2,
# cells (+-1,+1), honest dirTbBase) return NOT CERTIFIED uniformly (0/128)?
# Distinguish: inconclusive (KB no confluence) vs nontrivial-or-mixed.
# Then: same case with the (true, diagram-independent) completion relation R3
# added, and with a deeper enumeration.
LoadPackage("kbmag");
F := FreeGroup("x","y","r","s","A","B","M","N");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;;
A := F.5;; B := F.6;; M := F.7;; N := F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
base := [comm(x,y)*comm(r,s), A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
         B*x*B^-1*y, B*r*B^-1*r^-1];;
R3 := B*(s^-1*r^-1*y*x)*B^-1*(r^-1*s^-1*x)^-1;;
dTaB := A*r^-1;; dTaF := (r*x)^-1;; dTbF := s*r^-1*s^-1;;

mkrels := function(m, n, e3, e4, e5, e6, eA, eB)
  return [ A*s*A^-1*(N^e3*y)^-1,
           B*y*B^-1*(M^e4*y*x)^-1,
           B*s*B^-1*(N^e6 * r^-1*M^e5*r * s)^-1,
           M*(dTaB*dTaF^n)^eA,
           N*((r^-1*M^(-e5)*r*B)*dTbF^m)^eB ];
end;;

certifyVerbose := function(G)
  local H, rws, kb, gens;
  H := Image(IsomorphismSimplifiedFpGroup(G));
  Print("  simplified: gens=", Length(GeneratorsOfGroup(H)),
        " rels=", Length(RelatorsOfFpGroup(H)),
        " len=", Sum(RelatorsOfFpGroup(H), Length), "\n");
  rws := KBMAGRewritingSystem(H);
  kb := CALL_WITH_CATCH(KnuthBendix, [rws]);
  if kb[1] <> true then Print("  KB: ERROR\n"); return; fi;
  if not IsConfluent(rws) then Print("  KB: NOT CONFLUENT (inconclusive)\n"); return; fi;
  gens := GeneratorsOfGroup(FreeStructureOfRewritingSystem(rws));
  Print("  KB: confluent, Size=", Size(rws),
        ", all gens->1: ", ForAll(gens, z -> IsOne(ReducedForm(rws, z))), "\n");
end;;

for cse in [ [1,1, 1,1,1,1,1,1], [-1,1, 1,1,1,1,1,1], [1,1, -1,1,-1,1,1,-1] ] do
  Print("=== G2 case m=", cse[1], " n=", cse[2], " e=", cse{[3..8]}, " ===\n");
  Print(" without R3:\n");
  certifyVerbose(F / Concatenation(base, mkrels(cse[1],cse[2],cse[3],cse[4],cse[5],cse[6],cse[7],cse[8])));
  Print(" with R3:\n");
  certifyVerbose(F / Concatenation(base, [R3], mkrels(cse[1],cse[2],cse[3],cse[4],cse[5],cse[6],cse[7],cse[8])));
od;
QUIT_GAP(0);
