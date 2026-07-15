# vperiph2.g — V-PERIPH power upgrade (the 2026-07-15 audit round): low-index coset actions.
# Same must-hold/control words as vperiph.g, evaluated in the permutation
# representations on cosets of ALL subgroups of index <= 5 of the completed
# complement presentation (richer than epimorphisms onto fixed targets).
# Convention (1,1,1). Run: gap -q -A vperiph2.g > vperiph2_out.txt
F := FreeGroup("x","y","r","s","A","B","M","N");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;;
A := F.5;; B := F.6;; M := F.7;; N := F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
base := [ comm(x,y)*comm(r,s),
  A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
  B*x*B^-1*y,    B*r*B^-1*r^-1,
  B*(s^-1*r^-1*y*x)*B^-1*(r^-1*s^-1*x)^-1 ];;
rels := Concatenation(base,
  [ A*s*A^-1*(N*y)^-1, B*y*B^-1*(M*y*x)^-1, B*s*B^-1*(r^-1*M*r*s)^-1 ]);;
G := F / rels;;
words := [
  [ "w1 [M,TaB]",  comm(M, A*r^-1),                  true  ],
  [ "w2 [M,TaF]",  comm(M, (r*x)^-1),                true  ],
  [ "w3 [TaB,TaF]",comm(A*r^-1, (r*x)^-1),           true  ],
  [ "w4 [N,TbB]",  comm(N, r^-1*M^-1*r*B),           true  ],
  [ "w5 [N,TbF]",  comm(N, s*r^-1*s^-1),             true  ],
  [ "w6 [TbB,TbF]",comm(r^-1*M^-1*r*B, s*r^-1*s^-1), true  ],
  [ "w7 CTRL [N,B]",       comm(N, B),               false ],
  [ "w8 CTRL [M,A]",       comm(M, A),               false ],
  [ "w9 CTRL sign-flip",   comm(N, r^-1*M*r*B),      false ] ];;
wg := w -> MappedWord(w, GeneratorsOfGroup(F), GeneratorsOfGroup(G));;
li := LowIndexSubgroupsFpGroup(G, 5);;
Print("subgroups of index <= 5: ", Length(li), "\n");
mustFail := false;; caught := [];;
for H in li do
  hom := FactorCosetAction(G, H);
  for w in words do
    if not IsOne(Image(hom, wg(w[2]))) then
      if w[3] then
        Print("*** MUST-HOLD FAILS (index ", Index(G,H), "): ", w[1], " ***\n");
        mustFail := true;
      else
        if not w[1] in caught then Add(caught, w[1]);
          Print("control caught (index ", Index(G,H), "): ", w[1], "\n"); fi;
      fi;
    fi;
  od;
od;
Print("=== vperiph2 summary: must-hold failures: ", mustFail,
      "; controls caught: ", caught, " ===\n");
QUIT_GAP(0);
