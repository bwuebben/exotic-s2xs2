# vperiph.g — V-PERIPH: peripheral-structure battery (2026-07-15).
#
# NEW validation class, independent of the triviality verdict (hence immune to
# the E2 degeneracy "wrong words generically trivialize"): on each boundary
# 3-torus of the surgery-tori complement C, pi_1 is abelian, so — with the
# same-whisker basings the derivation uses — the meridian and BOTH pushoff
# classes must PAIRWISE COMMUTE in pi_1(C):
#
#   on dnu(T_alpha):  [M, dirTaBase] = [M, dirTaFib] = [dirTaBase, dirTaFib] = 1
#   on dnu(T_beta):   [N, dirTbBase(e5)] = [N, dirTbFib] = [dirTbBase(e5), dirTbFib] = 1
#
# These are NECESSARY conditions for the six input words to be correct
# boundary classes. By the completeness theorem (paper: completeness section: the completed
# presentation, base + R3 + corrections, presents pi_1(C) itself), a
# commutator with NONTRIVIAL image in ANY quotient of the completed complement
# presentation certifies that the corresponding word is NOT a boundary-torus
# class (or that completeness fails) — either way a discovery. Trivial image
# in every reachable quotient = survival evidence, not proof.
#
# POWER CONTROLS (must FAIL for the battery to mean anything):
#   w7 = [N, B]                      — the PRE-repair dirTbBase (the session-15 bug):
#                                      if the battery catches it, this test would have
#                                      found the bug without the T4 calibration.
#   w8 = [M, A]                      — drift dropped from dirTaBase (wrong word for
#                                      the least-tested input).
#   w9 = [N, r^-1*M^(+e5)*r*B]       — correction sign WRONGLY coupled (tests whether
#                                      the session-16 anti-coupling is detectable).
#
# Method: p-quotients (EpimorphismPGroup, core GAP) for p in {2,3,5}, class up
# to 6, of the completed complement presentation, in all 8 correction-sign
# conventions; plus a small GQuotients battery (wrapped in catch).
# Run: gap -q -A vperiph.g > vperiph_out.txt

F := FreeGroup("x","y","r","s","A","B","M","N");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;;
A := F.5;; B := F.6;; M := F.7;; N := F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
R0 := comm(x,y)*comm(r,s);;
base := [ R0,
  A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
  B*x*B^-1*y,    B*r*B^-1*r^-1 ];;
R3 := B*(s^-1*r^-1*y*x)*B^-1*(r^-1*s^-1*x)^-1;;
delta := r^-1;;
mkAs := function(e) return A*s*A^-1 * (N^e*y)^-1; end;;
mkBy := function(e) return B*y*B^-1 * (M^e*y*x)^-1; end;;
mkBs := function(e) return B*s*B^-1 * (delta*M^e*delta^-1*s)^-1; end;;

dirTaBase := A*r^-1;;
dirTaFib  := (r*x)^-1;;
dirTbFib  := s*r^-1*s^-1;;
dirTbBase := function(e5) return r^-1*M^(-e5)*r * B; end;;

# words to test, as functions of e5 (index, name, must-hold?)
mkWords := function(e5)
  return [
    [ "w1 [M,dirTaBase]",    comm(M, dirTaBase),            true  ],
    [ "w2 [M,dirTaFib]",     comm(M, dirTaFib),             true  ],
    [ "w3 [TaBase,TaFib]",   comm(dirTaBase, dirTaFib),     true  ],
    [ "w4 [N,dirTbBase]",    comm(N, dirTbBase(e5)),        true  ],
    [ "w5 [N,dirTbFib]",     comm(N, dirTbFib),             true  ],
    [ "w6 [TbBase,TbFib]",   comm(dirTbBase(e5), dirTbFib), true  ],
    [ "w7 CTRL [N,B(old)]",  comm(N, B),                    false ],
    [ "w8 CTRL [M,A(nodrift)]", comm(M, A),                 false ],
    [ "w9 CTRL [N,sign-flip]", comm(N, r^-1*M^(e5)*r*B),    false ]
  ];
end;;

wordInG := function(G, w)
  return MappedWord(w, GeneratorsOfGroup(F), GeneratorsOfGroup(G));
end;;

Print("=== V-PERIPH: boundary-torus commutation battery ===\n");
Print("completed complement presentation (base + R3 + corrections), 8 sign conventions\n\n");

anyMustFail := false;;
ctrlCaught := rec(w7:=false, w8:=false, w9:=false);;

for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
  rels := Concatenation(base, [R3], [mkAs(e3), mkBy(e4), mkBs(e5)]);
  G := F / rels;
  ab := AbelianInvariants(G);
  words := mkWords(e5);
  Print("--- convention e=(", e3, ",", e4, ",", e5, "): H1(C)=", ab, " ---\n");
  for p in [2,3,5] do
    for cl in [4,6] do
      epi := CALL_WITH_CATCH(EpimorphismPGroup, [G, p, cl]);
      if epi[1] <> true then
        Print("   p=", p, " cl=", cl, ": p-quotient FAILED to compute\n");
        continue;
      fi;
      epi := epi[2];
      Q := Image(epi);
      for w in words do
        img := Image(epi, wordInG(G, w[2]));
        if not IsOne(img) then
          if w[3] then
            Print("   p=", p, " cl=", cl, " |Q|=", Size(Q),
                  ": *** MUST-HOLD FAILS: ", w[1], " nontrivial ***\n");
            anyMustFail := true;
          else
            Print("   p=", p, " cl=", cl, " |Q|=", Size(Q),
                  ": control caught: ", w[1], " nontrivial (GOOD: power)\n");
            if w[1]{[1..2]} = "w7" then ctrlCaught.w7 := true;
            elif w[1]{[1..2]} = "w8" then ctrlCaught.w8 := true;
            elif w[1]{[1..2]} = "w9" then ctrlCaught.w9 := true; fi;
          fi;
        fi;
      od;
    od;
  od;
od; od; od;

Print("\n=== small nonsolvable-quotient battery (convention (1,1,1)) ===\n");
rels := Concatenation(base, [R3], [mkAs(1), mkBy(1), mkBs(1)]);;
G := F / rels;;
words := mkWords(1);;
for T in [ [SymmetricGroup(3),"S3"], [AlternatingGroup(4),"A4"],
           [SymmetricGroup(4),"S4"], [AlternatingGroup(5),"A5"] ] do
  q := CALL_WITH_CATCH(GQuotients, [G, T[1]]);
  if q[1] <> true then Print("  ", T[2], ": GQuotients errored, skipped\n"); continue; fi;
  Print("  ", T[2], ": ", Length(q[2]), " epimorphism(s)\n");
  for hom in q[2] do
    for w in words do
      img := Image(hom, wordInG(G, w[2]));
      if not IsOne(img) then
        if w[3] then
          Print("    *** MUST-HOLD FAILS in ", T[2], ": ", w[1], " ***\n");
          anyMustFail := true;
        else
          Print("    control caught in ", T[2], ": ", w[1], " (GOOD: power)\n");
          if w[1]{[1..2]} = "w7" then ctrlCaught.w7 := true;
          elif w[1]{[1..2]} = "w8" then ctrlCaught.w8 := true;
          elif w[1]{[1..2]} = "w9" then ctrlCaught.w9 := true; fi;
        fi;
      fi;
    od;
  od;
od;

Print("\n=== V-PERIPH SUMMARY ===\n");
if anyMustFail then
  Print("*** SOME MUST-HOLD COMMUTATOR FAILED — a derived word is NOT a\n");
  Print("*** boundary-torus class (or completeness fails). STOP AND THINK.\n");
else
  Print("All six must-hold commutators trivial in every computed quotient,\n");
  Print("all 8 conventions.\n");
fi;
Print("Power: w7 (pre-repair [N,B]) caught = ", ctrlCaught.w7, "\n");
Print("       w8 (no-drift [M,A])   caught = ", ctrlCaught.w8, "\n");
Print("       w9 (sign-flip)        caught = ", ctrlCaught.w9, "\n");
Print("(A control NOT caught means the battery lacks power against that\n");
Print("error class — report honestly; it does NOT invalidate the must-holds.)\n");
QUIT_GAP(0);
