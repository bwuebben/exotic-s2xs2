# decide_t4.g — T4 known-answer test of the pi1 pipeline (see the paper: known-answer calibration).
#
# Closed model: T4 = trivial T2-bundle over T2; Lagrangian tori
# T1 = X x A1, T2 = Y x A2 (Baldridge-Kirk arXiv:math/0703065 Sec.2 config).
# Our membrane derivation gives, after Tietze-eliminating the meridians:
#   gens x,y,a,b
#   rels [x,y], [a,x], [a,y], [a,b],
#        w1^e1 * x^q1   (filling of T1; w1 = correction word for mu1, placement L/R)
#        w2^e2 * y^q2   (filling of T2; w2 = correction word for mu2)
# Ground truths (paper: known-answer calibration; Baldridge-Kirk arXiv:math/0703065 Sec.2):
#   GT2: every double-surgery convention presents (Z^2 x|_A Z) x Z, A != I:
#        H1 = Z^2, non-abelian, fingerprint matches G_true(tr3) or G_true(tr1).
#   GT3: single surgery presents H_3(Z) x Z: H1 = Z^3, hits Heisenberg mod 3.
#   GT4: ANY abelian/trivial double-surgery output = pipeline bug.
# Run: gap -q -A decide_t4.g > t4_out.txt

start := Runtime();
Elapsed := function() return Concatenation("[", String(Int((Runtime()-start)/1000)), "s]"); end;

LI_DEPTH := 5;   # low-index fingerprint depth

# ---------- fingerprint machinery ----------
FP := function(G)
  local fp, i, subs, H, inv;
  fp := [ SortedList(AbelianInvariants(G)) ];
  for i in [2..LI_DEPTH] do
    subs := LowIndexSubgroupsFpGroup(G, i);
    inv := [];
    for H in subs do
      if Index(G, H) = i then
        Add(inv, SortedList(AbelianInvariants(H)));
      fi;
    od;
    Add(fp, SortedList(inv));
  od;
  return fp;
end;

# ---------- ground-truth groups ----------
# G_true(A) = <x,y,b,a | [x,y],[a,x],[a,y],[a,b], bxb^-1 = x^A11 y^A21, byb^-1 = x^A12 y^A22>
GTrue := function(A)
  local F, x, y, a, b, rels;
  F := FreeGroup("x","y","a","b");
  x := F.1; y := F.2; a := F.3; b := F.4;
  rels := [ Comm(x,y), Comm(a,x), Comm(a,y), Comm(a,b),
            b*x*b^-1 * (x^A[1][1]*y^A[2][1])^-1,
            b*y*b^-1 * (x^A[1][2]*y^A[2][2])^-1 ];
  return F / rels;
end;

Atr3 := [[2,1],[1,1]];      # hyperbolic, trace 3 (Sol x S1 type)
Atr1 := [[0,1],[-1,1]];     # elliptic, order 6, trace 1
Apar := [[1,1],[0,1]];      # parabolic: H3(Z) x Z (single surgery)
AparI := [[1,-1],[0,1]];    # parabolic, other sign

Print("computing ground-truth fingerprints ", Elapsed(), "\n");
fpTr3 := FP(GTrue(Atr3));;
fpTr1 := FP(GTrue(Atr1));;
fpPar := FP(GTrue(Apar));;
fpParI := FP(GTrue(AparI));;
Print("  fp(tr3) = fp(tr1)?  ", fpTr3 = fpTr1,
      "   (false required for discrimination)\n");
Print("  fp(par) = fp(parI)? ", fpPar = fpParI, "\n");

# ---------- explicit non-abelian finite targets (no SmallGroup needed) ----------
# affine group (Z/3)^2 x|_M <M mod 3> as permutations of the 9 vectors
VecIdx := function(v) return 3*(v[1] mod 3) + (v[2] mod 3) + 1; end;
AllVecs := [];
for i in [0..2] do for j in [0..2] do Add(AllVecs, [i,j]); od; od;
TransPerm := function(t)
  return PermList(List(AllVecs, v -> VecIdx([v[1]+t[1], v[2]+t[2]])));
end;
LinPerm := function(M)
  return PermList(List(AllVecs,
    v -> VecIdx([M[1][1]*v[1]+M[1][2]*v[2], M[2][1]*v[1]+M[2][2]*v[2]])));
end;
AffGroup := function(M)
  return Group(TransPerm([1,0]), TransPerm([0,1]), LinPerm(M));
end;
# Heisenberg mod 3 = extraspecial group of order 27, exponent 3 (core GAP)
H27 := ExtraspecialGroup(27, 3);
NonAbTargets := [ SymmetricGroup(3), AlternatingGroup(4),
                  AffGroup(Atr3), AffGroup(Atr1), H27 ];
Print("targets built: orders ", List(NonAbTargets, Size),
      "  nonabelian ", List(NonAbTargets, g -> not IsAbelian(g)), "\n");

HitsNonAb := function(G)
  local T, r;
  for T in NonAbTargets do
    r := CALL_WITH_CATCH(GQuotients, [G, T]);
    if r[1] = true and Length(r[2]) > 0 then return true; fi;
  od;
  return false;
end;

# ---------- the swept presentations (our derivation) ----------
F := FreeGroup("x","y","a","b");
x := F.1; y := F.2; a := F.3; b := F.4;
baseRels := [ Comm(x,y), Comm(a,x), Comm(a,y), Comm(a,b) ];

# correction word for mu^e as forced by placement:
#   L: b k b^-1 = mu^e k  ->  mu^e = b k b^-1 k^-1
#   R: b k b^-1 = k mu^e  ->  mu^e = k^-1 b k b^-1
CorrWord := function(k, placement)
  if placement = "L" then return b*k*b^-1*k^-1; else return k^-1*b*k*b^-1; fi;
end;

# Lasso-corrected correction value: b k b^-1 = d mu^e d^-1 * k (L)  =>  mu^e = d^-1 (b k b^-1 k^-1) d
# The membrane/slide derivation (run addendum) gives d1 = Id (T1: direct basing
# crosses nothing) and d2 = b^{+-1} (T2: the pushoff-circle slide past A1
# crosses T1 once; wrap-around basing => one b-conjugation). Sweep all 9
# (d1,d2) pairs to confirm exactly the derived pair passes.
Print("\n=== DOUBLE-SURGERY SWEEP with lasso dimension (9 x 64) ", Elapsed(), " ===\n");
lassoNames := ["1","b","B"];   # Id, b, b^-1
lassos := [ One(F), b, b^-1 ];
for iL1 in [1..3] do for iL2 in [1..3] do
  d1 := lassos[iL1]; d2 := lassos[iL2];
  nPass := 0; nAbelian := 0; nWrongFP := 0; nCases := 0;
  passTypes := [];
  for p1 in ["L","R"] do for p2 in ["L","R"] do
  for e1 in [1,-1] do for e2 in [1,-1] do
  for q1 in [1,-1] do for q2 in [1,-1] do
    nCases := nCases + 1;
    rels := Concatenation(baseRels,
              [ (d1^-1 * CorrWord(y,p1) * d1)^e1 * x^q1,     # filling of T1
                (d2^-1 * CorrWord(x,p2) * d2)^e2 * y^q2 ]);  # filling of T2
    G := F / rels;
    h1 := SortedList(AbelianInvariants(G));
    if h1 <> [0,0] then
      nAbelian := nAbelian + 1;
    elif not HitsNonAb(G) then
      nAbelian := nAbelian + 1;
    else
      fpG := FP(G);
      if fpG = fpTr3 then
        nPass := nPass + 1; Add(passTypes, "tr3");
      elif fpG = fpTr1 then
        nPass := nPass + 1; Add(passTypes, "tr1");
      else
        nWrongFP := nWrongFP + 1;
      fi;
    fi;
  od; od; od; od; od; od;
  Print("LASSO (d1,d2)=(", lassoNames[iL1], ",", lassoNames[iL2],
        "): pass=", nPass, "/", nCases,
        " abelian=", nAbelian, " wrong-fp=", nWrongFP,
        "  types=", Collected(passTypes), "  ", Elapsed(), "\n");
od; od;

# ---------- single surgery (GT3): only T1 drilled; [b,x] is clean ----------
Print("\n=== SINGLE-SURGERY SWEEP (8 conventions) ", Elapsed(), " ===\n");
sPass := 0; sFail := 0;
for p1 in ["L","R"] do for e1 in [1,-1] do for q1 in [1,-1] do
  rels := Concatenation(baseRels,
            [ Comm(b,x),                      # clean now: no T2
              CorrWord(y,p1)^e1 * x^q1 ]);
  G := F / rels;
  h1 := SortedList(AbelianInvariants(G));
  fpG := FP(G);
  ok := (h1 = [0,0,0]) and (fpG = fpPar or fpG = fpParI)
        and Length(GQuotients(G, H27)) > 0;
  if ok then sPass := sPass + 1; else sFail := sFail + 1; fi;
  Print("SINGLE ", p1, " e=", e1, " q=", q1, ": H1=", h1,
        " Heis-quotient=", Length(GQuotients(G, H27)) > 0,
        " fp-par-match=", (fpG = fpPar or fpG = fpParI),
        "  ", Elapsed(), "\n");
od; od; od;
Print("SINGLE SUMMARY: pass=", sPass, " fail=", sFail, "\n");

# ---------- complement probe: no fillings -> H1 = Z^4 ----------
Print("\n=== COMPLEMENT PROBE ===\n");
for p1 in ["L","R"] do for e1 in [1,-1] do
  rels := Concatenation(baseRels,
            [ Comm(CorrWord(y,p1)^e1, x), Comm(CorrWord(y,p1)^e1, a),
              Comm(CorrWord(x,p1)^e1, y), Comm(CorrWord(x,p1)^e1, a) ]);
  G := F / rels;
  Print("PROBE ", p1, " e=", e1, ": H1(complement) = ",
        SortedList(AbelianInvariants(G)), "  (expected [0,0,0,0])\n");
od; od;

# ---------- BK cross-check (corroborative): their words, same manifold ----------
# AUDIT NOTE (2026-07-15, session 16c): GAP's Comm(g,h) = g^-1 h^-1 g h, so
# Comm(b^-1,y^-1) below renders BK's mu1 = [b^-1,y^-1] as a CONJUGATE (cyclic
# permutation) of the paper-literal word, spliced into the composite relator
# without conjugating the pushoff. Machine-checked benign here (identical
# fingerprints in all 4 cells with paper-literal words; the two action
# matrices are conjugate) — but "their words" is exact only up to this
# rendering. Recorded 2026-07-15.
Print("\n=== BALDRIDGE-KIRK CROSS-CHECK (corroborative) ", Elapsed(), " ===\n");
for q1 in [1,-1] do for q2 in [1,-1] do
  rels := [ Comm(x,a), Comm(y,a), Comm(y, b*a*b^-1),
            Comm(Comm(x,y), b), Comm(x, Comm(a,b)), Comm(y, Comm(a,b)),
            Comm(x,y), Comm(a,b),                     # caps
            Comm(b^-1,y^-1) * x^q1,                   # mu1 m1^q1, their words
            Comm(x^-1,b) * y^q2 ];                    # mu2 m2^q2
  G := F / rels;
  h1 := SortedList(AbelianInvariants(G));
  fpG := FP(G);
  verdict := "neither";
  if fpG = fpTr3 then verdict := "tr3"; elif fpG = fpTr1 then verdict := "tr1"; fi;
  Print("BK q=(", q1, ",", q2, "): H1=", h1, " nonabelian=", HitsNonAb(F/rels),
        " fp=", verdict, "\n");
od; od;

Print("\nDONE ", Elapsed(), "\n");
QUIT_GAP(0);
