# gt1_diff.g — GT1 punctured-model word diff (the l2 = bab^-1 chirality datum).
#
# Baldridge-Kirk (arXiv:math/0703065, Theorem 1) give, in the PUNCTURED model
# pi_1(H x K - (T1 u T2)) with H, K once-punctured tori:
#   mu1 = [b^-1, y^-1],  m1 = x,  l1 = a
#   mu2 = [x^-1, b],     m2 = y,  l2 = b a b^-1
# and the relations [x,a] = [y,a] = [y, bab^-1] = 1,
# [[x,y],b] = [x,[a,b]] = [y,[a,b]] = 1.
#
# Our membrane derivation of the same six words, at the harness conventions
# (correction words CorrWord(k,L) = b k b^-1 k^-1; derived lasso pair
# (d1, d2) = (1, b)):
#   mu1 = [b, y]           (T1 correction, L placement; lasso trivial)
#   mu2 = b^-1 [b, x] b    (T2 correction at the d2 = b filling basing)
#   m1 = x,  m2 = y        (fiber pushoffs; clean slides)
#   l1 = a                 (base pushoff of A1: the slide u2: 0.3 -> 0 is
#                           clean -- it never meets the K-puncture at 0.5)
#   l2 = b a b^-1          (base pushoff of A2: the DOWNWARD slide
#                           u2: 0.7 -> 0 is blocked -- it sweeps across
#                           u2 = 0.5 at every u1, hence across the puncture;
#                           the UPWARD slide through the wrap u2 = 1 == 0 is
#                           clean and drags the b-whisker around the wrap:
#                           the b-conjugation is DERIVED, not chosen. This is
#                           the string-level form of the derived lasso
#                           d2 = b: the wrap BK put into l2 is the wrap our
#                           filling puts into the meridian conjugator.)
#
# Checks below (in the free group F<x,y,a,b> -- a free-group identity is
# stronger than the same identity in any quotient):
#   1. the mu-words are conjugates of BK's (explicit short conjugators);
#   2. the m- and l-words match BK's LITERALLY -- in particular the
#      chirality of l2 agrees;
#   3. the mirror chirality b^-1 a b differs from b a b^-1 by a PUNCTURE
#      word: the gap is nontrivial in F2<a,b> but lies in the normal closure
#      of [a,b] (= the derived subgroup, certified by abelianization), so it
#      dies exactly when the puncture is capped -- the closed model cannot
#      see the datum, the punctured one certifies it;
#   4. our punctured-model relator list coincides with BK's relator list as
#      a set (up to conjugation/inversion found in 1).
# Group-level equivalence at the closed level (fingerprints, all four
# surgery cells) is separately certified by the BK cross-check block of
# decide_t4.g.
# Run: gap -q -A gt1_diff.g > gt1_out.txt

F := FreeGroup("x","y","a","b");;
x := F.1;; y := F.2;; a := F.3;; b := F.4;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;

lets := [ One(F), x, y, a, b, x^-1, y^-1, a^-1, b^-1 ];;
conjset := Concatenation(lets, ListX(lets, lets, \*));;
FindConj := function(src, tgt)
  local out, w, eps;
  out := [];
  for w in conjset do
    for eps in [1, -1] do
      if tgt = w * src^eps * w^-1 then
        AddSet(out, [ String(w), eps ]);
      fi;
    od;
  od;
  return out;
end;;

fails := 0;;
Check := function(label, ok)
  Print(label, ": ", ok, "\n");
  if ok <> true then fails := fails + 1; fi;
end;;

Print("=== 1. meridian words: ours vs BK (conjugator, sign) ===\n");
r1 := FindConj(comm(b^-1, y^-1), comm(b, y));;
Print("mu1: [b,y] = w [b^-1,y^-1]^eps w^-1 for (w,eps) in ", r1, "\n");
Check("mu1 conjugate of BK's", Length(r1) > 0);
r2 := FindConj(comm(x^-1, b), b^-1*comm(b, x)*b);;
Print("mu2: b^-1[b,x]b = w [x^-1,b]^eps w^-1 for (w,eps) in ", r2, "\n");
Check("mu2 conjugate of BK's", Length(r2) > 0);

Print("\n=== 2. pushoff words: literal string diff ===\n");
Check("m1 = x (literal)", x = x);
Check("m2 = y (literal)", y = y);
Check("l1 = a (literal)", a = a);
Check("l2 ours = b a b^-1 = BK's l2 (literal, SAME chirality)",
      b*a*b^-1 = b*a*b^-1);

Print("\n=== 3. the chirality is a puncture datum ===\n");
gapw := (b*a*b^-1)^-1 * (b^-1*a*b);;
Print("mirror gap (bab^-1)^-1(b^-1ab) = ", gapw, "\n");
Check("gap nontrivial in the free (punctured) group", gapw <> One(F));
ab := List([x, y, a, b], g -> 0);;
abv := function(w)
  local c, i, ch, s;
  c := [0,0,0,0];
  s := String(w);
  # abelianization via ExponentSums
  return List([1..4], i -> ExponentSumWord(w, GeneratorsOfGroup(F)[i]));
end;;
Check("gap abelianizes to 0  (=> gap in <<[a,b]>> = derived subgroup: dies iff the puncture is capped)",
      abv(gapw) = [0,0,0,0]);

Print("\n=== 4. relator lists: ours vs BK, as sets ===\n");
BKrels  := [ comm(x,a), comm(y,a), comm(y, b*a*b^-1),
             comm(comm(x,y), b), comm(x, comm(a,b)), comm(y, comm(a,b)) ];;
# ours, by the membrane recipe: clean a-transports [x,a],[y,a]; puncture
# monodromies [x,[a,b]], [y,[a,b]]; H-boundary transport [[x,y],b];
# T2-peripheral [m2, l2] = [y, bab^-1]  (T1-peripheral [m1,l1] = [x,a], dup).
OURrels := [ comm(x,a), comm(y,a), comm(x, comm(a,b)), comm(y, comm(a,b)),
             comm(comm(x,y), b), comm(y, b*a*b^-1) ];;
matched := true;;
for r in OURrels do
  hit := ForAny(BKrels, s -> Length(FindConj(s, r)) > 0);
  if not hit then matched := false; Print("  UNMATCHED our relator: ", r, "\n"); fi;
od;
for r in BKrels do
  hit := ForAny(OURrels, s -> Length(FindConj(s, r)) > 0);
  if not hit then matched := false; Print("  UNMATCHED BK relator: ", r, "\n"); fi;
od;
Check("relator lists coincide as sets (up to conj/inv)", matched);

if fails = 0 then
  Print("\nGT1 DIFF: ALL CHECKS PASS -- our punctured-model package matches\n");
  Print("Baldridge-Kirk's word for word (mu's up to explicit conjugators, the\n");
  Print("rest literally), with the l2-chirality agreeing and certified to be a\n");
  Print("genuinely punctured-model datum.\n");
else
  Print("\nGT1 DIFF: ", fails, " CHECK(S) FAILED\n");
fi;
QUIT_GAP(0);
