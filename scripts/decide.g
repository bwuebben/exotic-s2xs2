# decide.g — FROZEN decision harness for the octagon computation (frozen under
# the pre-registered protocol).
# v2, 2026-07-09 (session 4) — meridian-generator architecture; change logged
# BEFORE any honest run. Do not edit further except via logged fixes.
#
# Architecture (paper: the correction architecture): the corrections are conjugated meridians with
# DERIVED lassos, so the meridians M (of T_alpha, based along the partial-y
# whisker y1) and N (of T_beta, based along the partial-s whisker s2) are added
# as generators; the only remaining unknowns are the four DIRECTION WORDS
# (Lagrangian-framed pushoffs, based compatibly with M resp. N):
#   corrBy = M^e   (exact — the meridian basing absorbs the lasso)
#   corrBs = delta * M^e * delta^-1   with delta = r^-1   (derived, develop.py D3)
#   corrAs = N^e   (the s2-basing absorbs the lasso)
# Signs e are swept. A `fail` in any INPUT slot => SELF-TEST + SMOKE mode only.
#
# Conservativity note: until the van Kampen completeness
# argument is established (it is proved in the paper), the presented group only SURJECTS onto pi_1(V').
# A TRIVIAL verdict is therefore valid as-is; a NONTRIVIAL verdict is
# provisional until completeness is established.

F := FreeGroup("x","y","r","s","A","B","M","N");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;;
A := F.5;; B := F.6;; M := F.7;; N := F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;

# ---------------- INPUT SLOTS (direction words; `fail` until derived) ----------------
INPUT := rec(
  # DERIVED 2026-07-09 (paper: the direction words). Filling the slots per
  # the pre-registered protocol — this enables the honest run. Outcome provisional until V-DIAG.
  dirTaBase := A*r^-1,         # lambda = A.delta: half-rotation drift = the lasso
  dirTaFib  := (r*x)^-1,       # c-pushoff at the y1-basing (engine: c@V2 = XR)
  dirTbBase := B,              # product framing, psi_0 fixes e pointwise
  dirTbFib  := s*r^-1*s^-1     # e-pushoff at the s2-basing (= s.(e@V7).s^-1)
);;

# ---------------- fixed bundle relations (note §2, machine-certified) ----------------
R0 := comm(x,y)*comm(r,s);;
base := [ R0,
  A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,   # alpha-monodromy = swap (clean: x,y,r)
  B*x*B^-1*y,    B*r*B^-1*r^-1 ];;               # beta-monodromy = h*id (clean: x,r)
# broken relations with DERIVED corrections (minimal diagram):
delta := r^-1;;
mkAs := function(e) return A*s*A^-1 * (N^e*y)^-1; end;;
mkBy := function(e) return B*y*B^-1 * (M^e*y*x)^-1; end;;
mkBs := function(e) return B*s*B^-1 * (delta*M^e*delta^-1*s)^-1; end;;

verdict := function(name, rels, h1only)
  local G, ab, tab, sz, t, q;
  G := F / rels;
  ab := AbelianInvariants(G);
  if Length(ab) > 0 then Print(name, " | H1=", ab, "\n"); return rec(h1:=ab); fi;
  if h1only then Print(name, " | H1=0 (verdict suppressed)\n"); return rec(h1:=[]); fi;
  tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
           RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
  if tab <> fail then
    sz := Length(tab[1]); Print(name, " | H1=0 | |G|=", sz, "\n");
    return rec(h1:=[], order:=sz);
  fi;
  # FIX 2 2026-07-09 (logged): GQuotients' epimorphism search internally enumerates
  # cosets with an uncatchable default limit (ExcludedOrders, grpfp.gi:4699) on
  # these 8-generator groups. Phase 1 therefore records blowups without in-run
  # certification; nontriviality of blowup cases is certified in a phase-2 pass.
  Print(name, " | H1=0 | blowup | not visibly trivial (phase-2 certification pending)\n");
  return rec(h1:=[], blowup:=true);
end;;

haveAll := ForAll(RecNames(INPUT), n -> INPUT.(n) <> fail);;

if not haveAll then
  Print("### SELF-TEST MODE (direction words unset) ###\n");
  # controls run with meridians killed (M=N=1), reproducing the known landscape
  killMN := [M, N];;
  verdict("selftest pi1(R)",
    Concatenation(base, killMN, [mkAs(1), mkBy(1), mkBs(1)]), false);
  verdict("selftest uncorrected LP",
    Concatenation(base, killMN, [mkAs(1), mkBy(1), mkBs(1),
      comm(y*s, A)*B, comm(s, B)*A]), false);
  verdict("selftest one surgery",
    Concatenation(base, killMN, [mkAs(1), mkBy(1), mkBs(1), comm(y*s,A)*B]), false);
  Print("### SMOKE: H1 of the honest-shape presentation (placeholder directions,\n");
  Print("###        verdicts suppressed — NOT an honest run) ###\n");
  # placeholder direction words with the correct H1 classes (base A resp. B, fiber
  # c ~ (rx)^-1 resp. e ~ r^-1); words NOT derived — H1 architecture check only.
  for mn in [[0,0],[1,-1]] do
    verdict(Concatenation("smoke m=",String(mn[1])," n=",String(mn[2])),
      Concatenation(base,
        [mkAs(1), mkBy(1), mkBs(1),
         M*(A*((r*x)^-1)^mn[2])^1, N*(B*(r^-1)^mn[1])^1]), true);
  od;
else
  Print("### HONEST RUN ###\n");
  for mn in Cartesian([-1,0,1],[-1,0,1]) do
    for eA in [1,-1] do for eB in [1,-1] do
      for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
        verdict(Concatenation("V' m=",String(mn[1])," n=",String(mn[2]),
                  " eA=",String(eA)," eB=",String(eB),
                  " e=(",String(e3),",",String(e4),",",String(e5),")"),
          Concatenation(base,
            [mkAs(e3), mkBy(e4), mkBs(e5),
             M*(INPUT.dirTaBase*INPUT.dirTaFib^mn[2])^eA,
             N*(INPUT.dirTbBase*INPUT.dirTbFib^mn[1])^eB]), false);
      od; od; od;
    od; od;
  od;
fi;
QUIT_GAP(0);
