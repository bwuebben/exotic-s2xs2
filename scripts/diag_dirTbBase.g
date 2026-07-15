# diag_dirTbBase.g — DIAGNOSTIC (2026-07-15). NOT a verdict run.
#
# Hypothesis under test (from the T4 known-answer calibration, decide_t4.g): the input
# word dirTbBase = B is missing one M-meridian conjugate. The based-pushoff
# formula used in the paper (the direction words) is a bundle identity; realizing it in the
# complement sweeps the whisker s2 around beta-bar, and that annulus crosses
# T_alpha once: beta-bar crosses the alpha-cut once, and s2 crosses c at P_s'
# (the position table's own "c-crossing precedes e-crossing on the s-edge").
# The algebraic intersection is +-1, so no routing avoids it. The honest
# s2-based pushoff is then B with a lasso-conjugated M^{+-1} inserted; lasso
# candidates are the two c-arc routes (delta = r^-1, delta' = x), insertion
# side pre/post B, sign +-1.
#
# This run: LP cell (m,n)=(0,0), all 32 sign conventions, original word plus
# the 8 corrected candidates. Question: does the corrected word flip the cell
# from trivial to blowup/nontrivial?  (Blowup = "not visibly trivial", never
# "nontrivial" — pre-registered language applies.)

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

quickVerdict := function(rels)
  local G, ab, tab;
  G := F / rels;
  ab := AbelianInvariants(G);
  if Length(ab) > 0 then return Concatenation("H1=", String(ab)); fi;
  tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
           RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
  if tab <> fail then return Concatenation("|G|=", String(Length(tab[1]))); fi;
  return "blowup";
end;;

# candidate list: [name, word]
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
  counts := rec(trivial:=0, blowup:=0, finite:=0, h1:=0);;
  for eA in [1,-1] do for eB in [1,-1] do
  for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
    v := quickVerdict(Concatenation(base,
           [mkAs(e3), mkBy(e4), mkBs(e5),
            M*(dirTaBase)^eA,           # (m,n)=(0,0): fib exponents vanish
            N*(cand[2])^eB]));
    if v = "|G|=1" then counts.trivial := counts.trivial+1;
    elif v = "blowup" then counts.blowup := counts.blowup+1;
    elif StartsWith(v, "H1=") then counts.h1 := counts.h1+1;
    else counts.finite := counts.finite+1;
    fi;
  od; od; od; od; od;
  Print("CAND ", cand[1], ":  trivial=", counts.trivial,
        "  blowup=", counts.blowup, "  finite>1=", counts.finite,
        "  H1<>0=", counts.h1, "  (of 32)\n");
od;
QUIT_GAP(0);
