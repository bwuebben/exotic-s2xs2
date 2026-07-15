# vsens.g — word-sensitivity map at the LP cell (2026-07-15).
#
# Question quantified: WHICH single-word errors does the decision layer
# (H1 + enumeration + KB triviality) even notice at (m,n)=(0,0), convention
# e=(1,1,1,1,1)? The companion's E2 says wrong words make triviality generic,
# so "verdict trivial" cannot certify words; this script maps the blind spot
# for OUR words' immediate neighborhood. Classification per perturbation:
#   LOUD-H1      H1 <> 0                (caught by homology bookkeeping)
#   LOUD-FIN     finite nontrivial      (caught by enumeration)
#   SILENT       trivial (enum or KB)   (verdict layer blind — must be guarded
#                                        by V-PERIPH / T4 / framing instead)
#   MURKY        blowup, KB inconclusive
# Legitimate-equivalent perturbations (arc swaps) are included and EXPECTED
# silent — they present the same group.
# Run: gap -q vsens.g > vsens_out.txt   (kbmag needed for KB fallback)

LoadPackage("kbmag");
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

# honest ingredients at e=(1,1,1,1,1)
hAs := A*s*A^-1*(N*y)^-1;;
hBy := B*y*B^-1*(M*y*x)^-1;;
hBs := B*s*B^-1*(delta*M*delta^-1*s)^-1;;
hTaB := A*r^-1;;  hTaF := (r*x)^-1;;
hTbB := r^-1*M^-1*r*B;;  hTbF := s*r^-1*s^-1;;

mkRels := function(rAs, rBy, rBs, taB, taF, tbB, tbF)
  return Concatenation(base, [R3, rAs, rBy, rBs,
    M*(taB*taF^0)^1, N*(tbB*tbF^0)^1]);
end;;

certify := function(G)
  local H, rws, kb, gens;
  H := Image(IsomorphismSimplifiedFpGroup(G));
  rws := KBMAGRewritingSystem(H);
  kb := CALL_WITH_CATCH(KnuthBendix, [rws]);
  if kb[1] <> true or not IsConfluent(rws) then return "KB-inconclusive"; fi;
  gens := GeneratorsOfGroup(FreeStructureOfRewritingSystem(rws));
  if Size(rws) = 1 and ForAll(gens, z -> IsOne(ReducedForm(rws, z))) then
    return "KB-trivial";
  fi;
  return "KB-other";
end;;

probe := function(name, rels)
  local G, ab, tab, v;
  G := F / rels;
  ab := AbelianInvariants(G);
  if Length(ab) > 0 then
    Print(name, ": LOUD-H1  H1=", ab, "\n"); return;
  fi;
  tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
           RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
  if tab <> fail then
    if Length(tab[1]) = 1 then Print(name, ": SILENT  (enum-trivial)\n");
    else Print(name, ": LOUD-FIN  |G|=", Length(tab[1]), "\n"); fi;
    return;
  fi;
  v := certify(G);
  if v = "KB-trivial" then Print(name, ": SILENT  (blowup, KB-trivial)\n");
  elif v = "KB-other" then Print(name, ": LOUD-KB  (", v, ")\n");
  else Print(name, ": MURKY  (blowup, KB-inconclusive)\n"); fi;
end;;

Print("=== V-SENS: single-word perturbations at (0,0), e=(1,1,1,1,1) ===\n");
probe("baseline (honest words)          ",
  mkRels(hAs,hBy,hBs, hTaB,hTaF, hTbB,hTbF));

Print("\n-- dirTaBase (the drift word) --\n");
probe("TaB := A       (drift dropped)   ", mkRels(hAs,hBy,hBs, A,        hTaF, hTbB,hTbF));
probe("TaB := A*r     (drift inverted)  ", mkRels(hAs,hBy,hBs, A*r,      hTaF, hTbB,hTbF));
probe("TaB := A*x     (other arc; legit)", mkRels(hAs,hBy,hBs, A*x,      hTaF, hTbB,hTbF));
probe("TaB := r^-1*A  (order swapped)   ", mkRels(hAs,hBy,hBs, r^-1*A,   hTaF, hTbB,hTbF));

Print("\n-- dirTaFib --\n");
probe("TaF := r*x     (inverted)        ", mkRels(hAs,hBy,hBs, hTaB, r*x,      hTbB,hTbF));
probe("TaF := (rx)^-2 (doubled)         ", mkRels(hAs,hBy,hBs, hTaB, (r*x)^-2, hTbB,hTbF));

Print("\n-- dirTbBase (the repaired word) --\n");
probe("TbB := B            (pre-repair) ", mkRels(hAs,hBy,hBs, hTaB,hTaF, B,             hTbF));
probe("TbB := r^-1*M*r*B   (sign flip)  ", mkRels(hAs,hBy,hBs, hTaB,hTaF, r^-1*M*r*B,    hTbF));
probe("TbB := B*r^-1*M^-1*r (post)      ", mkRels(hAs,hBy,hBs, hTaB,hTaF, B*r^-1*M^-1*r, hTbF));
probe("TbB := x*M^-1*x^-1*B (other arc) ", mkRels(hAs,hBy,hBs, hTaB,hTaF, x*M^-1*x^-1*B, hTbF));

Print("\n-- dirTbFib --\n");
probe("TbF := r^-1    (unwhiskered)     ", mkRels(hAs,hBy,hBs, hTaB,hTaF, hTbB, r^-1));
probe("TbF := s*r*s^-1 (core inverted)  ", mkRels(hAs,hBy,hBs, hTaB,hTaF, hTbB, s*r*s^-1));

Print("\n-- corrections --\n");
probe("By: correction dropped           ", mkRels(hAs, B*y*B^-1*(y*x)^-1, hBs, hTaB,hTaF, hTbB,hTbF));
probe("Bs: arc delta'=x (legit)         ", mkRels(hAs, hBy, B*s*B^-1*(x*M*x^-1*s)^-1, hTaB,hTaF, hTbB,hTbF));
probe("As: correction dropped           ", mkRels(A*s*A^-1*y^-1, hBy, hBs, hTaB,hTaF, hTbB,hTbF));
QUIT_GAP(0);
