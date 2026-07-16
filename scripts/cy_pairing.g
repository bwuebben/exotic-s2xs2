# cy_pairing.g — the c_y-resolved arc-pairing instrument 
# (pre-registered design: see the paper, Seifert-calibration section).
#
# Model: puncture the fiber AT the degenerate basing point c_y; split the
# y-generator at the puncture: Y = the reference side-perturbation, the
# other side = M^{e0} Y (junction meridian absorbed by the y1-whisker).
# pi1(F - c_y) = F4<x,r,s,Y>; the meridian M is a rotation-family word.
# The surgered manifolds are the SAME as the main run (lambda through c_y
# is isotopic to lambda through zeta0), so the ground-truth columns are
# rebuilt with the identical code and parameters.
#
# The pairing question: which side-choice in Theta(s) (clean Y vs decorated
# M^{e0} Y) is coherent with which drift word dir_lambda in {A r^-1, A x}.
# Pre-registered expectation E2: decorated <-> A x, clean <-> A r^-1;
# cross-pairings match nothing.
# Run: gap -q -A cy_pairing.g > cy_out.txt

start := Runtime();
Elapsed := function() return Concatenation("[", String(Int((Runtime()-start)/1000)), "s]"); end;
LI_DEPTH := 6;
FP := function(G)
  local fp, i, subs, H, inv;
  fp := [ SortedList(AbelianInvariants(G)) ];
  for i in [2..LI_DEPTH] do
    subs := LowIndexSubgroupsFpGroup(G, i);
    inv := [];
    for H in subs do
      if Index(G, H) = i then Add(inv, SortedList(AbelianInvariants(H))); fi;
    od;
    Add(fp, SortedList(inv));
  od;
  return fp;
end;
comm := function(u,v) return u*v*u^-1*v^-1; end;

# ---------------- truth columns (identical to the main run) ----------------
Print("=== truth columns (rebuilt, lm = k3.r-1) ", Elapsed(), " ===\n");
FT := FreeGroup("x","r","k3","A","Y");;
xT:=FT.1;; rT:=FT.2;; kT:=FT.3;; AT:=FT.4;; YT:=FT.5;;
Hrels := [ AT*xT*AT^-1*rT^-1, AT*rT*AT^-1*xT^-1, AT*kT*AT^-1*kT ];;
lpT := rT^-1;; lmT := kT*rT^-1;;
GTrue := function(k, n, d)
  local a1, b1, a2, b2;
  a1 := 1 + d*k*n*(1+2*n);  b1 := -d*k*n^2;
  a2 := d*k*(1+2*n)^2;      b2 := 1 - d*k*(1+2*n)*n;
  return FT / Concatenation(Hrels,
    [ YT*(AT*lpT)*YT^-1 * ((AT*lmT)^a1*AT^(2*b1))^-1,
      YT*AT^2*YT^-1 * ((AT*lmT)^a2*AT^(2*b2))^-1 ]);
end;;
TruthFP := [];;
for k in [1,-1,2] do for n in [-1,0,1] do for d in [1,-1] do
  Add(TruthFP, rec(k:=k, n:=n, d:=d, fp:=FP(GTrue(k,n,d))));
od; od; od;
Print("truth done ", Elapsed(), "\n");

# ---------------- the c_y-resolved fiber model ----------------
F4 := FreeGroup("x","r","s","Y");;
x4:=F4.1;; r4:=F4.2;; s4:=F4.3;; Y4:=F4.4;;
# meridian candidates: rotations of ([x,Y][r,s])^{+-1}
L := [x4, Y4, x4^-1, Y4^-1, r4, s4, r4^-1, s4^-1];;
Rot := function(j) local w,i; w:=One(F4);
  for i in [j..8] do w:=w*L[i]; od; for i in [1..j-1] do w:=w*L[i]; od;
  return w; end;;

# generous conjugator set for the transport-law pruner (length <= 2)
lets := [ One(F4), x4, r4, s4, Y4, x4^-1, r4^-1, s4^-1, Y4^-1 ];;
conjset := Set(Concatenation(lets, ListX(lets, lets, \*)));;

Print("\n=== coherence filter: (M-rotation, e0, Theta(s)-form) ", Elapsed(), " ===\n");
# Theta(s)-forms: clean Y; decorated M^{e}Y and YM^{e} (both placements).
Coherent := [];;
for j in [1..8] do for eM in [1,-1] do
  Mw := Rot(j)^eM;
  forms := [ rec(name:="clean Y",   w:=Y4,        dec:=false),
             rec(name:="M.Y",       w:=Mw*Y4,     dec:=true),
             rec(name:="M^-1.Y",    w:=Mw^-1*Y4,  dec:=true),
             rec(name:="Y.M",       w:=Y4*Mw,     dec:=true),
             rec(name:="Y.M^-1",    w:=Y4*Mw^-1,  dec:=true) ];
  for f in forms do
    Th := GroupHomomorphismByImages(F4, F4, [x4,r4,s4,Y4], [r4, x4, f.w, s4]);
    # NB images: Theta(x)=r, Theta(r)=x, Theta(s)=f.w, Theta(Y)=s
    ThM := Image(Th, Mw);
    hits := [];
    for dh in conjset do for eta in [1,-1] do
      if ThM = dh*Mw^eta*dh^-1 then Add(hits, [dh, eta]); fi;
    od; od;
    if Length(hits) > 0 then
      lasso := First(hits, h -> h[1] in [r4^-1, x4]);
      Print("COHERENT rot=", j, " eM=", eM, " Theta(s)=", f.name,
            "  #(dh,eta)=", Length(hits));
      if lasso <> fail then Print("  lasso-shaped dh=", lasso[1], " eta=", lasso[2]); fi;
      Print("\n");
      Add(Coherent, rec(j:=j, eM:=eM, form:=f, Mw:=Mw));
    fi;
  od;
od; od;
if Length(Coherent) = 0 then
  Print("### STOP: no coherent package at the c_y-basing ###\n"); QUIT_GAP(1);
fi;
Print(Length(Coherent), " coherent combos ", Elapsed(), "\n");

# ---------------- group stage: the pairing table ----------------
Print("\n=== pairing stage: coherent combos x dir in {A.r-1, A.x} x k ", Elapsed(), " ===\n");
FO := FreeGroup("x","r","s","Y","A");;
xO:=FO.1;; rO:=FO.2;; sO:=FO.3;; YO:=FO.4;; AO:=FO.5;;
toFO := GroupHomomorphismByImages(F4, FO, [x4,r4,s4,Y4], [xO,rO,sO,YO]);;
dirs := [ ["A.r-1", AO*rO^-1], ["A.x", AO*xO] ];;

Verdict := function(name, rels, k)
  local G, ab, fpG, hits, t;
  G := FO / rels;
  ab := SortedList(AbelianInvariants(G));
  fpG := FP(G);
  hits := Set(List(Filtered(TruthFP, t -> t.k = k and t.fp = fpG), t -> t.n));
  Print(name, ": H1=", ab, "  matches n in ", hits, "  ", Elapsed(), "\n");
  return hits;
end;;

# deduplicate combos by (Theta(s) word, M-word) at the presentation level
seen := [];;
rows := [];;
for cmb in Coherent do
  key := [ Image(toFO, cmb.form.w), Image(toFO, cmb.Mw) ];
  if key in seen then continue; fi;
  Add(seen, key);
  base := [ AO*xO*AO^-1*rO^-1, AO*rO*AO^-1*xO^-1, AO*YO*AO^-1*sO^-1,
            AO*sO*AO^-1*Image(toFO, cmb.form.w)^-1 ];
  MwO := Image(toFO, cmb.Mw);
  for dpair in dirs do
    ok := [];
    for k in [1,-1,2] do
      hits := Verdict(Concatenation("rot", String(cmb.j), " eM=", String(cmb.eM),
                " Th(s)=", cmb.form.name, " dir=", dpair[1], " k=", String(k)),
                Concatenation(base, [ MwO*dpair[2]^k ]), k);
      Add(ok, 0 in hits);
    od;
    Add(rows, rec(j:=cmb.j, eM:=cmb.eM, form:=cmb.form.name, dec:=cmb.form.dec,
                  dir:=dpair[1], pass:=ForAll(ok, IsBool) and ForAll(ok, b -> b = true)));
  od;
od;

Print("\n=== PAIRING TABLE (pass = matches the n=0 truth class at all k) ===\n");
for r in rows do
  Print("  Theta(s) ", String(r.form, -8), "  dir ", String(r.dir, -6),
        "  rot=", r.j, " eM=", String(r.eM, 2), "   ",
        r.pass, "\n");
od;
cleanAr  := ForAny(rows, r -> (not r.dec) and r.dir = "A.r-1" and r.pass);;
cleanAx  := ForAny(rows, r -> (not r.dec) and r.dir = "A.x"   and r.pass);;
decAr    := ForAny(rows, r -> r.dec and r.dir = "A.r-1" and r.pass);;
decAx    := ForAny(rows, r -> r.dec and r.dir = "A.x"   and r.pass);;
Print("\nSUMMARY:  clean<->A.r-1: ", cleanAr, "   decorated<->A.x: ", decAx,
      "\n          clean<->A.x  : ", cleanAx, "   decorated<->A.r-1: ", decAr, "\n");
if cleanAr and decAx and not cleanAx and not decAr then
  Print("E2 CONFIRMED: the pairing is {clean s-transport <-> A.r-1,\n");
  Print("decorated <-> A.x}, cross-pairings excluded. The display word Ar^-1\n");
  Print("is the certified coherent partner of the clean-transport convention.\n");
else
  Print("E2 NOT in the pre-registered pattern -- interpret per the pre-registered outcome branches.\n");
fi;
Print("\nDONE ", Elapsed(), "\n");
QUIT_GAP(0);
