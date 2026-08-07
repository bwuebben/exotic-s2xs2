# confirm_coherent.g — post-run localization follow-up (see the paper: Seifert calibration).
# The main run (seifert_out.txt) showed: P1-package + dir A.r^-1 matches NO
# truth column, while the SAME package + dir A.x (labelled W3) matches the
# {n=0, n=-1} fp-class at every k. Localization: dir_lambda must be whisker-
# COHERENT with the zeta0-basing of the Theta-package; the coherent P1 word
# is A.x (gamma-route), and by the rho-mirror the coherent P2 word is A.r.
# Predictions tested here (pre-stated):
#   C1: P2-package + dir A.r      -> matches the same {n=0,-1} class, all k.
#   C2: P1-package combo 2 + A.x  -> matches (sign-pairing invariance).
#   C3: P2-package + dir A.x      -> matches nothing (incoherent mix control).
# Run: gap -q -A confirm_coherent.g > confirm_out.txt

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

# --- truth reference columns (lm = k3.r-1 representative; n in {-1,0,1}, d=1;
# the main run showed all lm/d give the same fp within each (k, n-class)) ---
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
for k in [1,-1,2] do for n in [-1,0,1] do
  Add(TruthFP, rec(k:=k, n:=n, fp:=FP(GTrue(k,n,1))));
  Print("truth k=", k, " n=", n, " ", Elapsed(), "\n");
od; od;

# --- packages (coherent family from the main run: rot=2, (eM,es) anti-coupled,
# pl=L; combo1 = (eM=1,es=-1), combo2 = (eM=-1,es=1)) ---
F4 := FreeGroup("x","y","r","s");;
x4:=F4.1;; y4:=F4.2;; r4:=F4.3;; s4:=F4.4;;
L := [x4, y4, x4^-1, y4^-1, r4, s4, r4^-1, s4^-1];;
Rot := function(j) local w,i; w:=One(F4);
  for i in [j..8] do w:=w*L[i]; od; for i in [1..j-1] do w:=w*L[i]; od;
  return w; end;;
FO := FreeGroup("x","y","r","s","A");;
xO:=FO.1;; yO:=FO.2;; rO:=FO.3;; sO:=FO.4;; AO:=FO.5;;
toFO := GroupHomomorphismByImages(F4, FO, [x4,y4,r4,s4], [xO,yO,rO,sO]);;
swapF4 := GroupHomomorphismByImages(F4, F4, [x4,y4,r4,s4], [r4,s4,x4,y4]);;

Case := function(name, rels, k)
  local G, hits, t;
  G := FO / rels;
  hits := [];
  for t in TruthFP do
    if t.k = k and t.fp = FP(G) then Add(hits, t.n); fi;
  od;
  Print(name, ": H1=", SortedList(AbelianInvariants(G)),
        "  matches n=", hits, "  ", Elapsed(), "\n");
end;;

for combo in [ [1,-1], [-1,1] ] do   # (eM, es)
  Mw := Rot(2)^combo[1];;
  MwO := Image(toFO, Mw);;
  Mw2O := Image(toFO, Image(swapF4, Mw));;
  P1r := [ AO*xO*AO^-1*rO^-1, AO*rO*AO^-1*xO^-1, AO*yO*AO^-1*sO^-1,
           AO*sO*AO^-1*(Image(toFO, Mw)^combo[2]*yO)^-1 ];;
  P2r := [ AO*xO*AO^-1*rO^-1, AO*rO*AO^-1*xO^-1, AO*sO*AO^-1*yO^-1,
           AO*yO*AO^-1*(Mw2O^combo[2]*sO)^-1 ];;
  for k in [1,-1,2] do
    Case(Concatenation("C1 P2+A.r    (eM=", String(combo[1]), ") k=", String(k)),
         Concatenation(P2r, [ Mw2O*(AO*rO)^k ]), k);
    Case(Concatenation("C2 P1+A.x    (eM=", String(combo[1]), ") k=", String(k)),
         Concatenation(P1r, [ MwO*(AO*xO)^k ]), k);
    Case(Concatenation("C3 P2+A.x ctl(eM=", String(combo[1]), ") k=", String(k)),
         Concatenation(P2r, [ Mw2O*(AO*xO)^k ]), k);
  od;
od;
Print("DONE ", Elapsed(), "\n");
QUIT_GAP(0);
