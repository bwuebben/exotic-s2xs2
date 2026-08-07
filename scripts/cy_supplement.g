# cy_supplement.g — arbiter supplement to cy_pairing.g (design-slip fix;
# the main run let the coherence pruner (conjugator
# search capped at length 2) exclude the clean-Theta(s) rows before they
# reached the truth-match, although the pre-registration named the
# truth-match as the arbiter.  The clean-swap transport law needs a longer
# (cyclic-rotation) conjugator, so its exclusion could be an artifact.
# Here the clean rows for the two whisker-relevant meridian families
# (rot 2 = V2/y1-side, rot 5 = V5/y2-side) are sent to the arbiter directly.
# Run: gap -q -A cy_supplement.g > cy_supp_out.txt

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
  Add(TruthFP, rec(k:=k, n:=n, fp:=FP(GTrue(k,n,d))));
od; od; od;
Print("truth rebuilt ", Elapsed(), "\n");

F4 := FreeGroup("x","r","s","Y");;
x4:=F4.1;; r4:=F4.2;; s4:=F4.3;; Y4:=F4.4;;
L := [x4, Y4, x4^-1, Y4^-1, r4, s4, r4^-1, s4^-1];;
Rot := function(j) local w,i; w:=One(F4);
  for i in [j..8] do w:=w*L[i]; od; for i in [1..j-1] do w:=w*L[i]; od;
  return w; end;;

# also: does the clean swap satisfy the transport law with LONGER conjugators?
lets := [ One(F4), x4, r4, s4, Y4, x4^-1, r4^-1, s4^-1, Y4^-1 ];;
c2 := Set(Concatenation(lets, ListX(lets, lets, \*)));;
c4 := Set(Concatenation(c2, ListX(c2, c2, \*)));;   # length <= 4
for j in [2, 5] do for eM in [1,-1] do
  Mw := Rot(j)^eM;
  Th := GroupHomomorphismByImages(F4, F4, [x4,r4,s4,Y4], [r4, x4, Y4, s4]);
  ThM := Image(Th, Mw);
  hits := Filtered(c4, w -> ThM = w*Mw*w^-1 or ThM = w*Mw^-1*w^-1);
  Print("clean-swap transport, rot=", j, " eM=", eM,
        ": conjugators (len<=4) found: ", Length(hits) > 0);
  if Length(hits) > 0 then Print("   e.g. ", hits[1]); fi;
  Print("\n");
od; od;

Print("\n=== arbiter: clean-Theta(s) rows, rot in {2,5} ", Elapsed(), " ===\n");
FO := FreeGroup("x","r","s","Y","A");;
xO:=FO.1;; rO:=FO.2;; sO:=FO.3;; YO:=FO.4;; AO:=FO.5;;
toFO := GroupHomomorphismByImages(F4, FO, [x4,r4,s4,Y4], [xO,rO,sO,YO]);;
base := [ AO*xO*AO^-1*rO^-1, AO*rO*AO^-1*xO^-1, AO*YO*AO^-1*sO^-1,
          AO*sO*AO^-1*YO^-1 ];;   # clean: Theta(s) = Y
dirs := [ ["A.r-1", AO*rO^-1], ["A.x", AO*xO] ];;
for j in [2, 5] do for eM in [1,-1] do
  MwO := Image(toFO, Rot(j)^eM);
  for dpair in dirs do
    passall := true;
    for k in [1,-1,2] do
      G := FO / Concatenation(base, [ MwO*dpair[2]^k ]);
      fpG := FP(G);
      hits := Set(List(Filtered(TruthFP, t -> t.k = k and t.fp = fpG), t -> t.n));
      Print("clean rot", j, " eM=", String(eM,2), " dir=", dpair[1],
            " k=", String(k,2), ": H1=", SortedList(AbelianInvariants(G)),
            "  matches n in ", hits, "  ", Elapsed(), "\n");
      if not 0 in hits then passall := false; fi;
    od;
    Print("  => clean rot", j, " eM=", String(eM,2), " dir=", dpair[1],
          "  ALL-k n=0 match: ", passall, "\n");
  od;
od; od;
Print("\nDONE ", Elapsed(), "\n");
QUIT_GAP(0);
