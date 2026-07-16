# decide_seifert.g — Seifert single-surgery drift calibration.
# Pre-registered design and derivation: see the paper, Seifert-calibration section.
# Structure:
#   Phase 0: reference pi1(N) (convention-free mapping torus) + fingerprint.
#   Phase 1: derive-by-constraint — pin (M-word rotation, eps_M, eps_s,
#            placement, delta-hat, eta) by the Theta(M)-transport law as
#            free-group identities. STOP if no coherent family.
#   Phase 2: ground-truth family G_true(lm; k, n, d) — HNN over
#            H = <x,r,k3,A | swap, A k3 A^-1 = k3^-1>, wall (A^2, A*l),
#            twist along sigma_n = (1+2n)*lam - n*f. GT-0 filter of the
#            lm-candidates at k=0; H1 checks (GT-1); column fingerprints.
#   Phase 3: our side G_ours(k) for coherent packages P1/P2 + probes +
#            wrong-word panel W1-W8; verdict = which truth column matches.
# Run: gap -q -A decide_seifert.g > seifert_out.txt

start := Runtime();
Elapsed := function() return Concatenation("[", String(Int((Runtime()-start)/1000)), "s]"); end;

LI_DEPTH := 6;   # depth 6 from the start (t4 deviation lesson; the pre-registered design)

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

comm := function(u,v) return u*v*u^-1*v^-1; end;

# ================= Phase 0: reference pi1(N) =================
Print("=== PHASE 0: reference pi1(N) ", Elapsed(), " ===\n");
FN := FreeGroup("x","y","r","s","A");;
xN:=FN.1;; yN:=FN.2;; rN:=FN.3;; sN:=FN.4;; AN:=FN.5;;
GN := FN / [ comm(xN,yN)*comm(rN,sN),
             AN*xN*AN^-1*rN^-1, AN*yN*AN^-1*sN^-1,
             AN*rN*AN^-1*xN^-1, AN*sN*AN^-1*yN^-1 ];;
Print("H1(pi1 N) = ", SortedList(AbelianInvariants(GN)), "  (expected [0,0,0])\n");
fpN := FP(GN);;
Print("fpN computed ", Elapsed(), "\n");

# ================= Phase 1: derive-by-constraint =================
# F4 = pi1(F minus puncture); M = v R0^e v^-1 = a cyclic rotation of R0^e.
Print("\n=== PHASE 1: coherence constraint (free-group identities) ", Elapsed(), " ===\n");
F4 := FreeGroup("x","y","r","s");;
x4:=F4.1;; y4:=F4.2;; r4:=F4.3;; s4:=F4.4;;
L := [x4, y4, x4^-1, y4^-1, r4, s4, r4^-1, s4^-1];;   # R0 letters
Rot := function(j)
  local w, i;
  w := One(F4);
  for i in [j..8] do w := w*L[i]; od;
  for i in [1..j-1] do w := w*L[i]; od;
  return w;
end;;
Survivors := [];;
for j in [1..8] do for eM in [1,-1] do
  Mw := Rot(j)^eM;
  for es in [1,-1] do for pl in ["L","R"] do
    if pl = "L" then ThS := Mw^es*y4; else ThS := y4*Mw^es; fi;
    hom := GroupHomomorphismByImages(F4, F4, [x4,y4,r4,s4], [r4, s4, x4, ThS]);
    ThM := Image(hom, Mw);
    for dh in [ r4^-1, r4^-1*Mw, r4^-1*Mw^-1, Mw*r4^-1, Mw^-1*r4^-1 ] do
      for eta in [1,-1] do
        if ThM = dh*Mw^eta*dh^-1 then
          Add(Survivors, rec(j:=j, eM:=eM, es:=es, pl:=pl, dh:=dh, eta:=eta, Mw:=Mw, hom:=hom));
          Print("COHERENT: rot=", j, " eM=", eM, " es=", es, " pl=", pl,
                " dh=", dh, " eta=", eta, "\n");
        fi;
      od;
    od;
  od; od;
od; od;
if Length(Survivors) = 0 then
  Print("### STOP: no coherent (M-word, signs) package found — derivation error. ###\n");
  QUIT_GAP(1);
fi;
Print(Length(Survivors), " coherent P1 combos ", Elapsed(), "\n");
# Theta^2 shape report for the first survivor: conjugator of Theta^2(M)
S1 := Survivors[1];;
Th2M := Image(S1.hom, Image(S1.hom, S1.Mw));;
Print("Theta^2(M) conjugator check: Theta^2(M) = (rx)^-1 M (rx)?  ",
      Th2M = (r4*x4)^-1*S1.Mw*(r4*x4), "   (or with eta twists — informational)\n");

# ================= Phase 2: ground-truth family =================
Print("\n=== PHASE 2: ground truth ", Elapsed(), " ===\n");
FT := FreeGroup("x","r","k3","A","Y");;
xT:=FT.1;; rT:=FT.2;; kT:=FT.3;; AT:=FT.4;; YT:=FT.5;;
Hrels := [ AT*xT*AT^-1*rT^-1, AT*rT*AT^-1*xT^-1, AT*kT*AT^-1*kT ];;
lpT := rT^-1;;
GTrue := function(lm, k, n, d)
  local a1, b1, a2, b2, w1, w2;
  a1 := 1 + d*k*n*(1+2*n);  b1 := -d*k*n^2;
  a2 := d*k*(1+2*n)^2;      b2 := 1 - d*k*(1+2*n)*n;
  w1 := (AT*lm)^a1 * AT^(2*b1);
  w2 := (AT*lm)^a2 * AT^(2*b2);
  return FT / Concatenation(Hrels,
    [ YT*(AT*lpT)*YT^-1 * w1^-1,  YT*AT^2*YT^-1 * w2^-1 ]);
end;;

lmCands := [ rT^-1, kT*rT^-1, rT^-1*kT, kT^-1*rT^-1, rT^-1*kT^-1 ];;
lmNames := [ "r-1(control)", "k3.r-1", "r-1.k3", "k3-1.r-1", "r-1.k3-1" ];;
# GT-0 filter: k=0 closure must be pi1(N). ab-prediction: odd-k3 candidates
# give [0,0,0]; the r-1 control keeps the Z/2 and must FAIL.
lmOK := [];;
for i in [1..Length(lmCands)] do
  G0 := GTrue(lmCands[i], 0, 0, 1);
  ab := SortedList(AbelianInvariants(G0));
  Print("GT-0 lm=", lmNames[i], ": H1=", ab);
  if ab = [0,0,0] then
    if FP(G0) = fpN then
      Print("  fp=pi1(N) MATCH");
      Add(lmOK, i);
    else
      Print("  fp MISMATCH");
    fi;
  fi;
  Print("  ", Elapsed(), "\n");
od;
if Length(lmOK) = 0 then
  Print("### STOP: no lm-candidate closes to pi1(N) — wall-word derivation error. ###\n");
  QUIT_GAP(1);
fi;
Print("GT-0 survivors: ", lmNames{lmOK}, "\n");

# Truth columns: fingerprints over the panel. GT-1: H1 must be Z^2 (k=+-1),
# Z^2+Z/2 (k=2) in EVERY column (n-blind, as pre-registered in the addendum).
TruthFP := [];;
for i in lmOK do
  for k in [1,-1,2] do for n in [-1,0,1] do for d in [1,-1] do
    G := GTrue(lmCands[i], k, n, d);
    ab := SortedList(AbelianInvariants(G));
    fpG := FP(G);
    Add(TruthFP, rec(lm:=i, k:=k, n:=n, d:=d, ab:=ab, fp:=fpG));
    Print("TRUTH lm=", lmNames[i], " k=", k, " n=", n, " d=", d,
          ": H1=", ab, "  ", Elapsed(), "\n");
  od; od; od;
od;
# Distinctness within each k (needed for discrimination claims):
for k in [1,-1,2] do
  cols := Filtered(TruthFP, t -> t.k = k);
  for i in [1..Length(cols)] do for j in [i+1..Length(cols)] do
    if cols[i].fp = cols[j].fp then
      Print("NOTE k=", k, ": columns (", lmNames[cols[i].lm], ",n=", cols[i].n,
            ",d=", cols[i].d, ") and (", lmNames[cols[j].lm], ",n=", cols[j].n,
            ",d=", cols[j].d, ") have EQUAL fp\n");
    fi;
  od; od;
od;
Print("truth columns done ", Elapsed(), "\n");

# ================= Phase 3: our side =================
Print("\n=== PHASE 3: pipeline side ", Elapsed(), " ===\n");
FO := FreeGroup("x","y","r","s","A");;
xO:=FO.1;; yO:=FO.2;; rO:=FO.3;; sO:=FO.4;; AO:=FO.5;;
# transfer F4-words into FO letters:
toFO := GroupHomomorphismByImages(F4, FO, [x4,y4,r4,s4], [xO,yO,rO,sO]);;
swapF4 := GroupHomomorphismByImages(F4, F4, [x4,y4,r4,s4], [r4,s4,x4,y4]);;

# Base Theta-relations for a package: P1 decorates s; P2 = swap-mirror
# (decorates y). rels(package) with filling word appended by caller.
P1rels := function(S)   # S = a coherent survivor record
  local ThS;
  if S.pl = "L" then ThS := S.Mw^S.es*y4; else ThS := y4*S.Mw^S.es; fi;
  return [ AO*xO*AO^-1*rO^-1, AO*rO*AO^-1*xO^-1, AO*yO*AO^-1*sO^-1,
           AO*sO*AO^-1*Image(toFO, ThS)^-1 ];
end;;
P2rels := function(S)
  local Mw2, ThY;
  Mw2 := Image(swapF4, S.Mw);
  if S.pl = "L" then ThY := Mw2^S.es*s4; else ThY := s4*Mw2^S.es; fi;
  return [ AO*xO*AO^-1*rO^-1, AO*rO*AO^-1*xO^-1, AO*sO*AO^-1*yO^-1,
           AO*yO*AO^-1*Image(toFO, ThY)^-1 ];
end;;

Verdict := function(name, rels, k)
  local G, ab, fpG, hits, t;
  G := FO / rels;
  ab := SortedList(AbelianInvariants(G));
  fpG := FP(G);
  hits := [];
  for t in TruthFP do
    if t.k = k and t.fp = fpG then
      Add(hits, [ lmNames[t.lm], t.n, t.d ]);
    fi;
  od;
  Print(name, ": H1=", ab, "  matches=", hits, "  ", Elapsed(), "\n");
  return hits;
end;;

# Probes (package = first coherent survivor):
S1 := Survivors[1];;
MwO := Image(toFO, S1.Mw);;
Print("PROBE complement (no filling): H1=",
      SortedList(AbelianInvariants(FO / P1rels(S1))), "  (expected [0,0,0])\n");
GkillM := FO / Concatenation(P1rels(S1), [MwO]);;
Print("PROBE kill-M: H1=", SortedList(AbelianInvariants(GkillM)),
      "  fp=pi1(N)? ", FP(GkillM) = fpN, "  (must be true)  ", Elapsed(), "\n");

# Honest sweep: coherent packages x k in {1,-1,2}; ours-dictionary fixed
# (d-ambiguity lives in the truth columns).
dirP1 := AO*rO^-1;;
dirP2 := AO*xO^-1;;
for si in [1..Minimum(2, Length(Survivors))] do
  S := Survivors[si];
  for k in [1,-1,2] do
    Verdict(Concatenation("OURS P1(combo ", String(si), ") k=", String(k)),
            Concatenation(P1rels(S), [ Image(toFO, S.Mw)*dirP1^k ]), k);
    Verdict(Concatenation("OURS P2(combo ", String(si), ") k=", String(k)),
            Concatenation(P2rels(S), [ Image(toFO, Image(swapF4, S.Mw))*dirP2^k ]), k);
  od;
od;

# Wrong-word panel: P1 package of survivor 1, ONE slot replaced.
Print("\n=== W-PANEL (single-slot errors inside the P1 package) ", Elapsed(), " ===\n");
Wdirs := [ [ "W1 dir=A.r",        AO*rO ],
           [ "W2 dir=A.x-1",      AO*xO^-1 ],
           [ "W3 dir=A.x",        AO*xO ],
           [ "W4 dir=A.y-1r-1y",  AO*yO^-1*rO^-1*yO ],
           [ "W5 dir=A.r-1[x,y]", AO*rO^-1*comm(xO,yO) ],
           [ "W8 dir=A.r-1.M",    AO*rO^-1*MwO ] ];;
for wd in Wdirs do
  for k in [1,-1,2] do
    Verdict(Concatenation(wd[1], " k=", String(k)),
            Concatenation(P1rels(S1), [ MwO*wd[2]^k ]), k);
  od;
od;
# W6: arc-incoherent Theta(s)-decoration (x-route conjugate), derived dir.
W6rels := [ AO*xO*AO^-1*rO^-1, AO*rO*AO^-1*xO^-1, AO*yO*AO^-1*sO^-1,
            AO*sO*AO^-1*( (xO^-1*MwO*xO)^S1.es * yO )^-1 ];;
for k in [1,-1,2] do
  Verdict(Concatenation("W6 Theta(s) x-route k=", String(k)),
          Concatenation(W6rels, [ MwO*dirP1^k ]), k);
od;
# W7: drift-blind monodromy (forgot the push entirely).
W7rels := [ AO*xO*AO^-1*rO^-1, AO*rO*AO^-1*xO^-1, AO*yO*AO^-1*sO^-1,
            AO*sO*AO^-1*yO^-1 ];;
for k in [1,-1,2] do
  Verdict(Concatenation("W7 drift-blind k=", String(k)),
          Concatenation(W7rels, [ MwO*dirP1^k ]), k);
od;

Print("\nDONE ", Elapsed(), "\n");
QUIT_GAP(0);
