F := FreeGroup("x","y","r","s","A","B");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;; A := F.5;; B := F.6;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
R0 := comm(x,y)*comm(r,s);;
Ax := A*x*A^-1*r^-1;;  Ar := A*r*A^-1*x^-1;;
Bx := B*x*B^-1*y;;     Br := B*r*B^-1*r^-1;;
Ay := A*y*A^-1*s^-1;;  As := A*s*A^-1*y^-1;;
By := B*y*B^-1*x^-1*y^-1;;  Bs := B*s*B^-1*s^-1;;
S1 := function(z, e) return comm(z, A)^e * B; end;;
S2y := function(e) return comm(y, B)^e * A; end;;
S2s := function(e) return comm(s, B)^e * A; end;;
quick := function(name, rels)
  local G, ab, tab, sz;
  G := F / rels;
  ab := AbelianInvariants(G);
  if Length(ab) > 0 then Print(name, " | H1=", ab, "\n"); return "H1"; fi;
  tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
           RelatorsOfFpGroup(G), [] : max := 150000, silent := true);
  if tab <> fail then
    sz := Length(tab[1]); Print(name, " | |G|=", sz, "\n");
    if sz = 1 then return "TRIVIAL"; else return "FINITE"; fi;
  fi;
  Print(name, " | enum-blowup (not visibly trivial)\n"); return "BLOWUP";
end;;
Print("=== CONTROLS ===\n");
quick("pi1(R)", [R0,Ax,Ay,Ar,As,Bx,By,Br,Bs]);
quick("R+S1", [R0,Ax,Ay,Ar,As,Bx,By,Br,Bs,S1(y*s,1)]);
quick("R+S2", [R0,Ax,Ay,Ar,As,Bx,By,Br,Bs,S2s(1)]);
Print("=== A: corrected dual [s,B] ===\n");
for e1 in [1,-1] do for e2 in [1,-1] do
  quick(Concatenation("LPs (",String(e1),",",String(e2),")"),
        [R0,Ax,Ay,Ar,As,Bx,By,Br,Bs,S1(y*s,e1),S2s(e2)]);
od; od;
Print("=== B: G_clean (broken relations dropped) ===\n");
for e1 in [1,-1] do for e2 in [1,-1] do
  quick(Concatenation("Gclean-s (",String(e1),",",String(e2),")"),
        [R0,Ax,Ar,Bx,Br,S1(y*s,e1),S2s(e2)]);
od; od;
Print("=== C: correction sensitivity, 120 trials ===\n");
corrTa := [ x*r, (x*r)^-1, A, A^-1 ];;
corrTb := [ r, r^-1, B, B^-1 ];;
conjs := [ One(F), x, y, r, s, A, B, x*y, r*s ];;
cnt := rec(TRIVIAL:=0, FINITE:=0, H1:=0, BLOWUP:=0);;
rs := RandomSource(IsMersenneTwister, 20260709);;
for trial in [1..120] do
  u1 := corrTa[Random(rs,1,4)]; g1 := conjs[Random(rs,1,9)];
  u2 := corrTb[Random(rs,1,4)]; g2 := conjs[Random(rs,1,9)];
  u3 := corrTa[Random(rs,1,4)]; g3 := conjs[Random(rs,1,9)];
  u4 := corrTa[Random(rs,1,4)]; g4 := conjs[Random(rs,1,9)];
  Byc := B*y*B^-1 * (g1*u1*g1^-1) * y^-1*x^-1;
  Bsc := B*s*B^-1 * (g2*u2*g2^-1) * s^-1;
  Ayc := A*y*A^-1 * (g3*u3*g3^-1) * s^-1;
  Asc := A*s*A^-1 * (g4*u4*g4^-1) * y^-1;
  e1 := Random(rs,[1,-1]); e2 := Random(rs,[1,-1]);
  res := quick(Concatenation("C#",String(trial)),
     [R0,Ax,Ar,Bx,Br,Byc,Bsc,Ayc,Asc,S1(y*s,e1),S2s(e2)]);
  cnt.(res) := cnt.(res) + 1;
od;
Print("SUMMARY C: TRIVIAL=", cnt.TRIVIAL, " FINITE=", cnt.FINITE,
      " H1=", cnt.H1, " BLOWUP=", cnt.BLOWUP, "\n");
QUIT_GAP(0);
