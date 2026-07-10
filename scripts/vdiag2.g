# vdiag2.g — V-DIAG second diagram (2026-07-09).
# Diagram 2 = detour representative of beta-bar; only the Bs relation changes:
#   B s B^-1 = N^e6 . (r^-1 M^e5 r) . s     (extra N-term, lasso = N's whisker;
# M-part unchanged since the detour-transport loop D(u_c) = 1 by the position
# table u_c < u_{s_e}; the C2 crossing vanishes since s . phi0(e) = s . a = 0).
F := FreeGroup("x","y","r","s","A","B","M","N");;
x:=F.1;;y:=F.2;;r:=F.3;;s:=F.4;;A:=F.5;;B:=F.6;;M:=F.7;;N:=F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
R0 := comm(x,y)*comm(r,s);;
base := [R0, A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1, B*x*B^-1*y, B*r*B^-1*r^-1];;
dirTaBase := A*r^-1;;  dirTaFib := (r*x)^-1;;
dirTbBase := B;;       dirTbFib := s*r^-1*s^-1;;
cnt := rec(triv:=0, blow:=0, h1:=0, fin:=0);;
for mn in Cartesian([-1,0,1],[-1,0,1]) do
  ct := rec(triv:=0, blow:=0, h1:=0, fin:=0);
  for eA in [1,-1] do for eB in [1,-1] do
    for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do for e6 in [1,-1] do
      rels := Concatenation(base,
        [ A*s*A^-1*(N^e3*y)^-1,
          B*y*B^-1*(M^e4*y*x)^-1,
          B*s*B^-1*(N^e6 * r^-1*M^e5*r * s)^-1,
          M*(dirTaBase*dirTaFib^mn[2])^eA,
          N*(dirTbBase*dirTbFib^mn[1])^eB ]);
      G := F / rels;
      ab := AbelianInvariants(G);
      if Length(ab) > 0 then ct.h1 := ct.h1 + 1;
      else
        tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
                 RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
        if tab = fail then ct.blow := ct.blow + 1;
        elif Length(tab[1]) = 1 then ct.triv := ct.triv + 1;
        else ct.fin := ct.fin + 1;
          Print("FINITE>1 at mn=",mn," |G|=",Length(tab[1]),"\n"); fi;
      fi;
    od; od; od; od;
  od; od;
  Print("DIAG2 m=",mn[1]," n=",mn[2],": triv=",ct.triv," blow=",ct.blow,
        " fin=",ct.fin," h1=",ct.h1,"\n");
  cnt.triv := cnt.triv+ct.triv; cnt.blow := cnt.blow+ct.blow;
  cnt.fin := cnt.fin+ct.fin;   cnt.h1 := cnt.h1+ct.h1;
od;
Print("DIAG2 TOTAL: triv=",cnt.triv," blow=",cnt.blow," fin=",cnt.fin,
      " h1=",cnt.h1,"\n");
QUIT_GAP(0);
