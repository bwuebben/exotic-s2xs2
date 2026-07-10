# placement_check.g — logged robustness EXTENSION (2026-07-09, session 5).
# Not part of the frozen harness. Probes the membrane-orientation convention debt:
# each correction may multiply on the LEFT or the RIGHT of the transported word
# (the sign sweep does not cover this). Runs the LP cell (m=0, n=0) over all
# 2^3 placements x 2^5 signs = 256 cases, enumeration-only verdicts.

F := FreeGroup("x","y","r","s","A","B","M","N");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;;
A := F.5;; B := F.6;; M := F.7;; N := F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
R0 := comm(x,y)*comm(r,s);;
base := [ R0,
  A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
  B*x*B^-1*y,    B*r*B^-1*r^-1 ];;
delta := r^-1;;
dirTa := A*r^-1;;          # m=n=0: base directions only
dirTb := B;;

mk := function(conj, img, corr, side)   # relation  conj = [corr.]img[.corr]
  if side = 1 then return conj * (corr*img)^-1;
  else return conj * (img*corr)^-1; fi;
end;;

cnt := rec(triv:=0, blow:=0, h1:=0, fin:=0);;
for pAs in [1,2] do for pBy in [1,2] do for pBs in [1,2] do
  for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
    for eA in [1,-1] do for eB in [1,-1] do
      rels := Concatenation(base,
        [ mk(A*s*A^-1, y,   N^e3,               pAs),
          mk(B*y*B^-1, y*x, M^e4,               pBy),
          mk(B*s*B^-1, s,   delta*M^e5*delta^-1, pBs),
          M*dirTa^eA, N*dirTb^eB ]);
      G := F / rels;
      ab := AbelianInvariants(G);
      if Length(ab) > 0 then cnt.h1 := cnt.h1 + 1;
      else
        tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
                 RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
        if tab = fail then cnt.blow := cnt.blow + 1;
        elif Length(tab[1]) = 1 then cnt.triv := cnt.triv + 1;
        else cnt.fin := cnt.fin + 1;
          Print("FINITE>1 at p=(",pAs,",",pBy,",",pBs,") e=(",e3,",",e4,",",e5,
                ",",eA,",",eB,") |G|=",Length(tab[1]),"\n");
        fi;
      fi;
    od; od;
  od; od; od;
od; od; od;
Print("PLACEMENT SWEEP (m=0,n=0): TRIVIAL=", cnt.triv, " BLOWUP=", cnt.blow,
      " FINITE>1=", cnt.fin, " H1nonzero=", cnt.h1, "\n");
QUIT_GAP(0);
