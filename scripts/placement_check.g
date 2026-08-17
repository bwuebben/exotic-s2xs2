# placement_check.g — placement-robustness check, updated 2026-08-15.
# Probes the transport-annulus orientation convention:
# each correction may multiply on the LEFT or the RIGHT of the transported word
# (the sign sweep does not cover this). Runs the LP cell (m=0, n=0) over all
# 2^3 placements x 2^5 signs = 256 cases.  The true drilled-fiber relation R3
# is included, as it is in the relation system used in the proof; coset
# enumeration is followed by kbmag certification on any overflow.

LoadPackage("kbmag");

F := FreeGroup("x","y","r","s","A","B","M","N");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;;
A := F.5;; B := F.6;; M := F.7;; N := F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
R0 := comm(x,y)*comm(r,s);;
base := [ R0,
  A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
  B*x*B^-1*y,    B*r*B^-1*r^-1,
  B*(s^-1*r^-1*y*x)*B^-1*(r^-1*s^-1*x)^-1 ];;
delta := r^-1;;
dirTa := A*x;;             # m=n=0: y1-side base directions only
# LOGGED FIX 2026-07-15 (paper: the pushoff-basing correction): honest
# dirTbBase, sign anti-coupled to e5.
dirTb := function(e5) return delta*M^(-e5)*delta^-1 * B; end;;

mk := function(conj, img, corr, side)   # relation  conj = [corr.]img[.corr]
  if side = 1 then return conj * (corr*img)^-1;
  else return conj * (img*corr)^-1; fi;
end;;

certify := function(G)
  local H, rws, kb;
  H := Image(IsomorphismSimplifiedFpGroup(G));
  rws := KBMAGRewritingSystem(H);
  kb := CALL_WITH_CATCH(KnuthBendix, [rws]);
  if kb[1] <> true or not IsConfluent(rws) then return false; fi;
  return Size(rws) = 1 and
    ForAll(GeneratorsOfGroup(FreeStructureOfRewritingSystem(rws)),
           z -> IsOne(ReducedForm(rws, z)));
end;;

cnt := rec(enum:=0, kb:=0, inconclusive:=0, h1:=0, fin:=0);;
for pAs in [1,2] do for pBy in [1,2] do for pBs in [1,2] do
  for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
    for eA in [1,-1] do for eB in [1,-1] do
      rels := Concatenation(base,
        [ mk(A*s*A^-1, y,   N^e3,               pAs),
          mk(B*y*B^-1, y*x, M^e4,               pBy),
          mk(B*s*B^-1, s,   delta*M^e5*delta^-1, pBs),
          M*dirTa^eA, N*dirTb(e5)^eB ]);
      G := F / rels;
      ab := AbelianInvariants(G);
      if Length(ab) > 0 then cnt.h1 := cnt.h1 + 1;
      else
        tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
                 RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
        if tab = fail then
          if certify(G) then cnt.kb := cnt.kb + 1;
          else cnt.inconclusive := cnt.inconclusive + 1; fi;
        elif Length(tab[1]) = 1 then cnt.enum := cnt.enum + 1;
        else cnt.fin := cnt.fin + 1;
          Print("FINITE>1 at p=(",pAs,",",pBy,",",pBs,") e=(",e3,",",e4,",",e5,
                ",",eA,",",eB,") |G|=",Length(tab[1]),"\n");
        fi;
      fi;
    od; od;
  od; od; od;
od; od; od;
Print("PLACEMENT SWEEP WITH R3 (m=0,n=0): ENUM_TRIVIAL=", cnt.enum,
      " KB_TRIVIAL=", cnt.kb, " INCONCLUSIVE=", cnt.inconclusive,
      " FINITE>1=", cnt.fin, " H1nonzero=", cnt.h1, "\n");
QUIT_GAP(0);
