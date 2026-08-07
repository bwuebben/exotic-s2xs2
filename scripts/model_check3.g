# Hostile-verification checks for memo Section 2 claims NOT covered by monodromy_check2.g
# (written 2026-07-09 for the arXiv-note verification pass).
#
# Conventions: comm(u,v) = u v u^-1 v^-1 (as in all project scripts).
# CompositionMapping(g,f) in GAP is g after f: x |-> g(f(x)).

# ---------- (1) The trefoil monodromy is the composite of the stated twists ----------
F2 := FreeGroup("x","y");;
x := F2.1;; y := F2.2;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
Ta := GroupHomomorphismByImages(F2, F2, [x,y], [x, y*x]);;      # twist along a ~ x
Tb := GroupHomomorphismByImages(F2, F2, [x,y], [x*y^-1, y]);;   # twist along b ~ y
h  := GroupHomomorphismByImages(F2, F2, [x,y], [y^-1, y*x]);;
TaTb := CompositionMapping(Ta, Tb);;   # Tb first, then Ta
TbTa := CompositionMapping(Tb, Ta);;   # Ta first, then Tb
Print("h = Ta o Tb (Tb first): ",
  (Image(TaTb,x)=Image(h,x)) and (Image(TaTb,y)=Image(h,y)), "\n");
Print("h = Tb o Ta (Ta first): ",
  (Image(TbTa,x)=Image(h,x)) and (Image(TbTa,y)=Image(h,y)), "\n");
# h^-1 formulas used in the note: x |-> xy, y |-> x^-1
hinv := GroupHomomorphismByImages(F2, F2, [x,y], [x*y, x^-1]);;
Print("h o hinv = id and hinv o h = id: ",
  (Image(CompositionMapping(h,hinv),x)=x) and (Image(CompositionMapping(h,hinv),y)=y) and
  (Image(CompositionMapping(hinv,h),x)=x) and (Image(CompositionMapping(hinv,h),y)=y), "\n");

# ---------- (2) Relator preservation by the two structural automorphisms ----------
F4 := FreeGroup("xx","yy","rr","ss");;
xx := F4.1;; yy := F4.2;; rr := F4.3;; ss := F4.4;;
rel := comm(xx,yy)*comm(rr,ss);;
phi := GroupHomomorphismByImages(F4, F4, [xx,yy,rr,ss], [rr,ss,xx,yy]);;            # swap
psi := GroupHomomorphismByImages(F4, F4, [xx,yy,rr,ss], [yy^-1, yy*xx, rr, ss]);;    # h * id
Print("psi(rel) = rel exactly: ", Image(psi,rel) = rel, "\n");
Print("phi(rel) = [xx,yy]^-1 rel [xx,yy] (conjugate): ",
  Image(phi,rel) = comm(xx,yy)^-1 * rel * comm(xx,yy), "\n");

# ---------- (3) The block-commutator identity: psi phi psi^-1 phi = h * h^-1 ----------
psiinv := GroupHomomorphismByImages(F4, F4, [xx,yy,rr,ss], [xx*yy, xx^-1, rr, ss]);;
lhs := CompositionMapping(psi, CompositionMapping(phi, CompositionMapping(psiinv, phi)));;
Hclosed := GroupHomomorphismByImages(F4, F4, [xx,yy,rr,ss], [yy^-1, yy*xx, rr*ss, rr^-1]);;  # h * h^-1
ok := true;;
for g in [xx,yy,rr,ss] do
  if Image(lhs,g) <> Image(Hclosed,g) then ok := false; fi;
od;
Print("psi o phi o psi^-1 o phi = h*h^-1 on all generators: ", ok, "\n");
Print("  images: ", List([xx,yy,rr,ss], g -> Image(lhs,g)), "\n");
# sanity: phi is an involution
Print("phi^2 = id: ",
  ForAll([xx,yy,rr,ss], g -> Image(CompositionMapping(phi,phi),g) = g), "\n");

# ---------- (4) Fiber H1 dies in H1(rr): coinvariants of <Phi,Psi> on Z^4 are 0 ----------
# Basis (x,y,r,s); rows are images. Phi = swap; Psi = H (+) I with H = [[0,1],[-1,1]]
# (row-vector convention: x |-> -y gives row (0,-1,0,0); y |-> y+x gives (1,1,0,0)).
Id4 := IdentityMat(4);;
PhiM := [[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]];;
PsiM := [[0,-1,0,0],[1,1,0,0],[0,0,1,0],[0,0,0,1]];;
Mstack := Concatenation(PhiM - Id4, PsiM - Id4);;
sm := SmithNormalFormIntegerMat(Mstack);;
diag := List([1..4], i -> sm[i][i]);;
Print("Smith normal form diagonal of (Phi-I; Psi-I): ", diag,
      "  => coinvariants trivial: ", diag = [1,1,1,1], "\n");

QUIT_GAP(0);
