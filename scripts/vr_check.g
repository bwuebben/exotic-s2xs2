# vr_check.g — V-R style falsification probes with KNOWN answers (logged).
F := FreeGroup("x","y","r","s","A","B","M","N");;
x:=F.1;;y:=F.2;;r:=F.3;;s:=F.4;;A:=F.5;;B:=F.6;;M:=F.7;;N:=F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
R0 := comm(x,y)*comm(r,s);;
base := [R0, A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1, B*x*B^-1*y, B*r*B^-1*r^-1];;
corr := [A*s*A^-1*(N*y)^-1, B*y*B^-1*(M*y*x)^-1, B*s*B^-1*(r^-1*M*r*s)^-1];;
Print("no fillings (expect H1 = Z^2):        ",
  AbelianInvariants(F/Concatenation(base,corr)), "\n");
Print("fiber-only fillings (expect H1 = Z^2): ",
  AbelianInvariants(F/Concatenation(base,corr,
    [M*((r*x)^-1)^1, N*(s*r^-1*s^-1)^1])), "\n");
Print("one base filling (expect H1 = Z):      ",
  AbelianInvariants(F/Concatenation(base,corr, [N*(r^-1*M^-1*r*B)])), "\n");
# (2026-07-15, paper: the pushoff-basing correction — honest dirTbBase at the
#  corr sign e5 = +1 used above.)
QUIT_GAP(0);
