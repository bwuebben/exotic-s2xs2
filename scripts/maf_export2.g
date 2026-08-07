# maf_export2.g — session-16 version of maf_export.g: writes the 8
# representative (±1,+1)-cell presentations WITH THE HONEST dirTbBase
# (logged fix, logged fix 2026-07-15: r^-1*M^(-e5)*r*B, sign anti-coupled to e5)
# as GASP/KBMAG rewriting-system files for MAF, into maf_runs2/.
# Run:   gap -q -A -T maf_export2.g
# Then on a machine with MAF:  MAFDIR=maf_runs2 ./maf_certify.sh   (or copy
# maf_certify.sh's loop with cd maf_runs2). Expected: same certificates as
# kbmag (kb_certify_out2.txt G1 / kb_diag2_full_out.txt).
LoadPackage("kbmag");
F := FreeGroup("x","y","r","s","A","B","M","N");;
x:=F.1;;y:=F.2;;r:=F.3;;s:=F.4;;A:=F.5;;B:=F.6;;M:=F.7;;N:=F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
base := [comm(x,y)*comm(r,s), A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
         B*x*B^-1*y, B*r*B^-1*r^-1,
         B*(s^-1*r^-1*y*x)*B^-1*(r^-1*s^-1*x)^-1];;   # incl. R3
mkG := function(m, n, e3, e4, e5, eA, eB)
  return F / Concatenation(base,
    [ A*s*A^-1*(N^e3*y)^-1, B*y*B^-1*(M^e4*y*x)^-1,
      B*s*B^-1*(r^-1*M^e5*r*s)^-1,
      M*((A*r^-1)*((r*x)^-1)^n)^eA,
      N*((r^-1*M^(-e5)*r*B)*(s*r^-1*s^-1)^m)^eB ]);
end;;
CASES := [ [ [1,1],  [ 1, 1, 1, 1, 1] ], [ [1,1],  [-1,-1,-1,-1,-1] ],
           [ [1,1],  [ 1,-1, 1,-1, 1] ], [ [1,1],  [-1, 1,-1, 1,-1] ],
           [ [-1,1], [ 1, 1, 1, 1, 1] ], [ [-1,1], [-1,-1,-1,-1,-1] ],
           [ [-1,1], [ 1,-1, 1,-1, 1] ], [ [-1,1], [-1, 1,-1, 1,-1] ] ];;
if not IsDirectoryPath("maf_runs2") then Exec("mkdir maf_runs2"); fi;
for i in [1..8] do
  cs := CASES[i];;
  G := mkG(cs[1][1], cs[1][2], cs[2][1], cs[2][2], cs[2][3], cs[2][4], cs[2][5]);;
  H := Image(IsomorphismSimplifiedFpGroup(G));;
  rws := KBMAGRewritingSystem(H);;
  WriteRWS(rws, Concatenation("maf_runs2/case", String(i)));
  Print("wrote case ", i, " (gens=", Length(GeneratorsOfGroup(H)), ")\n");
od;
FS := FreeGroup("a","b","c","d");;
SG := FS / [ FS.1*FS.2*FS.1^-1*FS.2^-1*FS.3*FS.4*FS.3^-1*FS.4^-1 ];;
WriteRWS(KBMAGRewritingSystem(SG), "maf_runs2/surface");
Print("wrote surface control\n");
QUIT_GAP(0);
