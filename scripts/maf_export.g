# maf_export.g — write the 8 representative (±1,+1)-cell presentations (and
# the surface-group control) as GASP/KBMAG rewriting-system files for MAF.
LoadPackage("kbmag");
if IsReadableFile("phase2_common.g") then Read("phase2_common.g");
else Read("scripts/phase2_common.g"); fi;
CASES := [ [ [1,1],  [ 1, 1, 1, 1, 1] ], [ [1,1],  [-1,-1,-1,-1,-1] ],
           [ [1,1],  [ 1,-1, 1,-1, 1] ], [ [1,1],  [-1, 1,-1, 1,-1] ],
           [ [-1,1], [ 1, 1, 1, 1, 1] ], [ [-1,1], [-1,-1,-1,-1,-1] ],
           [ [-1,1], [ 1,-1, 1,-1, 1] ], [ [-1,1], [-1, 1,-1, 1,-1] ] ];;
for i in [1..8] do
  cs := CASES[i];;
  G := mkG(cs[1][1], cs[1][2], cs[2][1], cs[2][2], cs[2][3], cs[2][4], cs[2][5]);;
  H := Image(IsomorphismSimplifiedFpGroup(G));;
  rws := KBMAGRewritingSystem(H);;
  WriteRWS(rws, Concatenation("maf_runs/case", String(i)));
  Print("wrote case ", i, "\n");
od;
FS := FreeGroup("a","b","c","d");;
SG := FS / [ FS.1*FS.2*FS.1^-1*FS.2^-1*FS.3*FS.4*FS.3^-1*FS.4^-1 ];;
WriteRWS(KBMAGRewritingSystem(SG), "maf_runs/surface");
Print("wrote surface control\n");
QUIT_GAP(0);
