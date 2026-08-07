if IsReadableFile("phase2_common.g") then Read("phase2_common.g");
else Read("scripts/phase2_common.g"); fi;
for c in [ [1,1,[1,1,1,1,1]], [-1,1,[1,1,1,1,1]] ] do
  G := mkG(c[1], c[2], c[3][1], c[3][2], c[3][3], c[3][4], c[3][5]);;
  H := Image(IsomorphismSimplifiedFpGroup(G));;
  Print("deep TC cell(", c[1], ",", c[2], ") e=", c[3], " max=1e8 ", stamp(), "\n");
  tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(H),
           RelatorsOfFpGroup(H), [] : max := 100000000, silent := true);
  if tab = fail then Print("  blowup at 1e8 ", stamp(), "\n");
  else Print("  TERMINATED: |G| = ", Length(tab[1]), " ", stamp(), "\n"); fi;
od;
QUIT_GAP(0);
