# kb_diag2_full.g — Knuth–Bendix certification of the FULL diagram-2 grid
# (2026-07-15; logged extension).
#
# Why: with the honest dirTbBase (logged fix, 2026-07-15) diagram 2 becomes
# almost entirely enumeration-blind (vdiag2_out2.txt: 8/576 enum-trivial,
# 568 blowup, H1 = 0 everywhere, zero finite>1). The V-DIAG cross-check
# therefore needs rewriting-system certificates across the whole grid, not
# only the two former blowup cells (kb_certify.g G2). Same certificate
# standard: confluence + language size 1 + every generator reduces to the
# identity. Diagram-2 relators only (no R3): a triviality verdict needs no
# completeness.
#
# R3 note (diag_g2_probe.g): without a completion relation the honest
# diagram-2 hard cells are KB-INCONCLUSIVE (no confluence; kb_certify G2 (run record)
# = 0/128 "not certified", never an adverse verdict). R3 = B k3 B^-1 = psi(k3)
# is a TRUE relation of pi_1(V'), derived from bundle transport (paper: completeness section;
# engine D5) — diagram-INDEPENDENT — and triviality certificates need only
# true relations. It is therefore included here; with it every probed case
# certifies (probe log committed).
#
# Expected: 576/576 certified trivial. Run record: kb_diag2_full_out.txt.

LoadPackage("kbmag");

F := FreeGroup("x","y","r","s","A","B","M","N");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;;
A := F.5;; B := F.6;; M := F.7;; N := F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
base := [comm(x,y)*comm(r,s), A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
         B*x*B^-1*y, B*r*B^-1*r^-1,
         B*(s^-1*r^-1*y*x)*B^-1*(r^-1*s^-1*x)^-1];;  # R3 (true, diagram-independent)
dirTaBase := A*r^-1;;  dirTaFib := (r*x)^-1;;
# honest dirTbBase (logged fix, 2026-07-15): sign anti-coupled to e5
dirTbBase := function(e5) return r^-1*M^(-e5)*r * B; end;;
dirTbFib := s*r^-1*s^-1;;

certify := function(G)
  local H, rws, kb;
  H := Image(IsomorphismSimplifiedFpGroup(G));
  rws := KBMAGRewritingSystem(H);
  kb := CALL_WITH_CATCH(KnuthBendix, [rws]);
  if kb[1] <> true or not IsConfluent(rws) then return "inconclusive"; fi;
  if Size(rws) = 1 and
     ForAll(GeneratorsOfGroup(FreeStructureOfRewritingSystem(rws)),
            z -> IsOne(ReducedForm(rws, z))) then
    return "trivial";
  fi;
  return "nontrivial-or-mixed";
end;;

tot := 0;; nn := 0;;
for mn in Cartesian([-1,0,1],[-1,0,1]) do
  ct := 0;
  for eA in [1,-1] do for eB in [1,-1] do
  for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do for e6 in [1,-1] do
    nn := nn + 1;
    G := F / Concatenation(base,
      [ A*s*A^-1*(N^e3*y)^-1,
        B*y*B^-1*(M^e4*y*x)^-1,
        B*s*B^-1*(N^e6 * r^-1*M^e5*r * s)^-1,
        M*(dirTaBase*dirTaFib^mn[2])^eA,
        N*(dirTbBase(e5)*dirTbFib^mn[1])^eB ]);
    if certify(G) = "trivial" then
      tot := tot + 1; ct := ct + 1;
    else
      Print("DIAG2 NOT CERTIFIED: cell(", mn[1], ",", mn[2], ") e=",
            [e3,e4,e5,e6,eA,eB], "\n");
    fi;
  od; od; od; od; od; od;
  Print("DIAG2-KB cell(", mn[1], ",", mn[2], "): certified trivial ", ct, "/64\n");
od;
Print("DIAG2-KB TOTAL: ", tot, "/", nn, " certified trivial\n");
QUIT_GAP(0);
