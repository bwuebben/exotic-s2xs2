# kb_certify.g — Knuth–Bendix certification of the FULL admissible grid
# (session 14, 2026-07-10). FROZEN once run.
#
# Method: for each presented group, run Knuth–Bendix completion (kbmag 1.5.11,
# shortlex) on the Tietze-simplified presentation. A certificate of triviality
# consists of: (a) the completion reaches confluence, (b) the language of
# irreducible words has size 1, and (c) EVERY generator explicitly reduces to
# the identity. Note (c) alone already proves triviality independently of the
# confluence claim: each KB rule is a consequence of the relators by
# construction, so a rewriting derivation g ->* 1 is a proof that g = 1 in the
# group. Trivial verdicts need no completeness of the relation set (the
# presented group surjects onto pi_1), and are conservative under both
# diagrams.
#
# Battery:
#   C1  positive control: genus-2 surface group -> automatic, Size = infinity
#       (the pipeline reports infinity when that is the truth);
#   C2  negative control: partial relation set (clean relations only, an
#       infinite group) -> KB does NOT complete spuriously;
#   G1  diagram 1, completed presentation (adds R3): all 9 cells x 32 signs
#       = 288 cases;
#   G2  diagram 2 (independent Bs-derivation, no R3): the two former blowup
#       cells x 64 signs = 128 cases.
#
# Expected output: 288/288 and 128/128 certified trivial; controls behave.
# 2026-07-15 NOTE (paper: the pushoff-basing correction): with the honest
# dirTbBase (logged fix), G1 stays 288/288, but G2 WITHOUT a completion
# relation becomes KB-inconclusive (0/128 "not certified" = no confluence,
# never an adverse verdict; see diag_g2_probe.g). The diagram-2 certification
# now lives in kb_diag2_full.g (full 9-cell grid, R3 added — true and
# diagram-independent).

LoadPackage("kbmag");
if IsReadableFile("phase2_common.g") then Read("phase2_common.g");
else Read("scripts/phase2_common.g"); fi;

# LOGGED FIX 2026-07-15 (paper: the pushoff-basing correction): honest
# dirTbBase = r^-1*M^(-e5)*r*B (sign anti-coupled to e5). mkG overridden
# LOCALLY — phase2_common.g stays untouched as the quotient-sweep record.
mkG := function(m, n, e3, e4, e5, eA, eB)
  return F / Concatenation(base,
    [ A*s*A^-1*(N^e3*y)^-1, B*y*B^-1*(M^e4*y*x)^-1,
      B*s*B^-1*(r^-1*M^e5*r*s)^-1,
      M*((A*r^-1)*((r*x)^-1)^n)^eA,
      N*((r^-1*M^(-e5)*r*B)*(s*r^-1*s^-1)^m)^eB ]);
end;;

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

Print("GAP ", GAPInfo.Version, " + kbmag ",
      InstalledPackageVersion("kbmag"), "\n");

# --- C1: positive control ---
FS := FreeGroup("a","b","c","d");;
SG := FS / [ Comm(FS.1,FS.2)*Comm(FS.3,FS.4) ];;
rwsS := KBMAGRewritingSystem(SG);;
if CALL_WITH_CATCH(AutomaticStructure, [rwsS, true]) = [true, true] then
  Print("C1 surface group: automatic, Size = ", Size(rwsS),
        "  (must be infinity)\n");
else
  Print("C1 surface group: FAILED to build automatic structure\n");
fi;

# --- C2: negative control (clean relations only; infinite) ---
Gc := F / base;;
rwsc := KBMAGRewritingSystem(Image(IsomorphismSimplifiedFpGroup(Gc)));;
CALL_WITH_CATCH(KnuthBendix, [rwsc]);;
Print("C2 partial-relations group: confluent = ", IsConfluent(rwsc),
      "  (must be false: KB does not complete spuriously here)\n");

# --- G1: diagram 1, completed presentation, full grid ---
t1 := 0;; n1 := 0;;
for m in [-1,0,1] do for n in [-1,0,1] do
  ct := 0;
  for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
  for eA in [1,-1] do for eB in [1,-1] do
    n1 := n1 + 1;
    if certify(mkG(m,n,e3,e4,e5,eA,eB)) = "trivial" then
      t1 := t1 + 1; ct := ct + 1;
    else
      Print("G1 NOT CERTIFIED: cell(",m,",",n,") e=",[e3,e4,e5,eA,eB],"\n");
    fi;
  od; od; od; od; od;
  Print("G1 cell(", m, ",", n, "): certified trivial ", ct, "/32\n");
od; od;
Print("G1 TOTAL: ", t1, "/", n1, " certified trivial\n");

# --- G2: diagram 2, the two former blowup cells ---
# NB: comm(u,v) = u v u^-1 v^-1 (phase2_common.g convention, the one in which
# the correction words were derived) — NOT GAP's built-in Comm = u^-1 v^-1 u v.
# A first run of this script used Comm here by mistake: it presents a DIFFERENT
# group (the derived corrections live in the fixed convention), and KB duly
# refused to certify anything (0/128). Kept as a cautionary note; the sweep
# below is in the correct convention.
base2 := [comm(x,y)*comm(r,s), A*x*A^-1*r^-1, A*y*A^-1*s^-1, A*r*A^-1*x^-1,
          B*x*B^-1*y, B*r*B^-1*r^-1];;
dTaB := A*r^-1;; dTaF := (r*x)^-1;; dTbF := s*r^-1*s^-1;;  # dTbB inlined below (2026-07-15 fix)
t2 := 0;; n2 := 0;;
for m in [1,-1] do
  for e3 in [1,-1] do for e4 in [1,-1] do for e5 in [1,-1] do
  for e6 in [1,-1] do for eA in [1,-1] do for eB in [1,-1] do
    n2 := n2 + 1;
    G2 := F / Concatenation(base2,
      [ A*s*A^-1*(N^e3*y)^-1, B*y*B^-1*(M^e4*y*x)^-1,
        B*s*B^-1*(N^e6 * r^-1*M^e5*r * s)^-1,
        M*(dTaB*dTaF^1)^eA,
        # LOGGED FIX 2026-07-15: honest dirTbBase (dTbB slot), anti-coupled sign
        N*((r^-1*M^(-e5)*r*B)*dTbF^m)^eB ]);
    if certify(G2) = "trivial" then
      t2 := t2 + 1;
    else
      Print("G2 NOT CERTIFIED: m=", m, " e=", [e3,e4,e5,e6,eA,eB], "\n");
    fi;
  od; od; od; od; od; od;
od;
Print("G2 TOTAL: ", t2, "/", n2, " certified trivial\n");
QUIT_GAP(0);
