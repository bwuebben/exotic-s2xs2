# Candidate presentations for pi_1 of Luttinger-surgered LP piece V and variants.
# Fiber F = Sigma_2, pi_1(F) = <x,y,r,s | [x,y][r,s]> (topologist commutator).
# pi_1(R) = pi_1(F) x| F_2<A,B>:  A acts by swap (x<->r, y<->s),  B acts by h*id,
# h: x -> y^-1, y -> yx  (trefoil monodromy, verified).
# Surgery relations: mu_{T_beta} ~ [z, alpha-word],  mu_{T_alpha} ~ [y, beta-word],
# with z in {ys, sy};  filling: mu * (direction) = 1.

F := FreeGroup("x","y","r","s","A","B");;
x := F.1;; y := F.2;; r := F.3;; s := F.4;; A := F.5;; B := F.6;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;

baserels := [
  comm(x,y)*comm(r,s),          # fiber relation
  A*x*A^-1*r^-1,                # alpha-monodromy = swap
  A*y*A^-1*s^-1,
  A*r*A^-1*x^-1,
  A*s*A^-1*y^-1,
  B*x*B^-1*y,                   # beta-monodromy = h * id :  BxB^-1 = y^-1
  B*y*B^-1*x^-1*y^-1,           # ByB^-1 = yx
  B*r*B^-1*r^-1,
  B*s*B^-1*s^-1
];;

tryGroup := function(name, rels)
  local G, ab, tab, sz, q, targets, t, found;
  G := F / rels;
  ab := AbelianInvariants(G);
  if Length(ab) > 0 then
    Print(name, " | H1 = ", ab, "  (NONTRIVIAL via H1)\n");
    return;
  fi;
  tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
           RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
  if tab <> fail then
    sz := Length(tab[1]);
    Print(name, " | H1 = 0 | |G| = ", sz, "\n");
    return;
  fi;
  # enumeration blew past limit: certify nontriviality via perfect quotients
  targets := [ [AlternatingGroup(5),"A5"], [PSL(2,7),"PSL(2,7)"],
               [AlternatingGroup(6),"A6"], [PSL(2,8),"PSL(2,8)"],
               [PSL(2,11),"PSL(2,11)"], [PSL(2,13),"PSL(2,13)"],
               [SL(2,5),"SL(2,5)"] ];
  found := false;
  for t in targets do
    q := GQuotients(G, t[1]);
    if Length(q) > 0 then
      Print(name, " | H1 = 0 | enum>4e5 | NONTRIVIAL: onto ", t[2], "\n");
      found := true; break;
    fi;
  od;
  if not found then
    Print(name, " | H1 = 0 | enum>4e5 | INCONCLUSIVE (no small perfect quotient)\n");
  fi;
end;;

zwords := [ [y*s, "ys"], [s*y, "sy"] ];;
eps := [1,-1];;

# ---------- LP original V ----------
Print("=== LP original V:  B = [z,A]^-e1,  A = [y,B]^-e2 ===\n");
for zw in zwords do for e1 in eps do for e2 in eps do
  tryGroup(Concatenation("LP z=",zw[2]," e1=",String(e1)," e2=",String(e2)),
    Concatenation(baserels,
      [ comm(zw[1],A)^e1 * B,  comm(y,B)^e2 * A ]));
od; od; od;

# ---------- Scheme A: mixed directions  beta' e^m, alpha' c^n ----------
Print("=== Scheme A: B r^m = [z,A]^-e1,  A (xr)^n = [y,B]^-e2 ===\n");
for zw in zwords do for e1 in eps do for e2 in eps do
  for m in [-1,1] do for n in [-1,1] do
    tryGroup(Concatenation("A z=",zw[2]," e1=",String(e1)," e2=",String(e2),
                           " m=",String(m)," n=",String(n)),
      Concatenation(baserels,
        [ comm(zw[1],A)^e1 * B * r^m,  comm(y,B)^e2 * A * (x*r)^n ]));
  od; od;
od; od; od;

# ---------- Scheme B: LP surgeries + two parallel-torus surgeries killing e, c ----------
Print("=== Scheme B: LP rels + r = [z,A]^-e3, xr = [y,B]^-e4 ===\n");
for zw in zwords do for e1 in eps do for e2 in eps do
  for e3 in eps do for e4 in eps do
    tryGroup(Concatenation("B z=",zw[2]," e=(",String(e1),",",String(e2),",",
                           String(e3),",",String(e4),")"),
      Concatenation(baserels,
        [ comm(zw[1],A)^e1 * B,      comm(y,B)^e2 * A,
          comm(zw[1],A)^e3 * r,      comm(y,B)^e4 * (x*r) ]));
  od; od;
od; od; od;

QUIT_GAP(0);
