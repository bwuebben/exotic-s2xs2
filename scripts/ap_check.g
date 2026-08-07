# Independent verification of Akhmedov-Park arXiv:1005.3346, Theorem 9:
# pi_1(M_n^p) = Z/p from the presentation of Lemma 8.
# Generators: a1,b1,a2,b2,c1,d1,C2 (= c2-tilde), d2.
F := FreeGroup("a1","b1","a2","b2","c1","d1","C2","d2");;
a1:=F.1;; b1:=F.2;; a2:=F.3;; b2:=F.4;; c1:=F.5;; d1:=F.6;; C2:=F.7;; d2:=F.8;;
comm := function(u,v) return u*v*u^-1*v^-1; end;;
APrels := function(n, p) return [
  a2^-1 * C2^-1*a1*C2,
  b2^-1 * C2^-1*b1*C2,
  b1^-1 * C2^-1*b2*C2,
  comm(b2,d2),
  comm(a1^-1*b1^-1*a2, d2),
  comm(a2^-1*b2^-1*a1, d2),
  comm(b1^-1,d1^-1)^n * a1^-1,
  comm(a1^-1,d1) * b1^-1,
  comm(b2^-1,d1^-1) * c1^-1,
  comm(b2,c1^-1) * d1^-1,
  C2^-1*a1*a2*C2*a1^-1*a2^-1 * (d2^p)^-1,
  comm(a1, C2*d2^-1*C2^-1) * (C2*b2)^-1,
  comm(a1,c1), comm(b1,c1), comm(a2,c1), comm(a2,d1), comm(b1,d2),
  comm(a1,b1)*comm(a2,b2),
  comm(c1,d1)*comm(C2,d2)
]; end;;
check := function(n,p)
  local G, ab, tab;
  G := F / APrels(n,p);
  ab := AbelianInvariants(G);
  tab := CosetTableFromGensAndRels(FreeGeneratorsOfFpGroup(G),
           RelatorsOfFpGroup(G), [] : max := 400000, silent := true);
  if tab <> fail then
    Print("AP (n=",n,",p=",p,"): H1=",ab," |pi1|=",Length(tab[1]),"\n");
  else
    Print("AP (n=",n,",p=",p,"): H1=",ab," enum>4e5\n");
  fi;
end;;
check(1,1); check(2,1); check(3,1); check(1,2); check(1,3); check(1,0);
QUIT_GAP(0);
