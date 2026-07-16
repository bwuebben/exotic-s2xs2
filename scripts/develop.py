#!/usr/bin/env python3
"""
develop.py — combinatorial developing-map engine for the octagon model.

Words over {x,y,r,s}; capital = inverse. The genus-2 octagon has edge word
    E1..E8 = x y x^-1 y^-1 r s r^-1 s^-1   (relator K = "xyXYrsRS"),
all vertices identified to the basepoint p. Corner V_i sits between E_{i-1} and E_i;
corner prefixes w_i = the relator prefix ending at V_i.

Universal-cover bookkeeping:
  * the tile with its corner V_i AT the base lift p~ is w_i^-1 . D0;
  * crossing OUT of a tile through its edge E_i multiplies the deck element on the
    right by g_out(i) = w_i * inv(w_{partner(i)+1})   (E1~E3, E2~E4, E5~E7, E6~E8,
    identified with orientation reversal; w_9 = the full relator);
  * a based loop certificate = (depart corner i, [crossed edges in order], arrive
    corner j)  ==>  based word = inv(w_i) * prod(g_out(k)) * w_j.

Group arithmetic: free reduction + Dehn's algorithm for the surface relator
(valid: genus-2 surface groups have Dehn presentations); equality u == v tested by
Dehn-reducing u*inv(v) to the empty word.  Abelianization as a second check.
"""

GENS = "xyrs"
REL  = "xyXYrsRS"          # [x,y][r,s]

def inv(w):  return w[::-1].swapcase()

def free_reduce(w):
    out = []
    for ch in w:
        if out and out[-1] == ch.swapcase(): out.pop()
        else: out.append(ch)
    return "".join(out)

def _rotations(w):
    return {w[i:] + w[:i] for i in range(len(w))}

_DEHN_TABLE = None
def _dehn_table():
    # subwords of length >= 5 of cyclic rotations of REL and inv(REL),
    # mapped to their shorter complements
    global _DEHN_TABLE
    if _DEHN_TABLE is None:
        t = {}
        for R in _rotations(REL) | _rotations(inv(REL)):
            for L in range(5, 9):                 # more than half of 8
                sub, comp = R[:L], inv(R[L:])
                t[sub] = comp                     # sub == comp in the group
        _DEHN_TABLE = t
    return _DEHN_TABLE

def dehn_reduce(w):
    w = free_reduce(w)
    t = _dehn_table()
    changed = True
    while changed:
        changed = False
        for L in (8, 7, 6, 5):
            for i in range(0, len(w) - L + 1):
                sub = w[i:i+L]
                if sub in t:
                    w = free_reduce(w[:i] + t[sub] + w[i+L:])
                    changed = True
                    break
            if changed: break
    return w

def eq(u, v):   return dehn_reduce(u + inv(v)) == ""
def is1(u):     return dehn_reduce(u) == ""

def abel(w):
    from collections import Counter
    c = Counter()
    for ch in w: c[ch.lower()] += (1 if ch.islower() else -1)
    return tuple(c[g] for g in GENS)

# ---------------- octagon combinatorics ----------------
W = ["", "", "x", "xy", "xyX", "xyXY", "xyXYr", "xyXYrs", "xyXYrsR", REL]  # W[i]=w_i, W[9]=relator
PARTNER = {1: 3, 3: 1, 2: 4, 4: 2, 5: 7, 7: 5, 6: 8, 8: 6}

def g_out(i):
    return free_reduce(W[i] + inv(W[PARTNER[i] + 1]))

def based_word(depart, crossings, arrive):
    w = inv(W[depart])
    for k in crossings: w += g_out(k)
    w += W[arrive]
    return dehn_reduce(w)

def swap_letters(w):
    m = {"x":"r","y":"s","r":"x","s":"y","X":"R","Y":"S","R":"X","S":"Y"}
    return "".join(m[c] for c in w)

def report(label, val, expect=None):
    ok = "" if expect is None else ("   OK" if eq(val, expect) else f"   *** MISMATCH (expected {expect})")
    print(f"{label:34s} = {val or '1'}{ok}")
    if expect is not None and not eq(val, expect):
        raise SystemExit(f"VALIDATION FAILED: {label}")

print("=== crossing decks g_out(i) (Dehn-reduced) ===")
for i in range(1, 9):
    print(f"  g_out(E{i}) = {dehn_reduce(g_out(i)) or '1'}")

print("\n=== V1: edge-loop parallels reproduce the generators ===")
# a loop from p through wedge V_i, across the tile parallel to edge E_i, out at wedge V_{i+1}
report("parallel E1 (V1 -> V2)", based_word(1, [], 2), "x")
report("parallel E2 (V2 -> V3)", based_word(2, [], 3), "y")
report("parallel E3 (V3 -> V4)", based_word(3, [], 4), "X")
report("parallel E5 (V5 -> V6)", based_word(5, [], 6), "r")
report("parallel E6 (V6 -> V7)", based_word(6, [], 7), "s")
report("parallel E8 (V8 -> V1)", based_word(8, [], 1), "S")

print("\n=== V2: the vertex-link loop is contractible ===")
# link cyclic order of wedges: V1 V4 V3 V2 V5 V8 V7 V6, crossing OUT through
# E1, E4, E3, E2, E5, E8, E7, E6 in turn (the octagon link walk)
link = ""
for k in [1, 4, 3, 2, 5, 8, 7, 6]: link += g_out(k)
report("product around the vertex link", dehn_reduce(link), "")

print("\n=== V3: push-off curve certificates ===")
# a = parallel of the x-circle on the {V4,V3} side: strip along E3 + corner arc
#     V4 -> V3 crossing the ray E4-start (= crossing OUT through E4); whisker V3.
report("a  (whisker V3)", based_word(3, [4], 3), "X")
# b = parallel of the y-circle on the {V3,V2} side: strip along E2 + corner arc
#     V3 -> V2 crossing OUT through E3; whisker V2.
report("b  (whisker V2)", based_word(2, [3], 2), "y")
# e = rho(a): strip along E7 + corner arc V8 -> V7 crossing OUT through E8; whisker V7.
report("e  (whisker V7)", based_word(7, [8], 7), "R")
# d = rho(b): strip along E6 + corner arc V7 -> V6 crossing OUT through E7; whisker V6.
report("d  (whisker V6)", based_word(6, [7], 6), "s")
# equivariance: swap(word of a at V3) must equal word of e at V7, etc.
report("swap(a@V3) == e@V7", swap_letters(based_word(3, [4], 3)), based_word(7, [8], 7))
report("swap(b@V2) == d@V6", swap_letters(based_word(2, [3], 2)), based_word(6, [7], 6))

print("\n=== D1: the invariant curve c ===")
# c = arc1 (E2 -> E8) u rho(arc1) (E6 -> E4); traversal: cross OUT through E8,
# then OUT through E4; rho-symmetric certificate [8, 4].
for wedge in (2, 3, 5, 1):
    print(f"  c (whisker V{wedge}) = {based_word(wedge, [8, 4], wedge) or '1'}"
          f"   abel={abel(based_word(wedge, [8, 4], wedge))}")
report("c free class check (whisker V2)", based_word(2, [8, 4], 2), "XR")  # = (rx)^-1
# reversed traversal:
report("c reversed (whisker V2)", based_word(2, [2, 6], 2), "rx")
# rho-equivariance of the invariant curve: swap of word at V2 vs word at V6
report("swap(c@V2) == c@V6", swap_letters(based_word(2, [8, 4], 2)), based_word(6, [4, 8], 6))

print("\n=== D2: the invariant curve z — certificate search ===")
# z = delta u rho(delta), delta from the x-pair to the r-pair; the closed-up
# crossing certificates are [b, a+4-partner] patterns; enumerate all rho-symmetric
# 2-crossing certificates on {E1,E3} x {E5,E7} in both traversal directions and
# report free classes (whisker V2 for comparability).
from itertools import product as iproduct
print("  candidates (whisker V2):")
seen = {}
for aedge, bedge in iproduct((1, 3), (5, 7)):
    # arc1 from E_a to E_b: consistency needs partner(bedge) == aedge+4
    if PARTNER[bedge] != aedge + 4:  continue
    for cert in ([bedge, aedge + 4 - 0], ):  # traversal 1: out E_b then out rho(E_a)...
        pass
# systematic: certificate = [out E_b, out E_{a'}] where entering at E_a means the
# previous crossing exited the partner of E_a; closure requires the second exit's
# partner to be E_a: second crossing edge = PARTNER-inverse... enumerate directly:
for c1, c2 in iproduct(range(1, 9), range(1, 9)):
    # rho-symmetry: c2 == ((c1+4-1) % 8) + 1 ; restrict to x/r-pair crossings only
    if c2 != ((c1 + 3) % 8) + 1:  continue
    if c1 not in (1, 3, 5, 7):    continue
    wword = based_word(2, [c1, c2], 2)
    print(f"    crossings [E{c1},E{c2}]: {wword or '1'}  abel={abel(wword)}")
    seen[(c1, c2)] = wword

print("\n(select the z certificate matching class ys/sy up to whisker+orientation;")
print(" disjointness from b,d to be certified at the figure level)")

print("\n=== D3: lasso words for the corrections ===")
# Convention P1: edge-parallel path segments are pushed into the tile on the
# LOWER-index edge of each identified pair (E2 over E4, E6 over E8, etc.).
#
# delta = s1 . (c-arc from c_s to c_y) . y1^-1  relates the two T_alpha meridian
# basings (along the partial-s and partial-y whiskers).
# Route via arc2 = rho(arc1) (inside D0 from P_s' in E6 to P_y' in E4):
#   depart wedge V6 (E6-parallel to P_s'), arc2 in-tile, cross OUT through E4
#   (to the E2-side per P1), E2-parallel back to p: arrive wedge V2.
delta = based_word(6, [4], 2)
print(f"  delta  (route arc2)          = {delta or '1'}   abel={abel(delta)}")
# Route via arc1 (cross the s-circle first: out through E6, then arc1^-1 to E2):
deltap = based_word(6, [6], 2)
print(f"  delta' (route arc1)          = {deltap or '1'}   abel={abel(deltap)}")
# V4: the two routes must differ by (a whiskered copy of) the closed curve c:
# delta' . delta^-1 should be conjugate to c^{+-1} (free class xr / (xr)^-1).
u = dehn_reduce(deltap + inv(delta))
print(f"  delta'.delta^-1              = {u or '1'}   abel={abel(u)}")
found = None
CONJTRY = [""] + [g for g in "xyrsXYRS"] + [a+b for a in "xyrsXYRS" for b in "xyrsXYRS"]
for tgt in ("XR", "rx", "RX", "xr"):
    for w in CONJTRY:
        if eq(u, w + tgt + inv(w)):
            found = (tgt, w); break
    if found: break
if found:
    print(f"  V4 OK: delta'.delta^-1 ~ {found[0]} (conjugator '{found[1] or '1'}')")
else:
    raise SystemExit("V4 FAILED: routes do not differ by the curve c")

print("\n=== D4: position tables (normal-curve data, derived by disjointness) ===")
print("""  s-edge (based loop s = E6 from V6 to V7):
     order of marked points from the start:  [ c-crossing P_s' ] < [ e-crossing s_e ]
     (e's corner arc crosses ray E8-start, i.e. E6 near its END; c's endpoint must
      avoid e's strip along E7 and the V8-corner => P_s' in the middle.  c.e = 0
      pins the order.)
  y-edge (based loop y = E2 from V2 to V3):
     order:  [ c-crossing c_y (middle) ] < [ a-crossing (near END, ray E4-start) ]
     (only c matters for the lassos; a is not a torus fiber.)
  Realization checks recorded: c crosses d's strip immediately on leaving E6
  (c.d = 1) and b's strip immediately on leaving E2 (c.b = 1); c avoids a, e, O, p;
  z's E7-crossing in the middle of E7 meets e's strip once (z.e = 1) and avoids
  d's corner arc (z.d = 0).  Side-of-O choice for arc1 is immaterial for words
  (O is not removed in C) and for SS23 (c misses O either way).""")

# ---------------- D4CERT: chord-diagram certificate (machine check) ----------------
# The prose above is now certified combinatorially.  Cut along the four edge
# circles: the octagon is a disk, and each curve of {a, b, d, e, c, z} appears
# as straight chords whose endpoints sit on the edges, tied in partner pairs
# E_i(t) ~ E_partner(i)(1 - t).  Facts used:
#   * two chords in a disk cross (once) iff their endpoints interleave in the
#     boundary cyclic order, and straight chords realize ALL pairwise minimal
#     crossing numbers simultaneously;
#   * hence a family of arcs with prescribed pairwise geometric intersection
#     numbers is realizable iff the interleaving matrix of the endpoint order
#     equals the prescribed matrix -- a finite check.
# The free data are the orders of the two crossing points on each of the four
# edge circles (2^4 = 16 configurations); everything else (which edges each
# arc connects) is the certified crossing structure of D1/D2/V3.  The
# certificate below shows EXACTLY ONE configuration realizes the model's
# intersection table, and in it the D4 orders hold; each order is forced by
# one disjointness (c.e = 0 -> s-edge order; c.a = 0 -> y-edge order;
# z.b = 0 and z.d = 0 -> the x/r-circle orders).
print("\n=== D4CERT: chord-diagram planarity certificate ===")

REQUIRED = {  # curve-level geometric intersection numbers of the frozen model
    ("a","b"): 1, ("a","d"): 0, ("a","e"): 0, ("b","d"): 0, ("b","e"): 0,
    ("d","e"): 1, ("a","c"): 0, ("b","c"): 1, ("c","d"): 1, ("c","e"): 0,
    ("a","z"): 1, ("b","z"): 0, ("d","z"): 0, ("e","z"): 1, ("c","z"): 2,
    ("c","c"): 0, ("z","z"): 0,
}

def _d4_arcs(fy, fx, fr, fs):
    # One binary flag per edge circle: True = the c/z-crossing point precedes
    # the pushoff-curve's point in the parameter of the LOW-index edge.
    def pos(flag): return (0.35, 0.65) if flag else (0.65, 0.35)
    p_y, q_y = pos(fy)   # E2: c_y (gamma) vs a-copy      (E4 copies at 1 - t)
    p_x, q_x = pos(fx)   # E1: z-alpha    vs b-copy       (E3 copies)
    p_r, q_r = pos(fr)   # E5: z-beta     vs d-copy       (E7 copies)
    p_s, q_s = pos(fs)   # E6: P_s' (rho gamma) vs e-copy (E8 copies)
    return {
        ("a", 0): ((2, q_y), (4, 1 - q_y)),
        ("b", 0): ((1, q_x), (3, 1 - q_x)),
        ("d", 0): ((5, q_r), (7, 1 - q_r)),
        ("e", 0): ((6, q_s), (8, 1 - q_s)),
        ("c", 0): ((2, p_y), (8, 1 - p_s)),     # gamma:      c_y  -> P_s'
        ("c", 1): ((6, p_s), (4, 1 - p_y)),     # rho(gamma): P_s' -> c_y
        ("z", 0): ((1, p_x), (7, 1 - p_r)),     # alpha
        ("z", 1): ((5, p_r), (3, 1 - p_x)),     # beta = rho(alpha)
    }

def _d4_matrix(arcs):
    # global boundary coordinate: edge index + position in (0,1)
    def glob(ep): return ep[0] + ep[1]
    from collections import Counter
    tot = Counter()
    keys = sorted(arcs)
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            k1, k2 = keys[i], keys[j]
            if k1 == k2: continue
            a1, a2 = sorted(map(glob, arcs[k1])), sorted(map(glob, arcs[k2]))
            inside = sum(1 for q in a2 if a1[0] < q < a1[1])
            if inside == 1:   # interleaved <=> the straight chords cross once
                pair = tuple(sorted((k1[0], k2[0])))
                tot[pair] += 1
    return {p: tot.get(p, 0) for p in REQUIRED}

_passing = []
for fy in (True, False):
    for fx in (True, False):
        for fr in (True, False):
            for fs in (True, False):
                got = _d4_matrix(_d4_arcs(fy, fx, fr, fs))
                bad = [f"{p[0]}.{p[1]}={got[p]}(need {REQUIRED[p]})"
                       for p in REQUIRED if got[p] != REQUIRED[p]]
                tag = f"(c_y<a:{int(fy)} z<b:{int(fx)} z<d:{int(fr)} P_s'<e:{int(fs)})"
                if bad:
                    print(f"  config {tag}: FAIL  [{'; '.join(bad)}]")
                else:
                    print(f"  config {tag}: REALIZABLE")
                    _passing.append((fy, fx, fr, fs))

if _passing != [(True, True, True, True)]:
    raise SystemExit("D4CERT FAILED: expected exactly the D4 configuration "
                     f"to be realizable, got {_passing}")
print("""  D4CERT OK: exactly one of the 16 configurations is realizable, and it is
  the D4 table:  P_s' precedes the e-crossing on the s-edge (forced by c.e=0),
  c_y precedes the a-crossing on the y-edge (forced by c.a=0), and z's
  crossings precede b's resp. d's on the x/r-circles (forced by z.b=z.d=0).
  The position table is a certificate, not a drawing.""")

print("\n=== D5: the third basis element of pi_1(F - nu(c)) ===")
# kappa_3 = tree(p -> Psb) . arc2mid . tree(Pyb -> p): a boundary-collar path of
# the mid region R_mid from wedge V7 to wedge V3, entirely inside D0: certificate
# (depart V7, no crossings, arrive V3).
k3 = based_word(7, [], 3)
print(f"  kappa_3 = {k3}   abel={abel(k3)}")
# c-avoidance check: <kappa_3, [c]> with [c] = -x-r:
# <k, -x-r> = -k_y*(y.x) - k_s*(s.r) = k_y + k_s   (y.x = s.r = -1)
ka = abel(k3)
pair_c = ka[1] + ka[3]
print(f"  <kappa_3, c> = {pair_c}   (must be 0)")
if pair_c != 0: raise SystemExit("D5 FAILED: kappa_3 meets c algebraically")
def psi_letter(ch):
    m = {"x":"Y","y":"yx","r":"r","s":"s","X":"y","Y":"XY","R":"R","S":"S"}
    return m[ch]
def psi_word(w): return dehn_reduce("".join(psi_letter(ch) for ch in w))
print(f"  psi(kappa_3) = {psi_word(k3)}")
print(f"  (R3: B kappa_3 B^-1 = psi(kappa_3) — the T_alpha completion relation)")
