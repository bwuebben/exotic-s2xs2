# exotic-s2xs2

Scripts and run logs for

- **An exotic S²×S² and an exotic ℂP²#ℂP̄²**
  (`papers/exotic-s2xs2-and-cp2.pdf`), 32 pp: Part I the reduction
  (§§3–6), Part II the computation (§§7–12), machine verification in
  Appendix A, known-answer calibrations in Appendix B,

by Bernd Johannes Wuebben.

A step-by-step expository companion to the computation is
`papers/walkthrough.pdf` (11 pp). It follows the August 16 basing repair and
proof boundary: the marked-surface bridge, every corrected relation and direction word, the
drilled-fiber relation, and the epimorphism onto $\pi_1(V)$ are worked at
tutorial grain; the specified $y_1/Ax$ $4{,}096$-case decision is separated
from the adjacent $y_2/Ar^{-1}$ run and from the
nine-pair diagnostic grids; and the Baldridge–Kirk T⁴ configuration is a fully
worked finite-fingerprint calibration example. Table 1 of the manuscript
(§9.1) lists every
relation of the system used at the specified Lidman–Piccirillo surgery with
its verification pointer, and Appendix A.4 gives the protocol for diffing an
independently derived relation sheet against it. Everything needed to reproduce every number in the
manuscript is here; total runtime is minutes on a laptop (except the optional
finite-quotient sweep, ~56 CPU-hours).

> **Status.** The manuscript is a draft and has not been submitted. Earlier
> versions were reviewed in correspondence with Lidman and Piccirillo, whose
> own computation of π₁(V) disagrees; the disagreement is unresolved, and the
> verification guide (Appendix A.3–A.4) is written so that any reader can
> localize a concrete error to one crossing count, one conjugating path or one
> framing word.

## The result

Lidman–Piccirillo (arXiv:2505.14387) built spin rational homology 4-spheres $B$
and $W$ with the same integer cohomology ring, distinguished by whether the
figure-eight knot is slice. **Part I** proves that making that pair
*homeomorphic* is equivalent to the exotic $S^2\times S^2$ problem: for any
admissible surgery variant $V'$ of their key piece, $W'\cong_{\rm homeo}B$ iff
$\pi_1(V'\cup_\sigma V')=1$, in which case the double is an exotic
$S^2\times S^2$. **Part II** then computes, from a fully explicit based model,
the fundamental group of the piece itself:

> **π₁(V) = 1 for the Lidman–Piccirillo piece $V=V'_{0,0}$** — the result of
> their two $+1$ Luttinger surgeries when the two base-direction curves are
> the Lagrangian push offs derived in §§8.4–8.5.

Lidman and Piccirillo leave the parametrizing curves on the two surgery tori
arbitrary, so $V$ names a manifold only after those choices are made; the
manuscript fixes them and makes no π₁ claim for any other parametrization.

The geometric identifications are proved in the manuscript, not certified by
this package: the equivariant marked-surface normal form for the octagon
(Lemma 7.1), the transport-annulus arguments behind the three corrections
(§8.3), the framing identification (§8.7), and the drilled-fiber bases
(Lemma 10.1). What the scripts reproduce is the algebra: the development of
each based word from its stated combinatorial description, the run summaries,
the consistency tests, and the finite-presentation decisions. The relation
system is proved to **surject** onto π₁(V) (§10), so a trivial outcome is a
proof and a nontrivial one would concern the presented group alone.

The proof-level decision is `fixed_v_certify.g`: for the specified piece with
the $y_1$-based direction word $Ax$, coset enumeration terminates with order
one in **all 4,096** cases of the robustness
family, which independently varies the five signs of the relation sheet, the
left/right placement of all three transport corrections, both conjugating arcs
for the $Bs$ correction, and both arc/sign/placement choices in the $T_\beta$
push off (Remark 11.1). The script separately repeats the $4{,}096$ cases with
the adjacent $y_2$-based word $Ar^{-1}=Ax(rx)^{-1}$; that run is also entirely
trivial but is not used for the specified $n=0$ piece. It uses core GAP only.

The remaining grids — the nine formal exponent pairs $(m,n)\in\{-1,0,1\}^2$ in
two based diagrams, 1,408 relation systems in all — are **diagnostics**. Every
one of them has $H_1=0$ and no terminating enumeration anywhere produced a
finite nontrivial group; the enumeration overflows recorded below carry no
group-theoretic conclusion, and no simple connectivity is claimed for the
eight additional parametrizations.

The derivation is calibrated on the Baldridge–Kirk double-Luttinger tori in T⁴
(`decide_t4.g`). With correctly derived lasso data, all 64 conventions have
the expected abelianization, an explicit non-abelian finite quotient, and the
expected low-index subgroup fingerprint through index 5, including the exact
tr 3 / tr 1 split. Every wrong lasso pair fails this comparison in all 64
cases; half of six wrong families are abelian, the signature of a missed
whisker crossing. These finite fingerprints are discriminating checks, not
proofs that two finitely presented groups are isomorphic. The calibration
found (and led to the repair of)
exactly one such error in an earlier version of one direction word
(`dir_Tβ^base`; §8.5 of the manuscript documents the found-and-fixed
derivation), after which **every conclusion was re-established with the
corrected, sign-coupled word**.

A second calibration targets the one input T⁴ cannot see, the
monodromy-drift word: torus surgery along the drift section of
MT(φ₀)×S¹, compared with independently presented torus-twist graph
manifolds (`decide_seifert.g`). With the basing-coherent words, the
abelianizations and low-index subgroup fingerprints through index 6 match at
three surgery coefficients in both mirror packages, while seven wrong-word
families (including the
whisker-conjugate family no earlier instrument could see) fail visibly.
This test's own first run caught (and fixed) a basing incoherence in its
construction (`confirm_coherent.g` records the localization); the arc
convention it exposed shifts one sweep index and is covered by the
both-arc grid above.

Together: an exotic $S^2\times S^2$, the first pair of homeomorphic closed
4-manifolds distinguished by unconstrained knot slicing, and a simply connected
exotic $\mathbb{CP}^2\sharp \overline{\mathbb{CP}}^2$.

## Requirements

- **GAP 4** (developed on 4.16.0; core library only, except that
  `kb_certify.g` and `kb_diag2_full.g` additionally need the **kbmag**
  package, built from https://github.com/gap-packages/kbmag into GAP's `pkg/`
  directory). GAP is not in Homebrew or MacPorts; see `docs/INSTALL_GAP.md`
  for a minimal source build on macOS, or use your distribution's package on
  Linux (`apt install gap`).
- **Python 3** (standard library only) for the development algorithm.

Run any script with

```
gap -q -A scripts/<name>.g
python3 scripts/develop.py
```

The proof-level check is `gap -q -A scripts/fixed_v_certify.g`; it needs core
GAP only and must report `TOTAL=4096 TRIVIAL=4096` with no overflows, finite
nontrivial outcomes, or nonzero abelianizations first for `FIXED V (y1/Ax)`
and then for `ADJACENT T_ALPHA SECTION (y2/Ar^-1)`. The first run is the one
used by the theorem. The Knuth–Bendix scripts concern the diagnostic grids.

## What each script certifies

| script | certifies | expected output | time |
|---|---|---|---|
| `monodromy_check2.g` | the trefoil monodromy model (fixes [x,y] exactly; h⁶ = boundary twist; h³ = −id) | `true / true / false` + h³ images | <1 s |
| `model_check3.g` | h = T_a∘T_b (order pinned); relation preservation; ψφψ⁻¹φ = h∗h⁻¹; fiber-H₁ Smith form | all `true`, diagonal `[1,1,1,1]` | <1 s |
| `develop.py` | the development algorithm + validations V1–V5; derived words: a, b, d, e; c ≃ (rx)⁻¹; source homology check $[z]=[y]-[s]$; lasso δ = r⁻¹; κ₃ = s⁻¹r⁻¹yx | all `OK`; see `logs/develop_out.txt` | <1 s |
| `pi1_grid.g` | experiment E1: all 72 uncorrected candidate presentations collapse | 72 × `\|G\| = 1` | ~1 s |
| `pi1_v2b.g` | controls + E2 sensitivity (120 seeded trials): both nondegenerate outcomes are generic, so a trivial outcome obtained from guessed words is worthless | `SUMMARY C: TRIVIAL=69 FINITE=0 H1=7 BLOWUP=44` | ~5 s |
| `ap_check.g` | E3: Akhmedov–Park Lemma 8 ⇒ π₁(Mₙᵖ) = ℤ/p at five parameter values | ℤ/1, ℤ/1, ℤ/1, ℤ/2, ℤ/3 | ~6 min |
| **`fixed_v_certify.g`** | **the proof-level decision**: the specified $y_1/Ax$ piece V = V′₀,₀ with R₃ over all 4,096 sign, arc, and placement choices, followed by the adjacent $y_2/Ar^{-1}$ section over the same 4,096 choices; core GAP only | both lines report `TOTAL=4096 TRIVIAL=4096 OVERFLOW=0 FINITE>1=0 H1nonzero=0` (`logs/fixed_v_certify_out.txt`) | ~20 min |
| `decide_t4.g` | **finite-fingerprint calibration** on the Baldridge–Kirk T⁴ configuration: derived lassos (1, b) pass 64/64 in abelianization, explicit non-abelian quotients, and low-index subgroup invariants through index 5, including the tr 3/tr 1 split; every wrong lasso pair fails; BK's own words agree case by case | `logs/t4_out.txt` | ~20 min |
| `gt1_diff.g` | string-level word diff against Baldridge–Kirk Theorem 1 in the punctured model: meridians match up to explicit basing conjugators, push offs match literally including the ℓ₂ = bab⁻¹ chirality (derived, not chosen), the mirror chirality is certified to be a puncture word, and the two relation lists coincide as sets | `logs/gt1_out.txt` | <1 s |
| `decide_seifert.g` | **Seifert drift calibration**: surgery along the drift section λ×S¹ ⊂ MT(φ₀)×S¹ compared with independently presented torus-twist groups; coherent words match abelianization and low-index subgroup fingerprints through index 6 at k ∈ {1, −1, 2}; seven wrong-word families are separated at every k | `logs/seifert_out.txt` | ~3 min |
| `confirm_coherent.g` | localization follow-up to the above: the basing-coherent drift words pass (predictions C1–C2), the incoherent mix fails (C3); the arc goes with the basing | `logs/confirm_out.txt` | ~1 min |
| `cy_pairing.g` (+ `cy_supplement.g`) | the **degenerate-basing pairing derived**: split-generator model at the basing point itself; the whisker's approach side determines transport conjugator and push-off label (y₁-side ⟷ r⁻¹ ⟷ Ax; y₂-side ⟷ x ⟷ Ar⁻¹), cross-pairings fail, and no undecorated transport survives the finite-fingerprint comparison | `logs/cy_out.txt`, `logs/cy_supp_out.txt` | ~2 min |
| `decide.g` | phase-1 harness (self-test mode without inputs; production 288-case sweep of the initial relation system with the derived, sign-coupled direction words) | see `logs/honest_run_2.log` (earlier record: `honest_run_1.log`) | ~1 min |
| `vsens.g` | one-at-a-time perturbation of the derived input words at the surgery variant (drift dropped, push off inverted, pre-repair word restored, correction sign flipped): 14 of 15 perturbations still present the trivial group, H₁ = 0 throughout — the quantitative reason the decision procedures alone settle nothing about geometric correctness | `logs/vsens_out.txt` | ~1 min |
| `vperiph.g`, `vperiph2.g` | peripheral-class probes of the relation group (quotient battery, then low-index coset actions); diagnostic only, since the relation group is known to surject onto the complement group, not to present it | `logs/vperiph_out.txt`, `logs/vperiph2_out.txt` | ~2 min |
| `diag_dirTbBase.g`, `diag_kb.g`, `diag_fullgrid.g`, `diag_r3.g` | the found-and-fixed audit trail: all 8 candidate resolutions of the push-off basing correction at the LP cell (enumeration + Knuth–Bendix certification of every overflow, 256/256 trivial) and both lasso arcs across the full grid (576/576 certified trivial) | `logs/diag_*_out.txt` | ~30 min total |
| `placement_check.g` | correction-placement robustness at the LP cell with R₃ (2³ placements × signs) | `PLACEMENT SWEEP WITH R3 (m=0,n=0): ENUM_TRIVIAL=256 KB_TRIVIAL=0 INCONCLUSIVE=0 FINITE>1=0 H1nonzero=0` | ~1 min |
| `vdiag2.g` | the independently derived second based diagram (576 diagnostic systems; enumeration only, since the corrected word is enumeration-hostile here; H₁ = 0 in all, zero finite>1) | see `logs/vdiag2_out2.txt`; rewriting in `kb_diag2_full.g` | ~2 min |
| `vr_check.g` | known-answer probes (no fillings → ℤ²; fiber-only → ℤ²; one filling → ℤ) | `[0,0] / [0,0] / [0]` | <1 s |
| `decide2.g` | the R₃-augmented relation system: the current $y_1/Ax$ formal sweep | `logs/decide2_out2.txt`: 6 of the 9 formal exponent pairs trivial in all 32 conventions, the LP pair (0,0) among them; the n=−1 column overflows (decided by rewriting below) | ~2 min |
| `phase2_parallel.sh` (+ `phase2_common.g`, `phase2_worker.g`) | historical pre-reindex finite-quotient diagnostics on 8 representative sign cases; these retain the old $Ar^{-1}$ coordinate, with $n_{\rm current}=n_{\rm old}+1$, so the archived logs remain reproducible | `logs/phase2_par_out.txt` | ~6 min (8-wide) |
| `phase3_resume.sh` (+ `phase3_worker.g`; `phase3_parallel.sh` from scratch) | continuation of the same historical cases: no proper subgroup of index ≤ 7 and no nontrivial finite quotient of order ≤ 10⁵; diagnostic evidence, not a proof of triviality | `logs/phase3_out.txt`, manifest `logs/phase3_done.txt` (122 jobs, zero hits) | ~56 CPU-h |
| `kb_certify.g` (needs kbmag; run from `scripts/` or repo root) | the diagnostic rewriting certificates for diagram 1: Knuth–Bendix trivializes the full R₃-augmented grid (G1: 288/288), every generator reduced to the identity; positive control (surface group → Size ∞) and negative control (partial relations → no completion) | `logs/kb_certify_out.txt` | ~5 min |
| `kb_diag2_full.g` (needs kbmag) | the same for the second diagram, with the (true, diagram-independent) relation R₃: 576/576 certified trivial, all nine formal exponent pairs, all 64 sign conventions each | `logs/kb_diag2_full_out.txt` | ~10 min |
| `diag_g2_probe.g` (needs kbmag) | why diagram 2 needs R₃: representative hard cases are non-confluent without it and collapse (some to the empty presentation) with it | `logs/diag_g2_probe_out.txt` | ~1 min |
| `tc_deep.g` | historical pre-reindex deep-enumeration record: two formal systems exceed 10⁸ cosets without terminating; retained in the old coordinate for reproducibility | `logs/tc_deep_out.txt` | ~2 min |
| `maf_export.g` + `maf_certify.sh` | **independent-engine cross-check**: the 8 representative presentations (pre-repair word) re-decided by MAF (Alun Williams' Monoid Automata Factory; no shared code with GAP/kbmag): word acceptor = 1 word, every generator → IdWord; surface control infinite. Build MAF from https://sourceforge.net/projects/maffsa/ at `-O0` (an optimized arm64 build miscompiles) | `logs/maf_out.txt` | ~1 min |
| `maf_export2.g` + `maf_certify2.sh` | the current $y_1/Ax$ independent-engine run: 8 representative cases + surface control through MAF (author's official v2.2.1 binaries): word acceptor = 1 word, every generator → IdWord, surface control infinite | `logs/maf_out2.txt` | ~1 min |

All logs in `logs/` are the actual outputs of these scripts. Enumerations are
capped at 4×10⁵ cosets; an overflow means only that the enumeration exceeded
the cap and carries no group-theoretic conclusion.

(The phase-3 launcher is resumable: the manifest lists completed (case, target)
jobs and re-running skips them. `phase3_worker.g` documents, and soundly works
around, a GAP 4.16.0 quirk in `GQuotients`' `ExcludedOrders` preprocessing whose
hard-coded coset limit these presentations exceed.)

## How to verify

Appendix A.3 of the manuscript ("Verification guide") lists the geometric
inputs and the most delicate steps of the derivation chain, in order, with the
check to consult for each; and Appendix A.4 gives the protocol for diffing an
independently derived relation sheet against Table 1 (§9.1): quotient out
global generator inversions, the five swept signs, basing conjugation and the
arc identity, and any surviving mismatch is *the* disagreement, localized to
one crossing count, one conjugating path or one framing word. Readers wishing
to verify — or break — the computation should start there, then re-run the
scripts above against the expected outputs and the committed logs.

## Layout

```
papers/    the manuscript; the step-by-step walkthrough companion (PDF)
scripts/   GAP scripts + the development algorithm (Python)
logs/      outputs of every run referenced in the manuscript
docs/      GAP install notes
```

## License

The scripts are released under the MIT License (see `LICENSE`). The papers are
© the author.

## Contact

Bernd Johannes Wuebben, New York, NY — wuebben@gmail.com
