# exotic-s2xs2

Scripts, certificates and run logs for the papers

- **Distinguishing homeomorphic 4-manifolds by slicing and the exotic S²×S²
  problem** (`papers/reduction.pdf`), and
- **Simple connectivity of the Lidman–Piccirillo piece: a certificate-complete
  computation** (`papers/computation.pdf`),

by Bernd J. Wuebben (arXiv links to follow). A step-by-step expository
companion to the computation — every corrected relation and every direction
word worked end to end at tutorial grain, closing with the Baldridge–Kirk T⁴
configuration as a fully worked calibration example — is
`papers/walkthrough.pdf`. Everything needed to reproduce every
number in both papers is here; total runtime is minutes on a laptop (except the
optional finite-quotient sweep, ~56 CPU-hours).

## The result

Lidman–Piccirillo (arXiv:2505.14387) built spin rational homology 4-spheres $B$
and $W$ with the same integer cohomology ring, distinguished by whether the
figure-eight knot is slice. The **reduction paper** proves that making that pair
*homeomorphic* is equivalent to the exotic $S^2\times S^2$ problem: for any
admissible surgery variant $V'$ of their key piece, $W'\cong_{\rm homeo}B$ iff
$\pi_1(V'\cup_\sigma V')=1$, in which case the double is an exotic
$S^2\times S^2$. The **computation paper** then computes, from a fully explicit
based model with machine-checkable certificates, complete presentations of
$\pi_1(V')$ for the whole admissible family and decides them by coset
enumeration and Knuth–Bendix completion:

> **π₁(V′) = 1 for the Lidman–Piccirillo piece and for all eight variant
> cells** — the full admissible grid, in every sign, placement and diagram
> convention swept.

The derivation pipeline is calibrated on a configuration with an independently
known answer — the Baldridge–Kirk double-Luttinger tori in T⁴
(`decide_t4.g`): with correctly derived lasso data the harness reproduces the
published non-abelian fundamental groups in all 64 conventions, and every
wrong lasso pair fails all 64, half of them abelian — the signature of a
missed whisker crossing. That calibration found (and led to the repair of)
exactly one such error in an earlier version of one direction word
(`dir_Tβ^base`; the computation paper's "pushoff-basing correction" section
documents the found-and-fixed derivation), after which **every verdict was
re-established with the corrected, sign-coupled word**.

Six cells are decided by coset enumeration; the corrected word makes the rest
enumeration-blind — past 10⁸ cosets on the two hardest cells (`tc_deep.g`),
which also survived a 56-CPU-hour finite-quotient search in their pre-repair
form (the phase-2/3 scripts below: perfect, no subgroup of index ≤ 7, no
nontrivial finite quotient of order ≤ 10⁵). All of them are decided by
**Knuth–Bendix completion**: the confluent rewriting system of
each presented group **reduces every generator to the identity**, a triviality
certificate independent even of the confluence claim. The certificates
re-derive the entire completed-presentation grid (`kb_certify.g`, 288/288)
and the entire second diagram (`kb_diag2_full.g`, 576/576 — all nine cells,
all 64 conventions each), with positive and negative controls. A second,
independently implemented engine (MAF, `maf_certify.sh`) reproduced the eight
representative verdicts of the pre-repair word; the corrected-word exports
for the identical MAF run are staged (`maf_export2.g`).

Together: an exotic $S^2\times S^2$, the first pair of homeomorphic closed
4-manifolds distinguished by unconstrained knot slicing, and a simply connected
exotic $\mathbb{CP}^2\sharp \overline{\mathbb{CP}}^2$.

## Requirements

- **GAP 4** (developed on 4.16.0; core library only, except that
  `kb_certify.g` additionally needs the **kbmag** package, built from
  https://github.com/gap-packages/kbmag into GAP's `pkg/` directory). GAP is
  not in Homebrew or MacPorts; see `docs/INSTALL_GAP.md` for a minimal source
  build on macOS, or use your distribution's package on Linux
  (`apt install gap`).
- **Python 3** (standard library only) for the developing engine.

Run any script with

```
gap -q -A scripts/<name>.g
python3 scripts/develop.py
```

## What each script certifies

| script | certifies | expected output | time |
|---|---|---|---|
| `monodromy_check2.g` | the trefoil monodromy model (fixes [x,y] exactly; h⁶ = boundary twist; h³ = −id) | `true / true / false` + h³ images | <1 s |
| `model_check3.g` | h = T_a∘T_b (order pinned); relator preservation; ψφψ⁻¹φ = h∗h⁻¹; fiber-H₁ Smith form | all `true`, diagonal `[1,1,1,1]` | <1 s |
| `develop.py` | the developing engine + validations V1–V4; derived words: a, b, d, e; c ≃ (rx)⁻¹; z ≃ sy; lasso δ = r⁻¹; κ₃ = s⁻¹r⁻¹yx | all `OK`; see `logs/develop_out.txt` | <1 s |
| `pi1_grid.g` | experiment E1: all 72 uncorrected candidate presentations collapse | 72 × `\|G\| = 1` | ~1 s |
| `pi1_v2b.g` | controls + E2 sensitivity (120 seeded trials) | `SUMMARY C: TRIVIAL=58 FINITE=0 H1=7 BLOWUP=55` | ~5 s |
| `ap_check.g` | E3: Akhmedov–Park Lemma 8 ⇒ π₁(Mₙᵖ) = ℤ/p | ℤ/1, ℤ/1, ℤ/1, ℤ/2, ℤ/3 | ~6 min |
| `decide_t4.g` | **known-answer calibration** on the Baldridge–Kirk T⁴ configuration: derived lassos (1, b) pass 64/64 with the exact tr 3/tr 1 fingerprint split; every wrong lasso pair fails 0/64 (half abelian); single surgery = H₃(ℤ)×ℤ 8/8; BK's own words cross-check | `logs/t4_out.txt` | ~20 min |
| `decide.g` | phase-1 harness (self-test mode without inputs; honest 288-case sweep with the derived, sign-coupled direction words) | see `logs/honest_run_2.log` (pre-repair record: `honest_run_1.log`) | ~1 min |
| `diag_dirTbBase.g`, `diag_kb.g`, `diag_fullgrid.g`, `diag_r3.g` | the found-and-fixed audit trail: all 8 candidate resolutions of the pushoff-basing correction at the LP cell (enum + KB certification of every blowup, 256/256 trivial) and both lasso arcs across the full completed grid (576/576 certified trivial) | `logs/diag_*_out.txt` | ~30 min total |
| `placement_check.g` | correction-placement robustness at the LP cell | `TRIVIAL=192 BLOWUP=64 FINITE>1=0 H1nonzero=0` (blowups decided by the KB certificates) | ~1 min |
| `vdiag2.g` | the independently derived second diagram (576 cases; enumeration only — the corrected word is enum-hostile here, H₁ = 0 in all, zero finite>1) | see `logs/vdiag2_out2.txt`; certification in `kb_diag2_full.g` | ~2 min |
| `vr_check.g` | known-answer probes (no fillings → ℤ²; fiber-only → ℤ²; one filling → ℤ) | `[0,0] / [0,0] / [0]` | <1 s |
| `decide2.g` | the **completed** presentation (adds R₃): the final sweep | `logs/decide2_out2.txt`: 6 of 9 cells trivial in all 32 conventions incl. the LP cell; the n=+1 column blows up (decided by KB below) | ~2 min |
| `phase2_parallel.sh` (+ `phase2_common.g`, `phase2_worker.g`) | the two open cells, 8 representative sign cases: Tietze to 3 gen/7 rel; enumeration past 8×10⁶ cosets; no subgroup of index ≤ 6; no quotient onto the 17 simples ≤ 14880 (except U₃(3)) | `logs/phase2_par_out.txt` | ~6 min (8-wide) |
| `phase3_resume.sh` (+ `phase3_worker.g`; `phase3_parallel.sh` from scratch) | same 8 cases: no subgroup of index ≤ 7; no quotient onto the remaining simples ≤ 10⁵ (U₃(3), A₈, L₃(4), …, U₃(4), L₂(53), M₁₂) ⇒ **no nontrivial finite quotient of order ≤ 10⁵** | `logs/phase3_out.txt`, manifest `logs/phase3_done.txt` (122 jobs, zero hits) | ~56 CPU-h |
| `kb_certify.g` (needs kbmag; run from `scripts/` or repo root) | **the decisive certificates**: Knuth–Bendix trivializes the full completed-presentation grid (G1: 288/288), every generator reduced to the identity; positive control (surface group → Size ∞) and negative control (partial relations → no completion). Its G2 block (diagram 2 without a completion relation) is KB-inconclusive with the corrected word — 0/128 reach confluence, never an adverse verdict (`diag_g2_probe.g` diagnoses this; the diagram-2 certification is the next row) | `logs/kb_certify_out2.txt` (pre-repair: `kb_certify_out.txt`) | ~5 min |
| `kb_diag2_full.g` (needs kbmag) | **the second-diagram certificates**: the full diagram-2 grid with the (true, diagram-independent) completion relation R₃ — 576/576 certified trivial, all nine cells, all 64 sign conventions each | `logs/kb_diag2_full_out.txt` | ~10 min |
| `diag_g2_probe.g` (needs kbmag) | why diagram 2 needs R₃: representative hard cases are non-confluent without it and collapse (some to the empty presentation) with it | `logs/diag_g2_probe_out.txt` | ~1 min |
| `tc_deep.g` | the deep-enumeration record: both former blowup cells exceed 10⁸ cosets without terminating — on groups now known trivial | `logs/tc_deep_out.txt` | ~2 min |
| `maf_export.g` + `maf_certify.sh` | **independent-engine cross-check**: the 8 representative presentations (pre-repair word) re-decided by MAF (Alun Williams' Monoid Automata Factory — no shared code with GAP/kbmag): word acceptor = 1 word, every generator → IdWord; surface control infinite. Build MAF from https://sourceforge.net/projects/maffsa/ at `-O0` (an optimized arm64 build miscompiles) | `logs/maf_out.txt` | ~1 min |
| `maf_export2.g` | the corrected-word exports (8 representative cases + surface control, into `maf_runs2/`) staged for the identical MAF run | writes `maf_runs2/` | <1 min |

All logs in `logs/` are the actual outputs of these scripts; every presented
group in every sweep has H₁ = 0, and no run anywhere produced a finite
nontrivial group.

(The phase-3 launcher is resumable: the manifest lists completed (case, target)
jobs and re-running skips them. `phase3_worker.g` documents — and soundly works
around — a GAP 4.16.0 quirk in `GQuotients`' `ExcludedOrders` preprocessing whose
hard-coded coset limit these presentations exceed.)

## How to verify

The computation paper's final section ("Verification guide") lists the most
delicate steps of the derivation chain, in order, with the certificate to check
for each. Readers wishing to verify — or break — the computation should start
there, then re-run the scripts above against the expected outputs and the
committed logs.

## Layout

```
papers/    the two papers (PDF)
scripts/   GAP scripts + the developing engine (Python)
logs/      outputs of every run referenced in the papers
docs/      GAP install notes
```

## License

The scripts are released under the MIT License (see `LICENSE`). The papers are
© the author.

## Contact

Bernd J. Wuebben, New York, NY — wuebben@gmail.com
