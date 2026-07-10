# exotic-s2xs2

Scripts, certificates and run logs for the papers

- **Distinguishing homeomorphic 4-manifolds by slicing and the exotic S²×S²
  problem** (`papers/reduction.pdf`), and
- **Simple connectivity of the Lidman–Piccirillo piece: a certificate-complete
  computation** (`papers/computation.pdf`),

by Bernd J. Wuebben (arXiv links to follow). Everything needed to reproduce every
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
enumeration:

> **π₁(V′) = 1 for the Lidman–Piccirillo piece and for all eight variant
> cells** — the full admissible grid, in every sign, placement and diagram
> convention swept.

Seven cells are decided by coset enumeration. The two cells that resist it —
past 10⁸ cosets (`tc_deep.g`), having also survived a 56-CPU-hour
finite-quotient search (the phase-2/3 scripts below: perfect, no subgroup of
index ≤ 7, no nontrivial finite quotient of order ≤ 10⁵) — are decided by
**Knuth–Bendix completion** (`kb_certify.g`): the confluent rewriting system of
each presented group **reduces every generator to the identity**, a triviality
certificate independent even of the confluence claim. The same certificates
re-derive the entire grid (288/288) and pass the independent second-diagram
cross-check (128/128), with positive and negative controls — and a second,
independently implemented engine (MAF, `maf_certify.sh`) reproduces all eight
verdicts.

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
| `decide.g` | phase-1 harness (self-test mode without inputs; honest 288-case sweep with the derived direction words) | see `logs/honest_run_1.log` | ~1 min |
| `placement_check.g` | correction-placement robustness at the LP cell | `TRIVIAL=256 BLOWUP=0 FINITE>1=0 H1nonzero=0` | ~1 min |
| `vdiag2.g` | the independently derived second diagram (576 cases) | see `logs/vdiag2_out.txt`; LP cell trivial 64/64 | ~2 min |
| `vr_check.g` | known-answer probes (no fillings → ℤ²; fiber-only → ℤ²; one filling → ℤ) | `[0,0] / [0,0] / [0]` | <1 s |
| `decide2.g` | the **completed** presentation (adds R₃): the final sweep | `logs/complete_run.log`: 7 of 9 cells trivial in all 32 conventions | ~2 min |
| `phase2_parallel.sh` (+ `phase2_common.g`, `phase2_worker.g`) | the two open cells, 8 representative sign cases: Tietze to 3 gen/7 rel; enumeration past 8×10⁶ cosets; no subgroup of index ≤ 6; no quotient onto the 17 simples ≤ 14880 (except U₃(3)) | `logs/phase2_par_out.txt` | ~6 min (8-wide) |
| `phase3_resume.sh` (+ `phase3_worker.g`; `phase3_parallel.sh` from scratch) | same 8 cases: no subgroup of index ≤ 7; no quotient onto the remaining simples ≤ 10⁵ (U₃(3), A₈, L₃(4), …, U₃(4), L₂(53), M₁₂) ⇒ **no nontrivial finite quotient of order ≤ 10⁵** | `logs/phase3_out.txt`, manifest `logs/phase3_done.txt` (122 jobs, zero hits) | ~56 CPU-h |
| `kb_certify.g` (needs kbmag; run from `scripts/` or repo root) | **the decisive certificates**: Knuth–Bendix trivializes the full grid (288/288) and the second diagram (128/128), every generator reduced to the identity; positive control (surface group → Size ∞) and negative control (partial relations → no completion) | `logs/kb_certify_out.txt` | ~3 min |
| `tc_deep.g` | the deep-enumeration record: both former blowup cells exceed 10⁸ cosets without terminating — on groups now known trivial | `logs/tc_deep_out.txt` | ~2 min |
| `maf_export.g` + `maf_certify.sh` | **independent-engine cross-check**: the 8 representative presentations re-decided by MAF (Alun Williams' Monoid Automata Factory — no shared code with GAP/kbmag): word acceptor = 1 word, every generator → IdWord; surface control infinite. Build MAF from https://sourceforge.net/projects/maffsa/ at `-O0` (an optimized arm64 build miscompiles) | `logs/maf_out.txt` | ~1 min |

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
