# An exotic S²×S² and an exotic ℂP²#ℂP̄²

> **Author's note:** Issues have been identified in the geometric derivation of the fundamental-group relations in Part II, so the claimed simple connectivity and its exotic-manifold consequences are not established by the current proof.
>
> The existence of an exotic S²×S² is a major open problem in four-dimensional topology: it asks whether one of the simplest four-dimensional spaces admits a different smooth structure. The problem has resisted decades of work. Earlier claimed constructions include that of Akhmedov and Park ([arXiv:1005.3346](https://arxiv.org/abs/1005.3346), 2010), which has not provided an established resolution. The present paper is another attempt at this difficult problem, and its current shortcomings are acknowledged here explicitly.
>
> The author remains committed to pursuing the problem and is currently working on a repair.

This repository accompanies Bernd Johannes Wuebben's paper
[*An exotic S²×S² and an exotic ℂP²#ℂP̄²*](papers/exotic-s2xs2-and-cp2.pdf),
posted as [arXiv:2608.17267v1](https://arxiv.org/abs/2608.17267) on
August 18, 2026.

The linked PDF is the September 27, 2026 revision with the author's note
directly below the abstract; the arXiv link identifies the original posting.

The repository contains the 33-page paper, an 11-page expository walkthrough,
and the scripts and run logs for every computation reported in the paper. The
manuscript is an arXiv preprint and has not been submitted to a journal.

## The claimed result

Lidman and Piccirillo construct two closed smooth 4-manifolds, $B$ and $W$,
with isomorphic integer cohomology rings. The figure-eight knot is smoothly
slice in $B$ but not in $W$. Their construction leaves open whether $B$ and
$W$ are homeomorphic.

The manifold $W$ is the quotient $V/\sigma$. The piece $V$ is obtained from
a genus-2 surface bundle by two Luttinger surgeries, and $\sigma$ is a free
involution on its boundary. The surgery curves must be specified before $V$
denotes a unique manifold. The paper fixes one explicit choice, written
$V=V'_{0,0}$, and claims

$$
\pi_1(V)=1.
$$

If established, this simple-connectivity statement would yield three conclusions:

1. The double $V\cup_\sigma V$ is homeomorphic but not diffeomorphic to
   $S^2\times S^2$.
2. The manifolds $B$ and $W$ are homeomorphic but not diffeomorphic, and the
   smooth structures are distinguished by whether the figure-eight knot is
   slice.
3. The regluing construction of Lidman and Piccirillo produces a simply
   connected manifold homeomorphic but not diffeomorphic to
   $\mathbb{CP}^2\sharp\overline{\mathbb{CP}}^2$.

The fundamental-group claim is for the explicitly parametrized piece
$V=V'_{0,0}$. The additional surgery parameters evaluated by some of the
scripts are consistency tests and are not additional manifold theorems.

## The proposed proof at a glance

The paper is organized in two parts.

**Part I: the topological reduction.** Sections 3--6 analyze the allowable
Luttinger surgeries and show how simple connectivity of $V$ gives the three
conclusions above. This part uses the topological classification results of
Freedman and Hambleton--Kreck together with the slicing obstruction from the
Lidman--Piccirillo construction.

**Part II: the fundamental-group computation.** Sections 7--12 build an
explicit model of the surface bundle and its two surgery tori. The genus-2
fiber is represented by a marked octagon and the base by a cut square. In
this model the meridians and Lagrangian push offs are written as based group
words, including the paths that connect every loop to a common basepoint.

These words define a finitely presented group $G$. Section 10 claims that
there is a surjection

$$
G\twoheadrightarrow\pi_1(V).
$$

GAP coset enumeration proves that $G$ is trivial. The calculation is repeated
over 4,096 choices of signs, paths, and correction placements, including the
stated relation system, and every case gives the trivial
group. Since every quotient of the trivial group is trivial, a valid surjection
would imply $\pi_1(V)=1$; the current geometric derivation does not establish it.

The scripts reproduce the algebraic calculations from the relations stated
in the paper. The proposed geometric derivation of those relations appears in
Sections 7--10 and is the subject of the author's note above. Appendix A connects
each relation to its proposed geometric source and corresponding verification record.

## A reading guide

There are several useful ways into the project.

### For a first overview

Read Section 1 of the
[*paper*](papers/exotic-s2xs2-and-cp2.pdf). It states the four main theorems,
defines the particular manifold $V$, and explains how the two parts of the
proof fit together.

### For the topology

Read Sections 2--6. They introduce the Luttinger-surgery conventions, give an
explicit presentation of the unsurgered surface bundle, establish rigidity of
the surgery parameters, and prove the reduction theorem.

### For the fundamental-group calculation

The [step-by-step walkthrough](papers/walkthrough.pdf) develops the same
calculation more slowly than the journal-style paper. It explains the marked
surface, the based loops, every input relation, the additional drilled-fiber
relation, and the claimed surjection onto $\pi_1(V)$.

The corresponding proposed argument in the paper is Sections 7--11:

- Section 7 fixes the marked octagon and cut-square model.
- Section 8 derives the meridians, transport relations, and surgery
  directions.
- Section 9 assembles the complete relation sheet.
- Section 10 claims that the relation group maps onto $\pi_1(V)$.
- Section 11 decides that relation group.

### For verification

Appendix A.1 explains the word-development algorithm. Appendix A.2 gives the
complete run table. Appendix A.3 lists the geometric inputs and the checks
attached to each one. Appendix A.4 gives a protocol for comparing an
independently derived relation sheet with Table 1 of the paper. Appendix B
works through two calibration examples whose fundamental groups are known by
other methods.

## Reproducing the main computation

### Requirements

- **GAP 4**. The computation was developed with GAP 4.16.0. The main
  theorem-facing script uses the core GAP library only.
- **Python 3**, using only the standard library, for the word-development
  calculation.
- The optional Knuth--Bendix scripts require the GAP package
  [kbmag](https://github.com/gap-packages/kbmag).
- The optional independent rewriting check uses
  [MAF](https://sourceforge.net/projects/maffsa/).

For a minimal GAP installation on macOS, see
[`docs/INSTALL_GAP.md`](docs/INSTALL_GAP.md). On Linux, GAP is generally
available through the system package manager.

### The algebraic group decision

From the repository root, run:

```bash
gap -q -A scripts/fixed_v_certify.g
```

The script uses core GAP. Its two summary lines are:

```text
FIXED V (y1/Ax) WITH R3: TOTAL=4096 TRIVIAL=4096 OVERFLOW=0 FINITE>1=0 H1nonzero=0
ADJACENT T_ALPHA SECTION (y2/Ar^-1) WITH R3: TOTAL=4096 TRIVIAL=4096 OVERFLOW=0 FINITE>1=0 H1nonzero=0
```

The first line is the calculation used in the claimed theorem. In that label,
`y1/Ax` identifies the based surgery direction chosen for $V=V'_{0,0}$. The
label `R3` refers to the additional relation whose geometric derivation is
claimed in Section 10. The second line records a neighboring
choice included as a consistency check. Neither calculation validates the
geometric derivation. The
committed output is
[`logs/fixed_v_certify_out.txt`](logs/fixed_v_certify_out.txt).

### The word-development calculation

Run:

```bash
python3 scripts/develop.py
```

The output should agree with
[`logs/develop_out.txt`](logs/develop_out.txt). The program translates each
recorded path through the octagon into a group word and performs the five
internal validations described in Appendix A.1.

## The remaining verification files

The complete inventory and run tables are in Appendix A of the paper. The
principal files are grouped below by purpose.

| Purpose | Principal files |
|---|---|
| Claimed main theorem | `fixed_v_certify.g`, `develop.py` |
| Surface-bundle model and controls | `monodromy_check2.g`, `model_check3.g`, `vr_check.g` |
| Baldridge--Kirk calibration | `decide_t4.g`, `gt1_diff.g` |
| Seifert-fibered calibration | `decide_seifert.g`, `confirm_coherent.g` |
| Alternative relation systems and based diagram | `decide.g`, `decide2.g`, `vdiag2.g`, `placement_check.g` |
| Knuth--Bendix decisions | `kb_certify.g`, `kb_diag2_full.g` |
| Sensitivity experiments | `pi1_grid.g`, `pi1_v2b.g`, `vsens.g` |
| Independent MAF check | `maf_export*.g`, `maf_certify*.sh` |
| Preserved longer diagnostics | `diag_*`, `phase2_*`, `phase3_*`, `tc_deep.g` |

All files in `logs/` are the recorded outputs of the corresponding runs.
Coset enumerations are capped at 400,000 cosets. An `OVERFLOW` means only that
the enumeration exceeded the cap; it does not imply that the group is
nontrivial. The main checks take minutes on a laptop. The optional historical
finite-quotient search represented by the `phase2_*` and `phase3_*` files took
about 56 CPU-hours; its completed logs are included.

## Repository layout

```text
papers/    the manuscript and step-by-step walkthrough (PDF)
scripts/   GAP scripts and the Python word-development program
logs/      outputs of every run cited in the manuscript
docs/      GAP installation notes
```

The source manuscript of record is maintained in the author's research
monorepo. This repository is the public verification mirror.

## License

The scripts are released under the MIT License; see [`LICENSE`](LICENSE).
The papers are © the author.

## Contact

Bernd Johannes Wuebben, New York, NY — wuebben@gmail.com
