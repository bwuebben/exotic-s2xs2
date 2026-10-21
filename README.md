# Toward an exotic S²×S² and an exotic ℂP²#ℂP̄²

Bernd Johannes Wuebben · Working revision, October 21, 2026 · 53 pages

[Latest version (PDF)](paper/main.pdf) · [LaTeX source](paper/main.tex) · [arXiv version](https://arxiv.org/abs/2608.17267v3) · [Verification scripts](paper/anc/README.md)

This paper studies the specified Luttinger filling **V = W(+1,−1;0,0)**
in the construction of Lidman and Piccirillo. It derives the complete
fundamental-group presentation, gives reversible algebraic reductions,
and identifies the remaining question as a single normal-closure problem
in an iterated HNN extension containing the genus-two surface group.
It also determines the universal central extension with cyclic kernel,
constructs geometric realizations of it, and extends the specified
boundary involution over the filling.

The October 21 working revision adds a proposed capacity obstruction
from a Lagrangian torus, a simply connected comparison, and a conditional
discussion of symplectic homology after the auxiliary handles. These
additions in Sections 17 and 18 await independent mathematical
verification. The October 16 manuscript remains available as
[arXiv:2608.17267v3](https://arxiv.org/abs/2608.17267v3).

**The final filling problem remains unresolved.** The commutator J = [x,y]
normally generates the filled group G. Deciding whether the final filling
relation kills J would decide whether G is trivial. The conclusions about
exotic S²×S² and ℂP²#ℂP̄² remain conditional on triviality of G.

## Reading guide

- **Part I:** the marked geometric construction and its geometric consequences.
- **Part II:** the complete group presentation, intermediate groups and kernels,
  structural reductions, and the universal central extension.
- **Part III:** geometric models for the central extension, the proposed
  Liouville and symplectic-homology additions, and the precise remaining questions.
- **Appendices:** the explicit identities and finite algebraic data used in
  the proofs, together with supporting subgroup calculations.

## Build the paper

The source is self-contained: both figures are drawn with TikZ and the
bibliography is included in `main.tex`. With a standard TeX Live installation:

```sh
cd paper
latexmk -pdf main.tex
```

No scripts, external images, data files, or separate bibliography files are
needed to build or read the paper.

## Check the displayed algebra

Two optional Python scripts independently check the explicit finite
calculations printed in the manuscript. They require only Python 3's
standard library:

```sh
python3 -B paper/anc/verify.py
python3 -B paper/anc/verify_expanded.py
```

Both commands print JSON with `"verification": "PASS"`. Run without `-O`,
since the checks use assertions. The [verification guide](paper/anc/README.md)
describes their exact scope and the supplied expected outputs. These scripts
do not decide G or verify the new symplectic arguments.

## Contents and license

`paper/` contains the PDF, its source, and the ancillary verification files.
The scripts are available under the [MIT license](LICENSE). The paper is
copyright © 2026 Bernd Johannes Wuebben; the software license does not apply
to the manuscript.
