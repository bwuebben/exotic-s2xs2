# Paper

**Toward an exotic S²×S² and an exotic ℂP²#ℂP̄²**

Bernd Johannes Wuebben · Working revision, October 21, 2026 · 56 pages, 2 figures

- [Latest version (PDF)](main.pdf)
- [Self-contained LaTeX source](main.tex)
- [arXiv version, October 16](https://arxiv.org/abs/2608.17267v3)
- [Exact algebra verification](anc/README.md)

This is a partial-results paper. It develops the specified filling
W(+1,−1;0,0), its complete fundamental-group presentation, structural
reductions, and geometric models. The final filling question is open;
the exotic homeomorphism conclusions are conditional.

Sections 17--19 add a proposed Lagrangian-torus capacity obstruction,
a simply connected comparison, exact tori detecting a transfer kernel,
and a conditional handle-surgery argument for symplectic homology. These
additions await independent mathematical verification, as indicated on
the first page. Section 20 gives an elementary group-algebra criterion
for infinite order of the central element, without proving that condition.
The handle discussion retains separate localization and SFT/Hamiltonian
comparison assumptions. This working revision has
not replaced the October 16 arXiv manuscript.

Build from this directory with `latexmk -pdf main.tex`. The figures and
bibliography are included in the source. No ancillary script is needed
for the build.

The two scripts in `anc/` check finite identities and other printed
algebraic data. They require only the Python 3 standard library and do
not determine whether the filled group is trivial or verify the new
symplectic arguments.
