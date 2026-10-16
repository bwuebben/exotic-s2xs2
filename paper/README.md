# Paper

**Toward an exotic S²×S² and an exotic ℂP²#ℂP̄²**

Bernd Johannes Wuebben · October 16, 2026 · 49 pages, 2 figures

- [PDF](main.pdf)
- [Self-contained LaTeX source](main.tex)
- [Exact algebra verification](anc/README.md)

This is a partial-results paper. It develops the specified filling
W(+1,−1;0,0), its complete fundamental-group presentation, structural
reductions, and geometric models. The final filling question is open;
the exotic homeomorphism conclusions are conditional.

Build from this directory with `latexmk -pdf main.tex`. The figures and
bibliography are included in the source. No ancillary script is needed
for the build.

The two scripts in `anc/` check finite identities and other printed
algebraic data. They require only the Python 3 standard library and do
not determine whether the filled group is trivial.
