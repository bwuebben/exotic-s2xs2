# Exact algebra verification

Ancillary files for *Toward an exotic S²×S² and an exotic ℂP²#ℂP̄²*,
Bernd Johannes Wuebben. The October 16 algebra checks also apply to
the unchanged finite algebraic data in the October 21 working revision.

These checks supplement the proofs and explicit data in the manuscript.
They require Python 3 and its standard library only. Keep the two scripts
in the same directory: `verify_expanded.py` imports word and integer
arithmetic from `verify.py`.

Run from this directory:

```sh
python3 -B verify.py
python3 -B verify_expanded.py
```

Do not use Python's `-O` option: the checks use assertions. Each script
prints JSON to standard output and creates no files. Expected outputs
are supplied as `verification.json` and `verification_expanded.json`.

## What is checked

`verify.py` checks:

- The explicit identities underlying the fifteen-to-five presentation
  equivalence, in both directions, and the subsequent two-row replacement.
- The crossing substitution and its inverse.
- A unimodular abelianization minor and the rank-three homology triples.
- The cotangent normal-coordinate and reflection matrices.

`verify_expanded.py` checks:

- The nineteen-factor central-extension identity and lift matrix.
- Additional crossing and boundary identities, including the inverse crossing.
- The trefoil semidirect-product dictionary and mixed-sign conjugator.
- The printed 30-state labelled action and folded subgroup intersection.
- All 96 Whitehead cut differences.
- The exterior algebra of the symplectic-pair correction.

These are exact finite calculations, not searches for finite quotients.
The programs do not prove geometric inclusions, HNN injectivity, or
four-manifold classification, and do not decide triviality or
nontriviality of the final group. The geometric and group-theoretic
arguments establishing the roles of these calculations are in the paper.
The checks do not verify the proposed symplectic arguments added in
Sections 17--19 of the working revision or the new group-algebra argument
in Section 20.

The scripts are distributed under the MIT license in `LICENSE.txt`.
