# Barker's order-10 recurrence for OEIS A321614, proved

A proof, for all `n`, of the order-10 linear recurrence and rational generating
function that Colin Barker conjectured in 2018 for
[OEIS A321614](https://oeis.org/A321614) — the number of ways to place `2n`
nonattacking kings on a `4 × 2n` board attaining the maximum `2n`, counted up to
the symmetry group of the rectangle. The OEIS entry still records the recurrence
and generating function as conjectures.

The paper is `a321614.tex`.

## What this repository is for

Everything needed to reproduce every number in the paper, and to check the
paper against the code. It is intended to be run and read by someone who did not
write it.

## Requirements

Python 3, and [uv](https://docs.astral.sh/uv/) (or plain `pip`). One script
needs `sympy`; the rest use the standard library only. No LaTeX toolchain is
needed to run any check — only to typeset the paper.

## Reproduce

```
uv run python verify_model.py            # the combinatorial links, exhaustively
uv run python verify_burnside.py         # each Burnside term, separately
uv run python kings.py check             # all 22 published terms
uv run python crosscheck.py 26           # a second, independent implementation
uv run --with sympy python prove.py      # the generating-function identity
uv run --with sympy python verify_paper.py   # this repo vs the paper
uv run python check_tex.py               # LaTeX structure + journal style rules
```

Each of these prints its own verdict and exits non-zero on failure.

## The proof in one paragraph

The maximum condition forces exactly one king into each of the `2n` disjoint
`2 × 2` blocks of the board, which turns the placements into walks in a 12-state
transfer graph. Burnside's lemma over the Klein four-group expresses `a(n)` as
the average of four fixed-point counts, each a linear combination of entries of
powers of the transfer matrix. Summing the four rational generating functions in
exact arithmetic reproduces Barker's generating function identically, so the
recurrence holds for all `n ≥ 10`. The only machine-assisted step is that final
rational identity; every combinatorial step is proved by hand in the paper, and
each is also closed by exhaustive finite check here.

## Two checks that are easy to skip, and should not be

**The Burnside terms are checked individually, not just their average.**
Burnside's lemma is a statement about each fixed-point count on its own. Two
terms can be wrong in opposite directions and cancel, or all four can be wrong
together in a way that preserves the average. `verify_burnside.py` compares
`T`, `F_h`, `F_v` and `F_r` each against definition-level enumeration.

**The paper is checked against the code.** Every other script here reads the
code. A matrix mistyped while being copied into the paper would be invisible to
all of them, so `verify_paper.py` re-reads `a321614.tex`, parses the transfer
matrix, the state list and the weight vectors out of the LaTeX, and compares
them with what the programs compute. That check caught a genuine derivation
error in the write-up during preparation: the `σ_v` fixed-point condition had
been stated as `r₀ = r₁` instead of `r₁ = 1 − r₀`. The conclusion was unaffected;
the reasoning was not, and a reader following the wrong condition would have got
stuck on the state `(1,1,1,1)`.

## Files

| file | role |
|---|---|
| `a321614.tex` | the paper |
| `kings.py` | block transfer matrix; regression against the b-file; brute force |
| `verify_model.py` | closes the combinatorial links exhaustively |
| `verify_burnside.py` | each Burnside term separately vs brute force |
| `crosscheck.py` | a second implementation that never uses the block structure |
| `prove.py` | the generating-function identity (needs `sympy`) |
| `verify_paper.py` | the paper's transcriptions vs the code |
| `check_tex.py` | LaTeX structure and journal style rules |
| `dump_M.py` | prints the explicit matrix and vectors |
| `terms.txt` | `n, a(n)` for `n = 0..600` |
| `b321614.txt` | the published 22-term b-file, used as a regression gate |

## Scope, stated plainly

This proves a specific 2018 conjecture using standard tools — transfer matrix,
Burnside's lemma, exact rational generating functions. It is not a new method,
and the sequence and the conjecture are not the author's. The numerical
extension of the sequence to `n = 5000` was carried out independently and
earlier by others and is not claimed here.

## License

MIT. See `LICENSE`.
