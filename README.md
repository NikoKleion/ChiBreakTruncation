# ChiBreakTruncation

Python package `chibreak`. It builds peaked circuits on a line of qubits, simulates them with a dense
statevector, and simulates them as a matrix product state with a truncation rule you can swap. From that it measures
chi_break, the smallest bond dimension at which reading the signs of the single-qubit marginals
returns the peak bitstring.

Use it to

- try a truncation rule of your own against top-chi on the same circuits
- make peaked circuits with a known answer for testing an attack
- find how deep a random block has to be before its middle cut is near volume law
- see how close a recovery is to failing

## Why

A peaked circuit hides one bitstring in its output. A quantum computer returns it in one run. A classical
simulation has to search for it, and the cost of that search is the evidence for quantum advantage.

The published cost is chi_break: how large the bond dimension of a matrix product state has to be before
the simulation reads the hidden string off the signs of the single-qubit marginals. Every published
chi_break was measured with one truncation rule, top-chi, which keeps the largest Schmidt coefficients at
each step because that is the best choice for state fidelity.

Reading the string takes the sign of n marginals, a weaker requirement than state fidelity. A rule that
keeps components for that purpose might reach the string at a smaller chi. If it does, published chi_break
values overstate how hard these circuits are. If it does not, chi_break holds up under a change of rule,
and that is worth knowing too. The near-flat Schmidt spectrum of these circuits is what makes this an open
question: when the coefficients are close in size, keeping the largest ones has no special claim.

I built this tool to test that. It needs three things: circuits with a peak you can weaken by one
parameter and an exact answer to compare against, a random block deep enough that the spectrum is in the
near-flat regime, and an MPS where the truncation rule is one function you can swap. The saved runs cover
the baseline: how deep the block must be, where the peak collapses as the mirror is broken and how the
smallest marginal moves before recovery fails, and how top-chi behaves against chi with a random rule as
a control.

## Install

```bash
pip install -e .
```

Python 3.11 or later, numpy and scipy.

## Use

```bash
python -m chibreak                       # list the runs
python -m chibreak chi_break             # print one run
python -m chibreak chi_break --save      # write results/chi_break.{md,csv}
python tests/run_tests.py                # tests, no pytest needed
```

```python
from chibreak import build_peaked_circuit, simulate_exact, simulate_mps, bits_from_marginals, true_peak

c = build_peaked_circuit(n=16, depth=24, alpha=0.05, seed=0)
peak, weight = true_peak(simulate_exact(c), 16)
m = simulate_mps(c, chi=32, rule="top_chi")
R = (bits_from_marginals(m.z_marginals()) == peak).mean()
```

[USAGE.md](USAGE.md) has the terms, the truncation-rule interface, the output columns and the limits.
[RESULTS.md](RESULTS.md) describes the three saved runs in `results/`.

## Contents

- `chibreak/circuit.py` the circuit `C = X^s B(alpha) A`: block A is `depth` layers of Haar-random
  two-qubit gates on nearest neighbours, block B is A reversed and inverted, s is a planted bitstring.
  alpha in [0, 1] moves each gate of B along the unitary geodesic from its exact inverse toward an
  independent random gate, so the peak weakens as alpha grows. At alpha = 0 the output is |s>.
- `chibreak/statevector.py` exact simulation, the marginals `<Z_i>`, the peak and its probability, the
  Schmidt spectrum across a cut.
- `chibreak/mps.py` the matrix product state. Every two-qubit gate ends in an SVD at its bond and a
  truncation rule picks which Schmidt components stay under a cap chi. `top_chi` keeps the largest, the
  rule behind published chi_break values. `random` keeps a uniform random subset, a control.
- `chibreak/metrics.py` entanglement entropy, effective rank, the smallest marginal, recovery.
- `chibreak/experiments.py` the runs `depth_calibration`, `alpha_sweep` and `chi_break`.
- `results/` saved output of each run: a markdown table, a csv, and per-instance rows.
- `tests/` checks against Kronecker products and exact states, and the saved numbers.

## Limits

- One-dimensional chain. Published chi_break values come from all-to-all circuits, so numbers here are
  not comparable to them.
- chi runs on powers of two, so chi_break resolves factors of two.
- Exact simulation needs 16 * 2^n bytes: 16 MB at n = 20, 4 GB at n = 28.
- alpha changes the peak weight and the mirror structure together.
- Sizes up to about 26 qubits on a chain. The published challenge instances are out of reach.

## References

- Aaronson and Zhang, arXiv:2404.14493. Peaked circuits and the random-block plus peaking-block construction.
- Gharibyan et al., arXiv:2510.25838. chi_break, top-chi truncation, and recovery from marginal signs.
- Harrow, Lowe and Witteveen, arXiv:2510.08518. Randomized truncation; not implemented here.
- Kremer and Dupuis, arXiv:2604.21908. Classical simulation of the published instances by MPO unswapping.
- Ferreira, arXiv:2607.07816. Truncated sparse statevector simulation of peaked circuits.
- Reweighted TEBD, arXiv:2412.08730. Observable-weighted SVD truncation.
- Eckart and Young, Psychometrika 1 (1936). Optimality of keeping the largest singular values.
- Page, Phys. Rev. Lett. 71, 1291 (1993). Average entanglement entropy of a random state, used as the
  reference value in the depth calibration.

## Thanks

To my teammates from team MerQury at the 2026 BlueQubit Yale peaked-circuit hackathon, Roman Bagdasarian,
Amon Koike, Zeenat Mayoon and Leong Wei Chan: thank you for one long night and the ideas you shared. And
to Tan Jun Liang (poig) for his public solutions and notes from earlier peaked-circuit hackathons. (:

## License

PolyForm Noncommercial 1.0.0, see [LICENSE.md](LICENSE.md). Copyright 2026 Nikolas Klein.
