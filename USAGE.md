# Usage

## Install

```bash
pip install -e .            # numpy and scipy
pip install -e ".[dev]"     # adds pytest
```

## Terms

- **Peaked circuit, peak weight.** A circuit is delta-peaked when one output bitstring has probability at
  least delta (Aaronson and Zhang, arXiv:2404.14493). `peak_weight` is the probability of the most likely
  bitstring, the true peak. The planted string s is the peak at alpha = 0 and can differ from the true peak
  for alpha > 0; the runs report both `R` (against the true peak) and `R_planted` (against s).
- **Marginal readout, R.** Bit i of the recovered string is 1 when `<Z_i> < 0`. R is the fraction of
  recovered bits equal to the target (Gharibyan et al., arXiv:2510.25838, where it is the recovery ratio).
- **Bond dimension chi, top-chi truncation.** An MPS with bond dimension chi keeps chi Schmidt components at
  each bond. Keeping the chi largest is optimal in the 2-norm (Eckart and Young, 1936).
- **chi_break.** The smallest bond dimension at which the marginal readout returns the peak, R = 1
  (Gharibyan et al.). Here chi takes the values 2, 4, 8, ... so chi_break is reported on that ladder.
- **min_abs_Z.** The smallest `|<Z_i>|` over qubits. **frac_abs_Z_below_0.25** is the fraction of qubits
  with `|<Z_i>|` below 0.25.
- **Entanglement entropy, effective rank.** The von Neumann entropy in bits of the Schmidt weights at the
  balanced cut between qubits [0, n/2) and [n/2, n), and the participation ratio `1 / sum p_i^2` of those
  weights. The entropy has maximum n/2 bits at that cut; Page's value for a Haar-random state is
  `n/2 - 1/(2 ln 2)` bits (Page, 1993).
- **chance_peak_weight.** The expected largest of 2^n Porter-Thomas weights, `(ln 2^n + 0.5772) / 2^n`,
  the peak weight of a state with no peak.
- **alpha, depth, seed.** Parameters of this generator: the mirror-breaking parameter in [0, 1], the number
  of brickwork layers in block A, and the RNG seed.

## Build a circuit and take the exact readout

```python
import numpy as np
from chibreak import (build_peaked_circuit, simulate_exact, state_after_block_a, z_marginals,
                      bits_from_marginals, true_peak, schmidt_spectrum, entropy_bits, effective_rank,
                      min_abs_z, fraction_below, recovery)

n, depth = 16, 24
c = build_peaked_circuit(n, depth, alpha=0.1, seed=0)

psi = simulate_exact(c)
z = z_marginals(psi, n)
bits = bits_from_marginals(z)
peak, weight = true_peak(psi, n)

print(recovery(bits, peak), recovery(bits, c["peak"]), min_abs_z(z), fraction_below(z, 0.25))

sv = schmidt_spectrum(state_after_block_a(c), n // 2, n)
print(entropy_bits(sv), effective_rank(sv))
```

`build_peaked_circuit` returns a dict: `gates` is a list of `(4x4 matrix, (i, i+1))` applied left to
right, all of A then all of B; `x_layer` and `peak` hold s; `a_gate_count` marks the end of A.

## Simulate with a bond cap

```python
from chibreak import simulate_mps, chi_ladder
from chibreak.mps import fidelity

for chi in chi_ladder(n, max_chi=64):
    m = simulate_mps(c, chi, rule="top_chi")
    R = recovery(bits_from_marginals(m.z_marginals()), peak)
    print(chi, R, min_abs_z(m.z_marginals()), fidelity(psi, m), m.discarded)
```

`simulate_mps(circuit, chi, rule, seed)` runs A, B and the X layer with every bond capped at chi.
`m.discarded` sums the truncated Schmidt weight over the run. `fidelity` contracts the MPS to a dense
vector, so keep it to n of about 16 or less.

## Write a truncation rule

A rule receives the singular values of one bond in descending order, the cap chi and a numpy
`Generator`, and returns the indices to keep. Register it under a name and pass that name.

```python
import numpy as np
from chibreak import TRUNCATION_RULES, simulate_mps

def alternate(s, chi, rng):
    # every other component, an example that ignores magnitude
    return np.arange(0, s.size, 2)[:chi]

TRUNCATION_RULES["alternate"] = alternate
m = simulate_mps(c, chi=16, rule="alternate")
```

The state is in mixed canonical form at the bond being truncated, so for `top_chi` the discarded weight
equals the sum of the dropped squared singular values.

## Run the saved experiments

```bash
python -m chibreak                                   # list runs
python -m chibreak depth_calibration --save
python -m chibreak alpha_sweep --save
python -m chibreak chi_break --save
python -m chibreak alpha_sweep chi_break --quick     # small smoke run
```

`--save` writes `results/<run>.md` (setup line, environment, tables), `results/<run>.csv` (the
aggregated table), `results/<run>_instances.csv` (one row per instance and, for `chi_break`, per chi and
rule) and, for `chi_break`, `results/<run>_curves.csv`. The functions in `chibreak.experiments` take
keyword arguments for `ns`, `alphas`, `seeds`, `max_chi` and `rules`:

```python
from chibreak import experiments as E
out = E.chi_break(ns=(12,), alphas=(0.0, 0.2), seeds=3, max_chi=64, rules=("top_chi",))
print(E.md_table(out["table"]))
```

Columns in the saved files:

- `depth_calibration`: `entropy_bits`, `entropy_std`, `percent_of_max`, `effective_rank`, `full_rank`.
- `alpha_sweep`: `peak_weight`, `peak_weight_std`, `chance_peak_weight`, `R`, `R_planted`,
  `peak_is_planted`, `min_abs_Z`, `median_abs_Z`, `frac_abs_Z_below_0.25`, `mid_entropy_bits`,
  `mid_effective_rank`.
- `chi_break` table: `chi_break_per_seed`, `chi_break_median`, `reached` (instances that hit R = 1 on the
  ladder). Curves: `R`, `R_std`, `fraction_R1`, `min_abs_Z`, `fidelity`, `discarded_weight` per chi.

## Tests

```bash
python -m pytest tests
python tests/run_tests.py
```

## Conventions

- Bit i is 1 when `<Z_i> < 0`.
- Qubit 0 is the most significant axis of the dense vector; bitstrings print qubit 0 first.
- RNG order inside `build_peaked_circuit` is block A, then s, then the random targets of block B, so a
  seed fixes A and s across alpha.
- Block A depth in the saved runs is `ceil(1.5 n)`; shallower blocks leave the balanced cut far from its
  maximum entropy (see the depth calibration in RESULTS.md).

## Limits

- Connectivity is a one-dimensional chain with nearest-neighbour gates. chi_break here describes this
  architecture and is not comparable to published all-to-all values.
- chi_break sits on a powers-of-two ladder. Differences below a factor of two show in R and in
  `min_abs_Z`, not in chi_break.
- Exact simulation needs 16 * 2^n bytes: 16 MB at n = 20, 256 MB at n = 24, 4 GB at n = 28. On a 16-core
  workstation one two-qubit gate takes 0.013 s at n = 20 and 0.56 s at n = 24, and the balanced-cut SVD
  takes 0.8 s at n = 20 and 33 s at n = 24.
- One MPS run costs about chi^3 per gate. At n = 20 with block A depth 30, chi = 256 takes about 80 s on
  the same machine.
- alpha moves the peak weight and the mirror structure at the same time; compare cells at matched
  `peak_weight` when that matters.
- The `random` rule keeps a uniform random subset and serves as a control. The randomized truncation of
  Harrow, Lowe and Witteveen is not implemented.
