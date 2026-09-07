# Results

Three runs, one command each, output in `results/`. Block A depth is `ceil(1.5 n)` except where depth is
the variable, five seeds per cell (three at n = 24), R against the exact true peak.

## depth_calibration

Entanglement entropy and effective rank at the balanced cut after block A, against block A depth, at
alpha = 0.

| n | max bits | depth n/2 | depth 1.5 n | depth 32 | Page value |
|--:|--:|--:|--:|--:|--:|
| 12 | 6 | 3.12 (52%) at 8 | 4.97 (83%) at 18 | 5.27 (88%) | 5.28 |
| 16 | 8 | 3.20 (40%) at 8 | 6.95 (87%) at 24 | 7.22 (90%) | 7.28 |
| 20 | 10 | 4.61 (46%) at 12 | 8.91 (89%) at 30 | 9.05 (90%) | 9.28 |
| 24 | 12 | 4.41 (37%) at 12 | beyond the run | 10.38 (86%) | 11.28 |

The depth-1.5n entries at n = 12, 16, 20 are the mid-circuit columns of the alpha sweep, which uses that
depth. The Page value, `n/2 - 1/(2 ln 2)` bits, is the balanced-cut entropy of a Haar-random state on n
qubits and the ceiling this brickwork approaches: at n = 12 the depth-32 value sits 0.01 bits under it.

Depth n/2 leaves the cut at 37 to 52 percent of its maximum, with an effective rank of 6 to 15
out of 64 to 4096. Depth 1.5 n reaches 83 to 89 percent and an effective rank of 24 of 64 (n = 12), 91
of 256 (n = 16) and 345 of 1024 (n = 20). The runs that follow use depth 1.5 n.

---

## alpha_sweep

Peak weight, recovery, smallest marginal and marginal spread against alpha at block A depth 1.5 n, five seeds
per cell. The `chance_peak_weight` column is the expected
largest of 2^n Porter-Thomas weights, `(ln 2^n + 0.5772) / 2^n`: the peak weight of a state with no peak.

| n | alpha | peak_weight | chance | R | min abs Z | fraction below 0.25 |
|--:|--:|--:|--:|--:|--:|--:|
| 12 | 0 | 1 | 2.17e-3 | 1 | 1 | 0 |
| 12 | 0.05 | 0.515 | | 1 | 0.567 | 0 |
| 12 | 0.075 | 0.225 | | 1 | 0.276 | 0 |
| 12 | 0.1 | 0.0710 | | 1 | 0.092 | 0.93 |
| 12 | 0.125 | 0.0161 | | 0.98 | 0.014 | 1 |
| 12 | 0.15 | 0.00286 | | 0.82 | 0.0009 | 1 |
| 12 | 0.2 | 0.00217 | | 0.60 | 0.0013 | 1 |
| 16 | 0 | 1 | 1.78e-4 | 1 | 1 | 0 |
| 16 | 0.05 | 0.284 | | 1 | 0.355 | 0 |
| 16 | 0.075 | 0.0592 | | 1 | 0.097 | 0.98 |
| 16 | 0.1 | 0.00676 | | 1 | 0.012 | 1 |
| 16 | 0.125 | 0.000476 | | 0.85 | 0.0010 | 1 |
| 16 | 0.15 | 0.000196 | | 0.54 | 0.0004 | 1 |
| 20 | 0 | 1 | 1.38e-5 | 1 | 1 | 0 |
| 20 | 0.05 | 0.128 | | 1 | 0.199 | 0.59 |
| 20 | 0.075 | 0.00984 | | 1 | 0.027 | 1 |
| 20 | 0.1 | 0.000268 | | 0.99 | 0.0008 | 1 |
| 20 | 0.125 | 1.36e-5 | | 0.51 | 0.0001 | 1 |
| 20 | 0.15 | 1.33e-5 | | 0.50 | 0.0001 | 1 |

- At alpha = 0 every instance returns the planted string with weight 1, as the construction requires.
- The peak weight falls faster with alpha at larger n: at alpha = 0.05 it is 0.52, 0.28 and 0.13 for
  n = 12, 16, 20. It reaches the chance value at alpha = 0.2 for n = 12, 0.15 for n = 16 and 0.125 for
  n = 20. From there on R sits between 0.38 and 0.60 and the argmax string differs from the planted one
  in every instance: the state carries no peak.
- R stays at 1 down to a peak weight of 0.071 (n = 12, alpha = 0.1), 0.0068 (n = 16, alpha = 0.1) and
  0.0098 (n = 20, alpha = 0.075), then falls to chance level within the next 0.05 of alpha.
- The smallest marginal falls before R does. At n = 16, alpha = 0.1 all five instances still read the
  peak (R = 1) while the smallest `|<Z_i>|` is 0.012 and every qubit is below 0.25. At n = 20,
  alpha = 0.075 R = 1 with a smallest marginal of 0.027 and every qubit below 0.25.
- At n = 20, alpha = 0.1 one instance of five reads one bit wrong from the exact state (seed 2,
  R = 0.95). Its peak weight is 2.1e-4, fifteen times chance, so a peak exists that the marginal signs
  no longer resolve. That instance cannot reach R = 1 at any bond dimension.
- The mid-circuit entropy at these depths is 4.97 of 6, 6.95 of 8 and 8.91 of 10 bits.

---

## chi_break

MPS simulation with the bond capped at chi on the ladder 2, 4, ..., min(2^(n/2), 256), under the
`top_chi` rule and the `random` control, five seeds per cell, block A depth 1.5 n. R is measured against
the exact true peak.

chi_break per seed, top_chi rule:

| n | full rank | alpha | peak_weight | chi_break per seed | median |
|--:|--:|--:|--:|--:|--:|
| 12 | 64 | 0 | 1 | 4 4 2 4 4 | 4 |
| 12 | 64 | 0.05 | 0.515 | 4 4 4 4 4 | 4 |
| 12 | 64 | 0.1 | 0.071 | 8 16 32 8 4 | 8 |
| 16 | 256 | 0 | 1 | 4 8 8 8 8 | 8 |
| 16 | 256 | 0.05 | 0.284 | 8 8 16 16 8 | 8 |
| 16 | 256 | 0.1 | 0.0068 | 32 8 64 32 64 | 32 |
| 20 | 1024 | 0 | 1 | 16 8 16 16 8 | 16 |
| 20 | 1024 | 0.05 | 0.128 | 16 8 32 16 32 | 16 |
| 20 | 1024 | 0.1 | 0.00027 | 256 256 >256 64 >256 | 3 of 5 reached |

Under the `random` rule every instance at every n and alpha reaches R = 1 only at the full rank, where no
truncation happens.

Curves at n = 16, alpha = 0, averaged over five seeds:

| chi | top_chi R | fraction R = 1 | min abs Z | fidelity | discarded | random R | random fidelity | random discarded |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2 | 0.55 | 0 | 0.011 | 0.0002 | 21.2 | 0.49 | 3e-6 | 153.9 |
| 4 | 0.83 | 0.2 | 0.032 | 0.007 | 12.9 | 0.53 | 3e-6 | 128.9 |
| 8 | 1 | 1 | 0.112 | 0.072 | 7.46 | 0.53 | 7e-6 | 99.7 |
| 16 | 1 | 1 | 0.275 | 0.210 | 4.56 | 0.53 | 6e-6 | 77.7 |
| 32 | 1 | 1 | 0.476 | 0.426 | 2.59 | 0.60 | 3e-5 | 52.0 |
| 64 | 1 | 1 | 0.728 | 0.703 | 1.12 | 0.58 | 2e-5 | 32.8 |
| 128 | 1 | 1 | 0.934 | 0.932 | 0.21 | 0.49 | 8e-5 | 10.3 |
| 256 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 0 |

- Under top-chi truncation at alpha = 0 and 0.05 the median chi_break is 4, 8 and 16 for n = 12, 16 and
  20, with every seed within a factor of two of the median. The full rank is 64, 256 and 1024. Over
  these sizes chi_break doubles per four qubits while the full rank quadruples.
- At alpha = 0.1, where the peak weight is 0.071, 0.0068 and 0.00027, the medians rise to 8, 32 and at
  least 256. At n = 20 two instances do not reach R = 1 on the ladder: seed 2 misreads one bit from the
  exact state, so no bond dimension recovers it, and seed 4 needs chi above 256.
- Recovery happens long before the state is accurate. At chi_break the fidelity to the exact state is
  0.14 (n = 12) and 0.072 (n = 16) at alpha = 0, and 0.10 at n = 16, alpha = 0.05. Fidelity above 0.9
  needs chi = 32 at n = 12 and chi = 128 at n = 16, eight to sixteen times chi_break.
- At alpha = 0 the smallest marginal of the truncated state tracks its fidelity: at n = 16 the pairs (min
  abs Z, fidelity) run 0.11 and 0.072 at chi = 8, 0.48 and 0.43 at chi = 32, 0.93 and 0.93 at
  chi = 128. The exact state has smallest marginal 1. A truncated state with fidelity F to |s> has `<Z_i>` close to
  F times the sign, so the smallest marginal reads how much of the peak survived.
- The random control reads at chance, R between 0.40 and 0.61, at every chi below the full rank, with
  fidelity between 1e-6 and 1e-3. At chi equal to half the full rank it discards 8.1 units of Schmidt
  weight over the run at n = 12 and 10.3 at n = 16, against 0.20 and 0.21 for top-chi. At these depths
  the spectrum still has a tail that carries the peak, and which components a rule keeps decides the
  outcome.
