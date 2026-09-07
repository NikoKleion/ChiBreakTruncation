# Tests

```bash
python -m pytest tests
python tests/run_tests.py       # no pytest needed
```

- `test_circuit.py` Haar gates are unitary, the geodesic hits both endpoints and stays unitary, block B is
  the reversed inverse of A at alpha = 0, alpha = 0 returns the planted string, and a seed fixes A and s
  across alpha.
- `test_statevector.py` two-qubit and one-qubit gate application against Kronecker products, the readout
  and peak bit order on a basis state, zero marginals on |+>, Schmidt spectrum normalisation.
- `test_metrics.py` entropy, effective rank, smallest marginal, fraction below a threshold and recovery on hand-computable
  inputs.
- `test_mps.py` the MPS at full bond dimension reproduces the exact state and marginals, the discarded
  weight of a truncation equals the dropped Schmidt weight, truncation rules return the right sizes, and a capped
  run stays normalised.
- `test_experiments.py` each run executes in quick mode and returns the documented columns.
- `test_pinned.py` numbers from the saved runs in `results/`, so a change in the generator or the MPS shows
  up against the published tables.
