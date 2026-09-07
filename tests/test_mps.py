# MPS engine: exactness at full bond dimension, canonical form, truncation accounting, truncation rules
import numpy as np

from chibreak import (build_peaked_circuit, simulate_exact, simulate_mps, z_marginals, bits_from_marginals,
                      schmidt_spectrum, MPS, TRUNCATION_RULES, chi_ladder)
from chibreak.mps import fidelity, random_subset


def test_full_chi_mps_matches_exact_state():
    c = build_peaked_circuit(8, 12, 0.1, seed=0)
    exact = simulate_exact(c)
    m = simulate_mps(c, chi=16)
    assert abs(m.norm() - 1.0) < 1e-10
    assert fidelity(exact, m) > 1 - 1e-10
    assert np.allclose(m.z_marginals(), z_marginals(exact, 8), atol=1e-9)
    assert abs(m.discarded) < 1e-10


def test_to_dense_bit_order():
    m = MPS(3)
    m.apply_1q(np.array([[0, 1], [1, 0]], dtype=complex), 1)
    v = m.to_dense()
    assert abs(v[int("010", 2)] - 1.0) < 1e-12


def test_truncation_discards_smallest_schmidt_weight():
    # build an entangled 4-qubit state at full chi, then an identity gate at chi = 1 on the middle bond
    c = build_peaked_circuit(4, 3, 0.0, seed=5)
    m = MPS(4)
    for (u, (i, j)) in c["gates"][: c["a_gate_count"]]:
        m.apply_2q(u, i, chi=4)
    sv = schmidt_spectrum(m.to_dense(), 2, 4)
    m.apply_2q(np.eye(4, dtype=complex), 1, chi=1)
    assert abs(m.discarded - float((sv[1:] ** 2).sum())) < 1e-10
    assert m.bond_dimensions()[1] == 1 and abs(m.norm() - 1.0) < 1e-10


def test_truncation_rules_sizes_and_ladder():
    s = np.sort(np.random.default_rng(0).random(10))[::-1]
    assert list(TRUNCATION_RULES["top_chi"](s, 4)) == [0, 1, 2, 3]
    idx = random_subset(s, 4, np.random.default_rng(1))
    assert len(idx) == 4 and len(set(idx)) == 4
    assert chi_ladder(12) == [2, 4, 8, 16, 32, 64] and chi_ladder(20, 256)[-1] == 256


def test_alpha_zero_recovers_at_full_chi_and_capped_run_stays_normalised():
    c = build_peaked_circuit(10, 15, 0.0, seed=1)
    m = simulate_mps(c, chi=32)
    assert (bits_from_marginals(m.z_marginals()) == c["peak"]).all()
    m2 = simulate_mps(c, chi=2)
    assert all(b <= 2 for b in m2.bond_dimensions()) and abs(m2.norm() - 1.0) < 1e-10


def test_random_rule_at_full_chi_is_exact():
    c = build_peaked_circuit(6, 9, 0.2, seed=7)
    m = simulate_mps(c, chi=8, rule="random", seed=3)
    assert fidelity(simulate_exact(c), m) > 1 - 1e-10
