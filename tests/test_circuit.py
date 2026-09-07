# circuit generator: Haar gates, geodesic, mirror structure, seed behaviour
import numpy as np

from chibreak import build_peaked_circuit, brickwork_pairs, haar_2q, geodesic, simulate_exact, true_peak


def test_haar_gate_is_unitary():
    rng = np.random.default_rng(0)
    for _ in range(50):
        u = haar_2q(rng)
        assert np.allclose(u.conj().T @ u, np.eye(4), atol=1e-12)


def test_geodesic_endpoints_and_unitarity():
    rng = np.random.default_rng(1)
    a, b = haar_2q(rng), haar_2q(rng)
    assert np.allclose(geodesic(a, b, 0.0), a) and np.allclose(geodesic(a, b, 1.0), b)
    for alpha in (0.1, 0.5, 0.9):
        g = geodesic(a, b, alpha)
        assert np.allclose(g.conj().T @ g, np.eye(4), atol=1e-10)


def test_brickwork_pairs_alternate():
    assert brickwork_pairs(6, 0) == [(0, 1), (2, 3), (4, 5)]
    assert brickwork_pairs(6, 1) == [(1, 2), (3, 4)]
    assert brickwork_pairs(5, 0) == [(0, 1), (2, 3)]


def test_block_b_is_reversed_inverse_at_alpha_zero():
    c = build_peaked_circuit(6, 5, 0.0, seed=2)
    k = c["a_gate_count"]
    assert k == sum(len(brickwork_pairs(6, l)) for l in range(5)) and len(c["gates"]) == 2 * k
    for (u, q), (v, p) in zip(c["gates"][:k], reversed(c["gates"][k:])):
        assert q == p and np.allclose(v, u.conj().T)


def test_alpha_zero_returns_planted_string():
    c = build_peaked_circuit(8, 12, 0.0, seed=3)
    bits, w = true_peak(simulate_exact(c), 8)
    assert w > 1 - 1e-10 and (bits == c["peak"]).all()


def test_seed_fixes_block_a_and_planted_string_across_alpha():
    c0, c1 = build_peaked_circuit(8, 12, 0.0, seed=4), build_peaked_circuit(8, 12, 0.3, seed=4)
    k = c0["a_gate_count"]
    assert all(np.allclose(a[0], b[0]) for a, b in zip(c0["gates"][:k], c1["gates"][:k]))
    assert (c0["peak"] == c1["peak"]).all()
    assert not np.allclose(c0["gates"][k][0], c1["gates"][k][0])
