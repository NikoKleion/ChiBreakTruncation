# dense simulation against Kronecker references; readout conventions
import numpy as np

from chibreak import apply_2q, apply_1q, haar_2q, z_marginals, bits_from_marginals, true_peak, schmidt_spectrum
from chibreak.statevector import X


def _random_state(n, rng):
    v = rng.standard_normal(2 ** n) + 1j * rng.standard_normal(2 ** n)
    return v / np.linalg.norm(v)


def test_apply_2q_matches_kronecker_reference():
    rng = np.random.default_rng(0)
    n = 5
    psi = _random_state(n, rng)
    for i in range(n - 1):
        u = haar_2q(rng)
        ref = np.kron(np.kron(np.eye(2 ** i), u), np.eye(2 ** (n - i - 2))) @ psi
        assert np.allclose(apply_2q(psi, u, i, i + 1, n), ref, atol=1e-12)


def test_apply_1q_matches_kronecker_reference():
    rng = np.random.default_rng(1)
    n = 5
    psi = _random_state(n, rng)
    for i in range(n):
        ref = np.kron(np.kron(np.eye(2 ** i), X), np.eye(2 ** (n - i - 1))) @ psi
        assert np.allclose(apply_1q(psi, X, i, n), ref, atol=1e-12)


def test_marginals_readout_and_peak_order_on_basis_state():
    n = 5
    e = np.zeros(2 ** n, dtype=complex)
    e[int("01101", 2)] = 1.0
    z = z_marginals(e, n)
    assert list(bits_from_marginals(z)) == [0, 1, 1, 0, 1]
    bits, w = true_peak(e, n)
    assert list(bits) == [0, 1, 1, 0, 1] and w == 1.0


def test_marginals_of_plus_state_are_zero():
    n = 3
    psi = np.ones(2 ** n, dtype=complex) / np.sqrt(2 ** n)
    assert np.allclose(z_marginals(psi, n), 0.0)


def test_schmidt_spectrum_normalised_and_product_state_rank_one():
    rng = np.random.default_rng(2)
    psi = _random_state(6, rng)
    sv = schmidt_spectrum(psi, 3, 6)
    assert abs((sv ** 2).sum() - 1.0) < 1e-12 and len(sv) == 8
    e = np.zeros(64, dtype=complex)
    e[5] = 1.0
    sv = schmidt_spectrum(e, 3, 6)
    assert abs(sv[0] - 1.0) < 1e-12 and np.allclose(sv[1:], 0.0)
