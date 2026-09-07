# spectrum and readout quantities
import numpy as np

from chibreak import entropy_bits, effective_rank, flatness, min_abs_z, fraction_below, recovery


def test_flat_spectrum_entropy_and_rank():
    sv = np.ones(16) / 4.0
    assert abs(entropy_bits(sv) - 4.0) < 1e-12
    assert abs(effective_rank(sv) - 16.0) < 1e-9 and abs(flatness(sv) - 1.0) < 1e-9


def test_single_coefficient_spectrum():
    sv = np.array([1.0, 0.0, 0.0, 0.0])
    assert entropy_bits(sv) == 0.0 and abs(effective_rank(sv) - 1.0) < 1e-12


def test_margin_fraction_below_recovery():
    z = np.array([0.9, -0.1, 0.3, -0.8])
    assert min_abs_z(z) == 0.1 and fraction_below(z) == 0.25
    assert recovery([0, 1, 0, 1], [0, 1, 1, 1]) == 0.75
