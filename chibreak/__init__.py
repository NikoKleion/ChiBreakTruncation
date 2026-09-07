# chibreak: peaked circuits on a 1D brickwork, exact ground truth, and MPS truncation with a pluggable truncation rule
from .circuit import build_peaked_circuit, brickwork_pairs, haar_2q, geodesic
from .statevector import (zero_state, apply_2q, apply_1q, simulate_exact, state_after_block_a,
                          z_marginals, bits_from_marginals, true_peak, schmidt_spectrum)
from .metrics import entropy_bits, effective_rank, flatness, min_abs_z, fraction_below, recovery
from .mps import MPS, TRUNCATION_RULES, simulate_mps, chi_ladder

__all__ = ["build_peaked_circuit", "brickwork_pairs", "haar_2q", "geodesic",
           "zero_state", "apply_2q", "apply_1q", "simulate_exact", "state_after_block_a",
           "z_marginals", "bits_from_marginals", "true_peak", "schmidt_spectrum",
           "entropy_bits", "effective_rank", "flatness", "min_abs_z", "fraction_below", "recovery",
           "MPS", "TRUNCATION_RULES", "simulate_mps", "chi_ladder"]
