# chibreak.metrics: spectrum and readout diagnostics
import numpy as np


def _probs(sv):
    p = np.asarray(sv, dtype=float) ** 2
    return p / p.sum()


def entropy_bits(sv):
    p = _probs(sv)
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log2(p)))


def effective_rank(sv):
    # participation ratio 1 / sum p_i^2
    p = _probs(sv)
    return float(1.0 / np.sum(p ** 2))


def flatness(sv):
    # effective rank over the number of singular values; 1 for a flat spectrum
    return effective_rank(sv) / len(sv)


def min_abs_z(z):
    # smallest |<Z_i>| over qubits
    return float(np.min(np.abs(z)))


def fraction_below(z, threshold=0.25):
    # fraction of qubits with |<Z_i>| below the threshold
    z = np.abs(np.asarray(z))
    return float(np.mean(z < threshold))


def recovery(bits, target):
    # fraction of bits equal to the target string
    return float(np.mean(np.asarray(bits) == np.asarray(target)))
