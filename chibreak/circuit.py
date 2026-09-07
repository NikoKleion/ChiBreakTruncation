# chibreak.circuit: mirrored peaked circuit C = X^s B(alpha) A on a 1D brickwork
import numpy as np
from scipy.linalg import expm, logm


def haar_2q(rng):
    # Haar-random 4x4 unitary: QR of a complex Ginibre matrix with the phase fix
    z = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / np.sqrt(2.0)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r)
    return q * (d / np.abs(d))


def brickwork_pairs(n, layer):
    # nearest-neighbour pairs of one layer; even layers start at qubit 0, odd layers at qubit 1
    return [(i, i + 1) for i in range(layer % 2, n - 1, 2)]


def geodesic(u_from, u_to, alpha):
    # U(alpha) = u_from expm(alpha logm(u_from^dagger u_to)); u_from at alpha 0, u_to at alpha 1
    if alpha <= 0.0:
        return u_from
    if alpha >= 1.0:
        return u_to
    return u_from @ expm(alpha * logm(u_from.conj().T @ u_to))


def build_peaked_circuit(n, depth, alpha=0.0, seed=None):
    # block A: depth brickwork layers of Haar gates; planted string s; block B: A reversed, each
    # inverse moved toward a fresh Haar gate by alpha. RNG order is A, s, B, so a seed fixes A and s
    # at every alpha.
    rng = np.random.default_rng(seed)
    a_gates = []
    for layer in range(depth):
        for (i, j) in brickwork_pairs(n, layer):
            a_gates.append((haar_2q(rng), (i, j)))
    s = rng.integers(0, 2, size=n)
    b_gates = []
    for (u, (i, j)) in reversed(a_gates):
        u_inv = u.conj().T
        b_gates.append((geodesic(u_inv, haar_2q(rng), alpha) if alpha > 0.0 else u_inv, (i, j)))
    return {"gates": a_gates + b_gates, "x_layer": s.copy(), "peak": s, "n": n, "depth": depth,
            "alpha": alpha, "seed": seed, "a_gate_count": len(a_gates)}
