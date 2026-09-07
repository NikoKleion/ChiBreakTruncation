# chibreak.mps: matrix product state simulation of the circuit with a pluggable truncation rule
import numpy as np

from .statevector import X

Z = np.diag([1.0, -1.0]).astype(np.complex128)


# truncation rules: singular values arrive sorted descending; return the indices to keep
def top_chi(s, chi, rng=None):
    return np.arange(min(chi, s.size))


def random_subset(s, chi, rng):
    # uniform random subset of chi components, a control
    return np.sort(rng.choice(s.size, size=min(chi, s.size), replace=False))


TRUNCATION_RULES = {"top_chi": top_chi, "random": random_subset}


def chi_ladder(n, max_chi=None):
    # powers of two up to the exact bond dimension 2^(n/2), capped at max_chi
    chis = [2 ** k for k in range(1, n // 2 + 1)]
    return [c for c in chis if max_chi is None or c <= max_chi]


class MPS:
    # tensors A[k] of shape (left, 2, right); sites left of `center` are left-canonical, sites right
    # of it right-canonical, so a truncation at the center bond is the Eckart-Young optimum
    def __init__(self, n):
        self.n = n
        self.A = []
        for _ in range(n):
            a = np.zeros((1, 2, 1), dtype=np.complex128)
            a[0, 0, 0] = 1.0
            self.A.append(a)
        self.center = 0
        self.discarded = 0.0

    def _shift_right(self, k):
        l, d, r = self.A[k].shape
        q, rr = np.linalg.qr(self.A[k].reshape(l * d, r))
        self.A[k] = q.reshape(l, d, q.shape[1])
        self.A[k + 1] = np.tensordot(rr, self.A[k + 1], axes=([1], [0]))
        self.center = k + 1

    def _shift_left(self, k):
        l, d, r = self.A[k].shape
        q, rr = np.linalg.qr(self.A[k].reshape(l, d * r).T)
        self.A[k] = q.T.reshape(q.shape[1], d, r)
        self.A[k - 1] = np.tensordot(self.A[k - 1], rr.T, axes=([2], [0]))
        self.center = k - 1

    def move_center(self, k):
        while self.center < k:
            self._shift_right(self.center)
        while self.center > k:
            self._shift_left(self.center)

    def apply_2q(self, u, i, chi, rule=top_chi, rng=None):
        # gate on (i, i+1), SVD at that bond, the rule picks the components, bond capped at chi
        if self.center < i:
            self.move_center(i)
        elif self.center > i + 1:
            self.move_center(i + 1)
        theta = np.tensordot(self.A[i], self.A[i + 1], axes=([2], [0]))
        l, r = theta.shape[0], theta.shape[3]
        theta = np.tensordot(u.reshape(2, 2, 2, 2), theta, axes=([2, 3], [1, 2]))
        theta = theta.transpose(2, 0, 1, 3).reshape(l * 2, 2 * r)
        U, s, Vh = np.linalg.svd(theta, full_matrices=False)
        idx = np.asarray(rule(s, chi, rng))
        kept = s[idx]
        w = float(np.dot(kept, kept))
        self.discarded += float(np.dot(s, s)) - w
        kept = kept / np.sqrt(w)
        self.A[i] = U[:, idx].reshape(l, 2, idx.size)
        self.A[i + 1] = (kept[:, None] * Vh[idx]).reshape(idx.size, 2, r)
        self.center = i + 1

    def apply_1q(self, u, i):
        self.A[i] = np.tensordot(u, self.A[i], axes=([1], [1])).transpose(1, 0, 2)

    def bond_dimensions(self):
        return [a.shape[2] for a in self.A[:-1]]

    def norm(self):
        e = np.ones((1, 1), dtype=np.complex128)
        for a in self.A:
            e = np.einsum("ab,aic,bid->cd", e, a.conj(), a, optimize=True)
        return float(np.sqrt(e[0, 0].real))

    def z_marginals(self):
        # <Z_i> for every qubit through left and right environments
        n = self.n
        L = [np.ones((1, 1), dtype=np.complex128)]
        for a in self.A[:-1]:
            L.append(np.einsum("ab,aic,bid->cd", L[-1], a.conj(), a, optimize=True))
        R = [None] * n
        R[n - 1] = np.ones((1, 1), dtype=np.complex128)
        for k in range(n - 1, 0, -1):
            a = self.A[k]
            R[k - 1] = np.einsum("aic,bid,cd->ab", a.conj(), a, R[k], optimize=True)
        norm = np.einsum("ab,aic,bid,cd->", L[0], self.A[0].conj(), self.A[0], R[0], optimize=True).real
        out = np.empty(n)
        for k in range(n):
            a = self.A[k]
            out[k] = np.einsum("ab,aic,ij,bjd,cd->", L[k], a.conj(), Z, a, R[k], optimize=True).real / norm
        return out

    def to_dense(self):
        # full state vector, qubit 0 most significant; for small n only
        v = self.A[0]
        for a in self.A[1:]:
            v = np.tensordot(v, a, axes=([v.ndim - 1], [0]))
        return v.reshape(-1)


def simulate_mps(circuit, chi, rule="top_chi", seed=0):
    # run A, B, then X^s with the bond capped at chi under the named truncation rule
    fn = TRUNCATION_RULES[rule]
    rng = np.random.default_rng(seed)
    m = MPS(circuit["n"])
    for (u, (i, j)) in circuit["gates"]:
        m.apply_2q(u, i, chi, fn, rng)
    for i, bit in enumerate(circuit["x_layer"]):
        if bit:
            m.apply_1q(X, i)
    return m


def fidelity(dense, mps):
    # |<dense|mps>|^2 with both states normalised
    v = mps.to_dense()
    return float(abs(np.vdot(dense, v)) ** 2 / (np.vdot(dense, dense).real * np.vdot(v, v).real))
