# chibreak.statevector: dense simulation, Z marginals, peak readout, Schmidt spectrum
import numpy as np

X = np.array([[0, 1], [1, 0]], dtype=np.complex128)


def zero_state(n):
    v = np.zeros(2 ** n, dtype=np.complex128)
    v[0] = 1.0
    return v


def apply_2q(state, u, i, j, n):
    # 4x4 gate on qubits (i, j); qubit 0 is the most significant axis
    t = np.tensordot(u.reshape(2, 2, 2, 2), state.reshape((2,) * n), axes=([2, 3], [i, j]))
    perm = [None] * n
    perm[i], perm[j] = 0, 1
    for slot, ax in zip([a for a in range(n) if a not in (i, j)], range(2, n)):
        perm[slot] = ax
    return np.transpose(t, perm).reshape(-1)


def apply_1q(state, u, i, n):
    t = np.tensordot(u, state.reshape((2,) * n), axes=([1], [i]))
    perm = [None] * n
    perm[i] = 0
    for slot, ax in zip([a for a in range(n) if a != i], range(1, n)):
        perm[slot] = ax
    return np.transpose(t, perm).reshape(-1)


def simulate_exact(circuit):
    # A, B, then the X^s layer on |0...0>
    n = circuit["n"]
    state = zero_state(n)
    for (u, (i, j)) in circuit["gates"]:
        state = apply_2q(state, u, i, j, n)
    for i, bit in enumerate(circuit["x_layer"]):
        if bit:
            state = apply_1q(state, X, i, n)
    return state


def state_after_block_a(circuit):
    n = circuit["n"]
    state = zero_state(n)
    for (u, (i, j)) in circuit["gates"][: circuit["a_gate_count"]]:
        state = apply_2q(state, u, i, j, n)
    return state


def z_marginals(state, n):
    # <Z_i> = P(q_i = 0) - P(q_i = 1)
    p = np.abs(state.reshape((2,) * n)) ** 2
    out = np.empty(n)
    for i in range(n):
        pi = p.sum(axis=tuple(a for a in range(n) if a != i))
        out[i] = pi[0] - pi[1]
    return out


def bits_from_marginals(z):
    # bit_i = 1 iff <Z_i> < 0
    return (np.asarray(z) < 0).astype(int)


def true_peak(state, n):
    # argmax |amplitude|^2 as a bitstring (qubit 0 first) and its probability
    probs = np.abs(state) ** 2
    idx = int(np.argmax(probs))
    bits = np.array([(idx >> (n - 1 - k)) & 1 for k in range(n)], dtype=int)
    return bits, float(probs[idx])


def schmidt_spectrum(state, cut, n):
    # singular values across the bond between qubits [0, cut) and [cut, n)
    return np.linalg.svd(state.reshape(2 ** cut, 2 ** (n - cut)), compute_uv=False)
