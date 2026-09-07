# numbers pinned from the saved runs in results/, recomputed on the cheap cells (n <= 16)
from chibreak import experiments as E


def _row(rows, **key):
    out = [r for r in rows if all(r[k] == v for k, v in key.items())]
    assert len(out) == 1, key
    return out[0]


def test_depth_calibration_pins():
    # results/depth_calibration.md: n = 12 depth 20 and n = 16 depth 24, five seeds each
    out = E.depth_calibration(ns=(12, 16), seeds=(5, 5), depths=(20, 24))
    r12 = _row(out["table"], n=12, depth=20)
    r16 = _row(out["table"], n=16, depth=24)
    assert abs(r12["entropy_bits"] - 5.09445) < 1e-4 and abs(r12["effective_rank"] - 26.5738) < 1e-3
    assert abs(r16["entropy_bits"] - 6.95071) < 1e-4 and abs(r16["effective_rank"] - 91.1005) < 1e-3


def test_alpha_sweep_pins():
    # results/alpha_sweep.md: n = 12 at depth 18 and n = 16 at depth 24, five seeds
    out = E.alpha_sweep(ns=(12, 16), alphas=(0.0, 0.05, 0.1, 0.15), seeds=5)
    r = _row(out["table"], n=12, alpha=0.05)
    assert abs(r["peak_weight"] - 0.51453) < 1e-4 and r["R"] == 1.0 and abs(r["min_abs_Z"] - 0.567108) < 1e-4
    r = _row(out["table"], n=12, alpha=0.1)
    assert r["R"] == 1.0 and abs(r["min_abs_Z"] - 0.0917777) < 1e-4 and abs(r["frac_abs_Z_below_0.25"] - 0.933333) < 1e-4
    r = _row(out["table"], n=12, alpha=0.15)
    assert abs(r["R"] - 0.816667) < 1e-4 and abs(r["peak_weight"] - 0.00286492) < 1e-6
    assert abs(r["chance_peak_weight"] - 0.0021716) < 1e-6
    r = _row(out["table"], n=16, alpha=0.05)
    assert abs(r["peak_weight"] - 0.284282) < 1e-4 and r["R"] == 1.0 and abs(r["min_abs_Z"] - 0.355384) < 1e-4
    r = _row(out["table"], n=16, alpha=0.1)
    assert r["R"] == 1.0 and abs(r["min_abs_Z"] - 0.012425) < 1e-5 and r["frac_abs_Z_below_0.25"] == 1.0
    r = _row(out["table"], n=16, alpha=0.15)
    assert abs(r["R"] - 0.5375) < 1e-4 and abs(r["peak_weight"] - 0.00019568) < 1e-7
    assert abs(r["mid_entropy_bits"] - 6.95071) < 1e-4


def test_chi_break_pins():
    # results/chi_break.md: n = 12 at depth 18, both rules, five seeds
    out = E.chi_break(ns=(12,), alphas=(0.0, 0.1), seeds=5, max_chi=64)
    assert _row(out["table"], alpha=0.0, rule="top_chi")["chi_break_per_seed"] == "4 4 2 4 4"
    assert _row(out["table"], alpha=0.1, rule="top_chi")["chi_break_per_seed"] == "8 16 32 8 4"
    assert _row(out["table"], alpha=0.0, rule="random")["chi_break_per_seed"] == "64 64 64 64 64"
    c = _row(out["curves"], alpha=0.0, rule="top_chi", chi=4)
    assert c["R"] == 1.0 and abs(c["fidelity"] - 0.143393) < 1e-4 and abs(c["min_abs_Z"] - 0.192338) < 1e-4
    c = _row(out["curves"], alpha=0.0, rule="top_chi", chi=32)
    assert abs(c["fidelity"] - 0.929865) < 1e-4 and abs(c["discarded_weight"] - 0.20124) < 1e-3
    c = _row(out["curves"], alpha=0.0, rule="random", chi=32)
    assert abs(c["R"] - 0.6) < 1e-6 and abs(c["discarded_weight"] - 8.07729) < 1e-3
