# the runs execute end to end in quick mode and return the documented columns
from chibreak import experiments as E


def test_depth_calibration_quick():
    out = E.depth_calibration(quick=True)
    assert {"n", "depth", "entropy_bits", "percent_of_max", "effective_rank"} <= set(out["table"][0])
    assert all(0 < r["percent_of_max"] <= 100 for r in out["table"])


def test_alpha_sweep_quick():
    out = E.alpha_sweep(quick=True)
    rows = out["table"]
    assert {"alpha", "peak_weight", "R", "min_abs_Z", "frac_abs_Z_below_0.25", "mid_entropy_bits"} <= set(rows[0])
    assert rows[0]["alpha"] == 0.0 and rows[0]["peak_weight"] > 1 - 1e-9 and rows[0]["R"] == 1.0


def test_chi_break_quick():
    out = E.chi_break(quick=True)
    assert {"chi_break_per_seed", "chi_break_median", "reached"} <= set(out["table"][0])
    assert {"chi", "R", "fidelity", "min_abs_Z"} <= set(out["curves"][0])
    top = [r for r in out["curves"] if r["rule"] == "top_chi" and r["alpha"] == 0.0]
    assert top[-1]["R"] == 1.0 and top[-1]["fidelity"] > 1 - 1e-9


def test_block_a_depth_rule():
    assert E.block_a_depth(12) == 18 and E.block_a_depth(14) == 21 and E.block_a_depth(20) == 30
