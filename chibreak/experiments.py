# chibreak.experiments: the documented runs; python -m chibreak <run> [--save] writes results/<run>.{csv,md}
import csv
import math
import os
import platform
import sys
import time
from datetime import date

import numpy as np

from .circuit import build_peaked_circuit, brickwork_pairs
from .statevector import (zero_state, apply_2q, simulate_exact, state_after_block_a, z_marginals,
                          bits_from_marginals, true_peak, schmidt_spectrum)
from .metrics import entropy_bits, effective_rank, min_abs_z, fraction_below, recovery
from .mps import simulate_mps, chi_ladder, fidelity

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
ALPHAS = (0.0, 0.02, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.5, 1.0)
EULER_GAMMA = 0.5772156649015329


def chance_peak_weight(n):
    # expected largest of 2^n Porter-Thomas weights, (ln 2^n + gamma) / 2^n: the peak weight of a state with no peak
    N = 2.0 ** n
    return (math.log(N) + EULER_GAMMA) / N


def block_a_depth(n):
    # depth rule: block A at least 1.5 n layers
    return int(math.ceil(1.5 * n))


def _mean(rows, key):
    return float(np.mean([r[key] for r in rows]))


def _std(rows, key):
    return float(np.std([r[key] for r in rows]))


# depth calibration: balanced-cut entropy and effective rank after block A, read at checkpoints of one deep build
def depth_calibration(ns=(12, 16, 20, 24), seeds=(5, 5, 5, 3), depths=tuple(range(4, 33, 4)), quick=False):
    if quick:
        ns, seeds, depths = (8,), (2,), (4, 8)
    inst, agg = [], []
    for n, n_seeds in zip(ns, seeds):
        for seed in range(n_seeds):
            c = build_peaked_circuit(n, max(depths), 0.0, seed)
            psi, g = zero_state(n), 0
            for layer in range(max(depths)):
                for (i, j) in brickwork_pairs(n, layer):
                    psi = apply_2q(psi, c["gates"][g][0], i, j, n)
                    g += 1
                if (layer + 1) in depths:
                    sv = schmidt_spectrum(psi, n // 2, n)
                    inst.append({"n": n, "depth": layer + 1, "seed": seed,
                                 "entropy_bits": entropy_bits(sv), "effective_rank": effective_rank(sv)})
        for d in depths:
            rows = [r for r in inst if r["n"] == n and r["depth"] == d]
            agg.append({"n": n, "depth": d, "seeds": len(rows), "max_entropy_bits": n // 2,
                        "entropy_bits": _mean(rows, "entropy_bits"), "entropy_std": _std(rows, "entropy_bits"),
                        "percent_of_max": 100.0 * _mean(rows, "entropy_bits") / (n // 2),
                        "effective_rank": _mean(rows, "effective_rank"), "full_rank": 2 ** (n // 2)})
    setup = (f"alpha = 0; block A built once at depth {max(depths)} per seed and read at each checkpoint depth, "
             f"so the depth-d row is the depth-d circuit of that seed; balanced cut n/2 after block A; "
             f"seeds per n: {dict(zip(ns, seeds))}.")
    return {"table": agg, "instances": inst, "setup": setup}


# alpha sweep: peak weight, recovery, smallest marginal and marginal spread at block A depth 1.5 n
def alpha_sweep(ns=(12, 16, 20), alphas=ALPHAS, seeds=5, quick=False):
    if quick:
        ns, alphas, seeds = (8,), (0.0, 0.1, 1.0), 2
    inst, agg = [], []
    for n in ns:
        depth = block_a_depth(n)
        mid = {}
        for alpha in alphas:
            for seed in range(seeds):
                c = build_peaked_circuit(n, depth, alpha, seed)
                if seed not in mid:
                    sv = schmidt_spectrum(state_after_block_a(c), n // 2, n)
                    mid[seed] = (entropy_bits(sv), effective_rank(sv))
                psi = simulate_exact(c)
                z = z_marginals(psi, n)
                bits = bits_from_marginals(z)
                peak_bits, w = true_peak(psi, n)
                inst.append({"n": n, "depth": depth, "alpha": alpha, "seed": seed, "peak_weight": w,
                             "R": recovery(bits, peak_bits), "R_planted": recovery(bits, c["peak"]),
                             "peak_is_planted": int(bool((peak_bits == c["peak"]).all())),
                             "min_abs_Z": min_abs_z(z), "median_abs_Z": float(np.median(np.abs(z))),
                             "frac_abs_Z_below_0.25": fraction_below(z), "mid_entropy_bits": mid[seed][0],
                             "mid_effective_rank": mid[seed][1]})
        for alpha in alphas:
            rows = [r for r in inst if r["n"] == n and r["alpha"] == alpha]
            agg.append({"n": n, "depth": depth, "alpha": alpha, "seeds": len(rows),
                        "peak_weight": _mean(rows, "peak_weight"), "peak_weight_std": _std(rows, "peak_weight"),
                        "chance_peak_weight": chance_peak_weight(n),
                        "R": _mean(rows, "R"), "R_planted": _mean(rows, "R_planted"),
                        "peak_is_planted": _mean(rows, "peak_is_planted"),
                        "min_abs_Z": _mean(rows, "min_abs_Z"), "median_abs_Z": _mean(rows, "median_abs_Z"),
                        "frac_abs_Z_below_0.25": _mean(rows, "frac_abs_Z_below_0.25"),
                        "mid_entropy_bits": _mean(rows, "mid_entropy_bits"), "max_entropy_bits": n // 2,
                        "mid_effective_rank": _mean(rows, "mid_effective_rank"), "full_rank": 2 ** (n // 2)})
    depths = ", ".join(f"n={n} depth={block_a_depth(n)}" for n in ns)
    setup = (f"block A depth = ceil(1.5 n): {depths}; {seeds} seeds per (n, alpha); a seed fixes block A and the "
             f"planted string across alpha; R is against the exact true peak, R_planted against the planted "
             f"string; frac_abs_Z_below_0.25 is the fraction of qubits with |<Z_i>| < 0.25; chance_peak_weight is the expected largest of 2^n "
             f"Porter-Thomas weights; mid-circuit columns are the balanced cut after block A.")
    return {"table": agg, "instances": inst, "setup": setup}


# chi_break: MPS at capped bond dimension under each truncation rule; R, smallest marginal, fidelity against chi
def chi_break(ns=(12, 16, 20), alphas=(0.0, 0.05, 0.1), seeds=5, max_chi=256, rules=("top_chi", "random"),
              random_max_n=16, fidelity_max_n=16, quick=False):
    if quick:
        ns, alphas, seeds, max_chi = (8,), (0.0, 0.1), 2, 16
    inst, curves, summary = [], [], []
    for n in ns:
        depth = block_a_depth(n)
        chis = chi_ladder(n, max_chi)
        for alpha in alphas:
            for seed in range(seeds):
                c = build_peaked_circuit(n, depth, alpha, seed)
                exact = simulate_exact(c)
                peak_bits, w = true_peak(exact, n)
                for rule in rules:
                    if rule == "random" and n > random_max_n:
                        continue
                    first = None
                    for chi in chis:
                        t0 = time.perf_counter()
                        m = simulate_mps(c, chi, rule, seed)
                        z = m.z_marginals()
                        R = recovery(bits_from_marginals(z), peak_bits)
                        if first is None and R == 1.0:
                            first = chi
                        inst.append({"n": n, "depth": depth, "alpha": alpha, "seed": seed, "rule": rule,
                                     "chi": chi, "peak_weight": w, "R": R, "min_abs_Z": min_abs_z(z),
                                     "fidelity": fidelity(exact, m) if n <= fidelity_max_n else float("nan"),
                                     "discarded_weight": m.discarded, "seconds": time.perf_counter() - t0})
                    summary.append({"n": n, "depth": depth, "alpha": alpha, "seed": seed, "rule": rule,
                                    "peak_weight": w, "chi_break": first if first is not None else f">{chis[-1]}"})
        for alpha in alphas:
            for rule in rules:
                for chi in chis:
                    rows = [r for r in inst if r["n"] == n and r["alpha"] == alpha and r["rule"] == rule
                            and r["chi"] == chi]
                    if not rows:
                        continue
                    curves.append({"n": n, "alpha": alpha, "rule": rule, "chi": chi, "seeds": len(rows),
                                   "R": _mean(rows, "R"), "R_std": _std(rows, "R"),
                                   "fraction_R1": float(np.mean([r["R"] == 1.0 for r in rows])),
                                   "min_abs_Z": _mean(rows, "min_abs_Z"),
                                   "fidelity": _mean(rows, "fidelity"),
                                   "discarded_weight": _mean(rows, "discarded_weight")})
    table = []
    for n in ns:
        for alpha in alphas:
            for rule in rules:
                rows = [r for r in summary if r["n"] == n and r["alpha"] == alpha and r["rule"] == rule]
                if not rows:
                    continue
                vals = [r["chi_break"] for r in rows]
                nums = [v for v in vals if isinstance(v, int)]
                table.append({"n": n, "depth": block_a_depth(n), "alpha": alpha, "rule": rule, "seeds": len(rows),
                              "max_chi": chi_ladder(n, max_chi)[-1], "peak_weight": _mean(rows, "peak_weight"),
                              "chi_break_per_seed": " ".join(str(v) for v in vals),
                              "chi_break_median": float(np.median(nums)) if len(nums) == len(vals) else "n/a",
                              "reached": f"{len(nums)}/{len(vals)}"})
    setup = (f"block A depth = ceil(1.5 n); chi on powers of two up to min(2^(n/2), {max_chi}); rules {rules} "
             f"(random only for n <= {random_max_n}); {seeds} seeds per (n, alpha); R against the exact true peak; "
             f"chi_break = smallest chi on the ladder with R = 1; fidelity |<exact|mps>|^2 for n <= {fidelity_max_n}; "
             f"discarded_weight is the summed truncated Schmidt weight over the whole run.")
    return {"table": table, "instances": inst, "setup": setup, "curves": curves}


RUNS = {"depth_calibration": depth_calibration, "alpha_sweep": alpha_sweep, "chi_break": chi_break}


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(v) for k, v in r.items()})


def md_table(rows):
    keys = list(rows[0].keys())
    out = ["| " + " | ".join(keys) + " |", "|" + "|".join("---:" for _ in keys) + "|"]
    for r in rows:
        out.append("| " + " | ".join(_fmt(r[k]) for k in keys) + " |")
    return "\n".join(out)


def environment():
    import scipy
    return (f"Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, "
            f"{platform.system()} {platform.machine()}")


def save(name, out, seconds):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    write_csv(os.path.join(RESULTS_DIR, f"{name}.csv"), out["table"])
    write_csv(os.path.join(RESULTS_DIR, f"{name}_instances.csv"), out["instances"])
    parts = [f"# {name}", "", f"Produced by `python -m chibreak {name} --save` on {date.today().isoformat()} "
             f"in {seconds:.0f} s. {environment()}.", "", f"Setup: {out['setup']}", "", md_table(out["table"])]
    if "curves" in out:
        write_csv(os.path.join(RESULTS_DIR, f"{name}_curves.csv"), out["curves"])
        parts += ["", "Curves, averaged over seeds:", "", md_table(out["curves"])]
    with open(os.path.join(RESULTS_DIR, f"{name}.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(parts) + "\n")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    names = [a for a in argv if not a.startswith("--")]
    if not names or "--help" in argv or "-h" in argv:
        print("python -m chibreak <run> [<run> ...] [--save] [--quick]\n")
        for k in RUNS:
            print(f"  {k}")
        return 0
    for name in names:
        if name not in RUNS:
            print(f"unknown run {name}; choose from {', '.join(RUNS)}")
            return 1
    for name in names:
        t0 = time.perf_counter()
        out = RUNS[name](quick="--quick" in argv)
        seconds = time.perf_counter() - t0
        print(f"----- {name} ({seconds:.0f} s) -----")
        print(f"Setup: {out['setup']}\n")
        print(md_table(out["table"]))
        if "curves" in out:
            print("\nCurves:\n")
            print(md_table(out["curves"]))
        print(flush=True)
        if "--save" in argv:
            save(name, out, seconds)
    return 0
