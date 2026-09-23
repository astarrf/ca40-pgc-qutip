#!/usr/bin/env python3
"""Check Monte Carlo sampling convergence for the 397-nm Joshi PGC model.

Run several independent random seeds at increasing trajectory counts while
holding every physical parameter and the Fock cutoff fixed. This diagnoses
trajectory sampling only; repeat a separate cutoff comparison afterward.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import qutip as qt

from joshi_397nm import Parameters, build_model


def run_check(p: Parameters, counts: list[int], repeats: int, seed: int,
              tolerance: float, prefix: Path, progress: bool = True) -> dict:
    if repeats < 2:
        raise ValueError("At least two independent seeds per count are required")
    if seed < 0 or tolerance <= 0:
        raise ValueError("Seed must be nonnegative and tolerance positive")
    if len(counts) < 2 or counts != sorted(set(counts)):
        raise ValueError("Use at least two distinct, increasing trajectory counts")
    minimum = 2 * p.fock_cutoff if p.initial_nbar > 0 else 2
    if counts[0] < minimum:
        raise ValueError(f"At least {minimum} trajectories are needed for the "
                         "nonzero components of this initial mixture")

    runs_path = prefix.with_name(prefix.name + "_runs.csv")
    summary_path = prefix.with_name(prefix.name + "_summary.csv")
    diagnostics_path = prefix.with_name(prefix.name + "_diagnostics.json")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for path in (runs_path, summary_path, diagnostics_path):
        if path.exists():
            raise FileExistsError(f"{path} already exists; choose another --prefix")

    print("Building model once for all trajectory counts...", flush=True)
    H, rho0, c_ops, observables, details = build_model(p)
    times = np.linspace(0, p.duration_us, p.points)
    e_ops = list(observables.values())
    keys = list(observables)
    options = {"store_states": False, "keep_runs_results": False,
               "nsteps": 10000, "atol": 1e-8, "rtol": 1e-6,
               "progress_bar": "text" if progress else ""}
    all_results: dict[int, list[dict[str, np.ndarray]]] = {}

    # Write each completed run immediately, so an interrupted long scan still
    # leaves usable per-seed traces. A later complete scan needs a new prefix.
    with runs_path.open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["trajectories", "replicate", "seed", "time_us",
                         *keys])
        for count in counts:
            all_results[count] = []
            for replicate in range(repeats):
                run_seed = seed + counts.index(count) * repeats + replicate
                print(f"{count} trajectories, repeat {replicate + 1}/{repeats}, "
                      f"seed {run_seed}", flush=True)
                result = qt.mcsolve(
                    H, rho0, times, c_ops=c_ops, e_ops=e_ops,
                    ntraj=count, seeds=run_seed, options=options,
                )
                curves = {key: np.real(np.asarray(values))
                          for key, values in zip(keys, result.expect)}
                all_results[count].append(curves)
                for i, t in enumerate(times):
                    writer.writerow([count, replicate + 1, run_seed, t,
                                     *(curves[key][i] for key in keys)])
                stream.flush()

    summaries = {}
    with summary_path.open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["trajectories", "time_us", "mean_nbar",
                         "seed_sd_nbar", "seed_sem_nbar", "mean_p_top5"])
        for count in counts:
            nbar = np.stack([run["nbar"] for run in all_results[count]])
            top = np.stack([run["p_top5"] for run in all_results[count]])
            summaries[count] = {
                "mean": nbar.mean(axis=0),
                "sd": nbar.std(axis=0, ddof=1),
                "sem": nbar.std(axis=0, ddof=1) / math.sqrt(repeats),
                "top": top.mean(axis=0),
            }
            s = summaries[count]
            for i, t in enumerate(times):
                writer.writerow([count, t, s["mean"][i], s["sd"][i],
                                 s["sem"][i], s["top"][i]])

    comparisons = []
    for low, high in zip(counts[:-1], counts[1:]):
        drift = np.abs(summaries[high]["mean"] - summaries[low]["mean"])
        combined_sem = np.hypot(summaries[low]["sem"],
                                summaries[high]["sem"])
        comparisons.append({
            "low": low, "high": high,
            "max_trace_difference": float(drift.max()),
            "final_difference": float(drift[-1]),
            "max_two_combined_sem": float((2 * combined_sem).max()),
        })
    last = comparisons[-1]
    high_sem = float((2 * summaries[counts[-1]]["sem"]).max())
    # Conservative practical screen, not a statistical proof: the full trace
    # must be stable against doubling and repeated seeds must agree tightly.
    provisional_pass = (last["max_trace_difference"] <= tolerance
                        and last["max_two_combined_sem"] <= tolerance
                        and high_sem <= tolerance)
    diagnostics = {
        "parameters": asdict(p), "counts": counts, "repeats": repeats,
        "base_seed": seed, "tolerance_phonons": tolerance,
        "initial_nbar_after_truncation": details["initial_nbar_after_truncation"],
        "initial_tail_omitted": details["initial_tail_omitted"],
        "max_top5_at_largest_count": float(summaries[counts[-1]]["top"].max()),
        "comparisons": comparisons,
        "largest_count_max_two_seed_sem": high_sem,
        "provisional_sampling_pass": provisional_pass,
        "note": "Repeated-seed SEM is estimated from a small number of repeats; "
                "this is a practical screen, not proof of convergence. "
                "Fock-cutoff convergence must be checked separately.",
    }
    diagnostics_path.write_text(json.dumps(diagnostics, indent=2) + "\n")
    print("\nTrajectory-count comparison (phonons):")
    for item in comparisons:
        print(f"  {item['low']} -> {item['high']}: max |trace shift|="
              f"{item['max_trace_difference']:.4g}, final shift="
              f"{item['final_difference']:.4g}, max 2×combined seed SEM="
              f"{item['max_two_combined_sem']:.4g}")
    print(f"Largest-count max 2×seed SEM={high_sem:.4g}; requested "
          f"tolerance={tolerance:.4g}; provisional pass={provisional_pass}")
    print(f"Initial truncated nbar={details['initial_nbar_after_truncation']:.4f}; "
          f"omitted tail={details['initial_tail_omitted']:.3g}; "
          f"largest-count max top-five population="
          f"{diagnostics['max_top5_at_largest_count']:.3g}")
    print(f"Saved {runs_path}, {summary_path}, {diagnostics_path}")
    if details["initial_tail_omitted"] > 1e-3 or diagnostics["max_top5_at_largest_count"] > 1e-3:
        print("Fock cutoff may be inadequate even if sampling passes.")
    return diagnostics


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts", help="increasing counts, e.g. 128,256,512; "
                    "default is 2N,4N,8N")
    ap.add_argument("--repeats", type=int, default=3,
                    help="independent seeds per count (default: 3)")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--tolerance", type=float, default=0.25,
                    help="target absolute nbar difference across the full trace")
    ap.add_argument("--prefix", type=Path, default=Path("mc_convergence"))
    ap.add_argument("--no-progress", action="store_true")
    ap.add_argument("--benchmark", action="store_true")
    for name in Parameters.__dataclass_fields__:
        default = getattr(Parameters(), name)
        ap.add_argument("--" + name.replace("_", "-"),
                        type=int if isinstance(default, int) else float)
    if argv is None and Path(sys.argv[0]).stem == "ipykernel_launcher":
        argv = []
    args = ap.parse_args(argv)
    p = Parameters()
    if args.benchmark:
        p = replace(p, trap_khz=1088.0, detuning_mhz=210.0,
                    saturation=3 * 1.088 * 1.35 / 210, delta_khz=60.0)
    updates = {name: getattr(args, name) for name in Parameters.__dataclass_fields__
               if getattr(args, name) is not None}
    p = replace(p, **updates)
    minimum = 2 * p.fock_cutoff if p.initial_nbar > 0 else 2
    counts = ([minimum, 2 * minimum, 4 * minimum] if args.counts is None
              else [int(item.strip()) for item in args.counts.split(",")])
    run_check(p, counts, args.repeats, args.seed, args.tolerance,
              args.prefix, progress=not args.no_progress)


if __name__ == "__main__":
    main()
