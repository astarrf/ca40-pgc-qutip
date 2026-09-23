#!/usr/bin/env python3
"""Single-ion 397-nm lin-perp-lin PGC, four electronic levels + one axial mode.

Based on Joshi et al., New J. Phys. 22, 103013 (2020), Appendix A:
https://arxiv.org/html/2006.15104#A1

All Hamiltonian frequencies and decay rates are angular radians/microsecond;
time is in microseconds and hbar=1. This is a four-level approximation: decay
into D_3/2, the 866-nm repumper, Zeeman shifts, micromotion, and other modes
are omitted. In particular it is not a quantitative prediction for a specific
apparatus until these and the beam geometry/intensity are calibrated.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import qutip as qt

HBAR = 1.054571817e-34
AMU = 1.66053906660e-27
MASS_CA40 = 39.96259098 * AMU


@dataclass(frozen=True)
class Parameters:
    trap_khz: float = 690.0           # omega_z / 2pi
    initial_nbar: float = 15.0
    detuning_mhz: float = 30.0        # Delta / 2pi, positive is blue
    delta_khz: float = 60.0           # (omega_1 - omega_2) / 2pi
    saturation: float = 0.063         # Joshi's s, per-beam stretched reference
    linewidth_mhz: float = 21.6      # Gamma / 2pi (total P_1/2 linewidth)
    wavelength_nm: float = 397.0
    axial_projection: float = 1.0     # |k_z|/k; both beams assumed symmetric
    phase_rad: float = 0.0            # static gradient phase at trap center
    fock_cutoff: int = 128            # keeps states |0> through |N-1>
    duration_us: float = 100.0
    points: int = 101


def angular_mhz(value: float) -> float:
    """MHz (ordinary cycles/us) -> rad/us."""
    return 2 * np.pi * value


def eta_for(p: Parameters) -> float:
    omega_si = 2 * np.pi * p.trap_khz * 1e3
    z0 = math.sqrt(HBAR / (2 * MASS_CA40 * omega_si))
    return 2 * np.pi / (p.wavelength_nm * 1e-9) * z0 * p.axial_projection


def rabi_for(p: Parameters) -> float:
    """Stretched-reference Ω (rad/us): s=(Ω²/2)/(Γ²/4+Δ²)."""
    gamma = angular_mhz(p.linewidth_mhz)
    detuning = angular_mhz(p.detuning_mhz)
    return math.sqrt(2 * p.saturation * (gamma**2 / 4 + detuning**2))


def validate(p: Parameters) -> None:
    if p.trap_khz <= 0 or p.linewidth_mhz <= 0 or p.wavelength_nm <= 0:
        raise ValueError("Trap frequency, linewidth, and wavelength must be positive")
    if p.initial_nbar < 0 or p.saturation < 0 or p.fock_cutoff < 2:
        raise ValueError("nbar and saturation must be nonnegative; cutoff >= 2")
    if not 0 <= p.axial_projection <= 1:
        raise ValueError("axial_projection must be between 0 and 1")
    if p.duration_us <= 0 or p.points < 2:
        raise ValueError("duration_us > 0 and points >= 2 are required")


def build_model(p: Parameters):
    """Return H(t), rho(0), collapse operators, observables, and diagnostics.

    Tensor order is electronic x motion. Electronic basis is
    [g+, g-, e+, e-]. The two driven sigma transitions are g- <-> e+
    (sigma+) and g+ <-> e- (sigma-). We use Joshi Appendix A Eq. (15):
    theta_+ = phi-pi/4 and theta_- = phi+pi/4. A frequency difference
    delta moves theta by delta*t/2. Reversing delta reverses its travel.
    """
    validate(p)
    N = p.fock_cutoff
    e = [qt.basis(4, j) for j in range(4)]
    eg, ig = qt.qeye(4), qt.qeye(N)
    a = qt.destroy(N)
    n = a.dag() * a
    eta = eta_for(p)
    # Exponentiate the full projected position operator, with no LD series.
    Dplus = (1j * eta * (a + a.dag())).expm()
    Dminus = Dplus.dag()
    tensor = qt.tensor
    lowering_plus = e[1] * e[2].dag()    # e+ -> g-
    lowering_minus = e[0] * e[3].dag()   # e- -> g+
    prefactor = rabi_for(p) / (2 * math.sqrt(3))
    theta_plus = p.phase_rad - np.pi / 4
    theta_minus = p.phase_rad + np.pi / 4
    # Coefficients of e^{+i delta t/2} and e^{-i delta t/2}.
    # Includes the Hermitian conjugates of both optical traveling waves.
    B = prefactor * (
        tensor(np.exp(1j * theta_plus) * lowering_plus, Dplus)
        + tensor(np.exp(1j * theta_minus) * lowering_minus, Dplus)
        + tensor(np.exp(1j * theta_plus) * lowering_plus.dag(), Dplus)
        + tensor(np.exp(1j * theta_minus) * lowering_minus.dag(), Dplus)
    )
    omega = 2 * np.pi * p.trap_khz / 1000  # rad/us
    detuning = angular_mhz(p.detuning_mhz)
    projector_e = e[2] * e[2].dag() + e[3] * e[3].dag()
    H0 = tensor(eg, omega * n) - detuning * tensor(projector_e, ig)
    delta = 2 * np.pi * p.delta_khz / 1000
    if delta == 0:
        H = H0 + B + B.dag()
    else:
        H = [H0, [B, lambda t, **kwargs: np.exp(0.5j * delta * t)],
             [B.dag(), lambda t, **kwargs: np.exp(-0.5j * delta * t)]]

    # Joshi Appendix A Table 1: recoil quadrature for sigma (m=+/-1)
    # and pi (m=0). Squared weights each sum to one. Each electronic
    # transition has three separate, incoherent Lindblad jump operators.
    # The second recoil moments are 2/5 for sigma emission and 1/5 for
    # pi emission; CG branching gives (2/3)(2/5)+(1/3)(1/5)=1/3 overall.
    # This is the paper's axial recoil factor alpha. Check geometry before use.
    recoil_sigma = ((-1, .2), (0, .6), (+1, .2))
    recoil_pi = ((-1, .1), (0, .8), (+1, .1))
    transitions = (
        (e[1] * e[2].dag(), 2 / 3, recoil_sigma),
        (e[0] * e[3].dag(), 2 / 3, recoil_sigma),
        (e[0] * e[2].dag(), 1 / 3, recoil_pi),
        (e[1] * e[3].dag(), 1 / 3, recoil_pi),
    )
    gamma = angular_mhz(p.linewidth_mhz)
    recoil = {-1: Dminus, 0: ig, +1: Dplus}
    c_ops = [math.sqrt(gamma * branching * weight) * tensor(op, recoil[q])
             for op, branching, choices in transitions for q, weight in choices]

    # Truncate and renormalize the input thermal distribution. This differs
    # from the infinite thermal state when N is too small, so report its tail.
    ratio = p.initial_nbar / (p.initial_nbar + 1) if p.initial_nbar else 0.0
    weights = (1 - ratio) * ratio ** np.arange(N)
    omitted_initial_probability = ratio**N
    rho_m = qt.Qobj(np.diag(weights / weights.sum()), dims=[[N], [N]])
    rho_g = .5 * (e[0] * e[0].dag() + e[1] * e[1].dag())
    rho0 = tensor(rho_g, rho_m)
    top = sum((qt.basis(N, j) * qt.basis(N, j).dag()
               for j in range(max(0, N - 5), N)), qt.Qobj(np.zeros((N, N))))
    observables = {
        "nbar": tensor(eg, n),
        "p_excited": tensor(projector_e, ig),
        "p_top5": tensor(eg, top),
        "p_n0": tensor(eg, qt.basis(N, 0) * qt.basis(N, 0).dag()),
    }
    details = {"eta": eta, "omega_stretched_mhz": rabi_for(p) / (2 * np.pi),
               "xi": p.detuning_mhz * p.saturation / (3 * p.trap_khz / 1000),
               "initial_tail_omitted": omitted_initial_probability,
               "initial_nbar_after_truncation": float(qt.expect(observables["nbar"], rho0))}
    return H, rho0, c_ops, observables, details


def simulate(p: Parameters):
    H, rho0, c_ops, obs, details = build_model(p)
    times = np.linspace(0, p.duration_us, p.points)
    # Store expectations only, never the (potentially large) density matrices.
    options = {"store_states": False, "nsteps": 10000, "atol": 1e-8, "rtol": 1e-6}
    result = qt.mesolve(H, rho0, times, c_ops=c_ops, e_ops=list(obs.values()),
                       args={}, options=options)
    columns = {name: np.real(np.asarray(values))
               for name, values in zip(obs, result.expect)}
    return times, columns, details


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark", action="store_true", help="Joshi-like +210 MHz, 1.088 MHz, xi=1.35, 60 kHz moving gradient")
    ap.add_argument("--trap-khz", type=float)
    ap.add_argument("--initial-nbar", type=float)
    ap.add_argument("--detuning-mhz", type=float)
    ap.add_argument("--delta-khz", type=float)
    ap.add_argument("--saturation", type=float)
    ap.add_argument("--linewidth-mhz", type=float)
    ap.add_argument("--wavelength-nm", type=float)
    ap.add_argument("--axial-projection", type=float)
    ap.add_argument("--phase-rad", type=float)
    ap.add_argument("--fock-cutoff", type=int)
    ap.add_argument("--duration-us", type=float)
    ap.add_argument("--points", type=int)
    ap.add_argument("--output", type=Path, default=Path("pgc_trace.csv"))
    args = ap.parse_args()
    p = Parameters()
    if args.benchmark:
        p = replace(p, trap_khz=1088.0, detuning_mhz=210.0,
                    saturation=3 * 1.088 * 1.35 / 210, delta_khz=60.0)
    updates = {name: getattr(args, name) for name in Parameters.__dataclass_fields__
               if getattr(args, name) is not None}
    p = replace(p, **updates)
    times, columns, details = simulate(p)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_us", *columns])
        writer.writerows(zip(times, *(columns.values())))
    print(f"Saved {args.output.resolve()}")
    print("eta={eta:.4f}, Ω_stretched/2π={omega_stretched_mhz:.3f} MHz, "
          "xi={xi:.3f}".format(**details))
    print("Initial truncated nbar={initial_nbar_after_truncation:.4f}, "
          "omitted thermal probability={initial_tail_omitted:.3g}".format(**details))
    print(f"Final nbar={columns['nbar'][-1]:.4f}, excited={columns['p_excited'][-1]:.4g}, "
          f"top-five Fock population={columns['p_top5'][-1]:.3g}")
    if details["initial_tail_omitted"] > 1e-3 or columns["p_top5"].max() > 1e-3:
        print("CUTOFF WARNING: increase --fock-cutoff and compare traces.")
    if args.benchmark:
        print("Joshi reports ~66,000/s cooling near this setting; this is a "
              "qualitative cross-check, not an asserted reproduction. Fit the "
              "transient and test Fock/phase convergence before comparison.")


if __name__ == "__main__":
    main()
