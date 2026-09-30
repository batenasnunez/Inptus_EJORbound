"""
frontier_experiments.py
========================

Computes and plots attainable graduation-performance frontiers for the
stage-structured progression model of

    "Attainable Graduation Frontiers at Fixed Annual Performance:
    Stage-Structured Optimization and Evidence from Spanish Higher
    Education"

Edit the CONFIG block below to run a concrete experiment: choose the
programme duration d and a persistence profile (uniform, linearly
varying across stages, or a custom array). Running the script produces
a figure with:

  - The true outer envelope of the attainable region: the MAXIMUM
    frontier under full persistence (lambda=1) and the MINIMUM
    frontier under zero persistence (lambda=0). Solid lines.
  - The two "inner" reference curves: the MINIMUM under full
    persistence and the MAXIMUM under zero persistence. Dotted lines.
  - The homogeneous-configuration benchmark curves (p_k = a for every
    stage) at lambda=0 and lambda=1. Solid lines. These are NOT
    frontiers -- they mark the special homogeneous configuration, not
    an extremum -- but are informative reference curves.
  - The persistence-specific frontier (min and max) for the chosen
    profile (uniform / linear / custom), computed numerically.
  - A random sample of feasible (p_A, OTGR(+1)) points generated under
    the chosen profile, illustrating the range of attainable outcomes
    (as in Figure 2 of the paper). This is the ONLY scattered element
    in the figure; every frontier/reference curve is a plain line.

Caching
-------
The full-persistence reference curves (max and min) require numerical
optimization and have no closed form for general d, so they are the
expensive part of this script. They are cached to
"full_persistence_cache_d<D>.npz" next to this script. On a repeat run
with the same d, grid, and solver settings, the cache is loaded instead
of recomputed, so subsequent runs spend their time only on the
persistence-specific (intermediate) profile, which is what actually
changes between experiments. Delete the cache file, or change
N_RESTARTS_REF/SEED/D/A_GRID, to force a recomputation.

The zero-persistence frontiers and the homogeneous benchmark curves are
closed-form (no optimization, no caching needed).

Numerical protocol (for reproducibility)
-----------------------------------------
Global minimization/maximization of OTGR(+1) subject to p_A(p) = a is
performed with scipy.optimize.minimize, method='SLSQP', an equality
constraint on p_A, bounds p_k in [1e-4, 1.0], ftol=1e-14, maxiter=1000.
For each (a, lambda-profile) pair, N_RESTARTS independent random
starting points (numpy Generator with a fixed seed) are tried; only
solutions satisfying the bounds and the equality constraint to within
1e-6 are accepted, and the best feasible objective value found across
restarts is reported. This is a global-search heuristic (multi-start
local SLSQP), not a certificate of global optimality; RUN_VALIDATION
checks it against the exact zero-persistence solution, where the true
global optimum is known analytically.

Parallelization
---------------
Each point of A_GRID is an independent optimization problem, so the
frontier computations (validation, full-persistence reference, and the
profile-specific frontier) are parallelized across grid points with
concurrent.futures.ProcessPoolExecutor, using N_WORKERS processes (set
in CONFIG; default is all logical cores minus one). Set N_WORKERS = 1
to fall back to plain sequential execution.

Requires: numpy, scipy, matplotlib.
"""

import os
import numpy as np
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# CONFIG -- edit this block to define an experiment
# ----------------------------------------------------------------------

D = 4                       # programme duration (number of stages)

# Persistence profile: choose one of 'uniform', 'linear', 'custom'
PERSISTENCE_MODE = 'uniform'

LAMBDA_VALUE = 0.5                       # used if mode == 'uniform'
LAMBDA_START, LAMBDA_END = 0.2, 0.8      # used if mode == 'linear'
LAMBDA_CUSTOM = np.array([0.9, 0.2, 0.6, 0.4])   # used if mode == 'custom'

A_GRID = np.linspace(0.02, 0.99, 30)     # annual-performance values to scan

N_RESTARTS_PROFILE = 150    # restarts for the profile-specific frontier
                             # (recomputed every run)
N_RESTARTS_REF = 300        # restarts for the full-persistence reference
                             # curves (cached after the first run)
SEED = 0                    # RNG seed, for reproducibility
BOUNDS_LOW = 1e-4
FTOL = 1e-14

N_WORKERS = max(1, (os.cpu_count() or 2) - 1)   # parallel processes;
                                                 # set to 1 for sequential

RUN_VALIDATION = True        # sanity-check the optimizer against the
                              # exact zero-persistence solution
N_SAMPLES = 3000              # random points illustrating attainable cases
SAMPLE_SEED = 1

CACHE_FILE = f"full_persistence_cache_d{D}.npz"
OUTPUT_FIGURE = "frontier_experiment.png"

# ----------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------

def build_lambda(d, mode, lambda_value=None, lambda_start=None,
                  lambda_end=None, lambda_custom=None):
    """Build a persistence-probability vector of length d."""
    if mode == 'uniform':
        return np.full(d, lambda_value, dtype=float)
    elif mode == 'linear':
        return np.linspace(lambda_start, lambda_end, d)
    elif mode == 'custom':
        lam = np.asarray(lambda_custom, dtype=float)
        if lam.shape[0] != d:
            raise ValueError(f"lambda_custom must have length d={d}, got {lam.shape[0]}")
        return lam
    else:
        raise ValueError(f"Unknown PERSISTENCE_MODE: {mode!r}")


def survival_factors(p, lam):
    """Pi_k = prod_{i<=k} p_i / (1 - lambda_i * (1 - p_i))."""
    q = 1.0 - p
    denom = 1.0 - lam * q
    return np.cumprod(p / denom)


def p_A(p, lam):
    """Stationary-population-weighted annual performance rate."""
    Pi = survival_factors(p, lam)
    return np.sum(Pi) / np.sum(Pi / p)


def OTGR1(p, lam):
    """OTGR(+1): graduation within nominal duration plus one year."""
    q = 1.0 - p
    return np.prod(p) * (1.0 + np.sum(lam * q))


def zero_persistence_exact(a, d):
    """Closed-form Theorem 2 frontiers at lambda = 0."""
    g_max = a / (d - (d - 1) * a)
    g_min = max(d * a - (d - 1), 0.0)
    return g_min, g_max


def homogeneous_curve(a, d, lam):
    """OTGR(+1) at the homogeneous configuration p_k = a for all k,
    with constant persistence lam. Note p_A = a exactly at this
    configuration, for any lam (Eq. homogeneous_performance)."""
    return a**d * (1.0 + lam * d * (1.0 - a))


# ----------------------------------------------------------------------
# Numerical optimization (documented protocol, see module docstring)
# ----------------------------------------------------------------------

def optimize_G(a, lam, sense, d, n_restarts, seed):
    """Solve one (a, lambda) point directly (sequential). Kept for
    single-point experimentation; frontier_over_grid uses the
    picklable, process-parallel _solve_for_a below for the grid case."""
    sign = 1.0 if sense == 'min' else -1.0
    bounds = [(BOUNDS_LOW, 1.0)] * d
    cons = [{'type': 'eq', 'fun': lambda p: p_A(p, lam) - a}]
    rng = np.random.default_rng(seed)
    best_val, best_p = None, None
    for _ in range(n_restarts):
        p0 = rng.uniform(0.05, 1.0, size=d)
        try:
            res = minimize(lambda p: sign * OTGR1(p, lam), p0, method='SLSQP',
                            bounds=bounds, constraints=cons,
                            options={'ftol': FTOL, 'maxiter': 1000})
        except Exception:
            continue
        if not res.success:
            continue
        p = res.x
        if np.any(p < BOUNDS_LOW - 1e-6) or np.any(p > 1.0 + 1e-6):
            continue
        if abs(p_A(p, lam) - a) > 1e-6:
            continue
        val = OTGR1(p, lam)
        if best_val is None or (sense == 'min' and val < best_val - 1e-10) or \
           (sense == 'max' and val > best_val + 1e-10):
            best_val, best_p = val, p.copy()
    return best_val, best_p


def _solve_for_a(task):
    """Worker for one grid point, run in a separate process. `task` is a
    plain tuple so it can be pickled and sent to a worker process; the
    logic is identical to optimize_G (kept in sync manually)."""
    a, lam, sense, d, n_restarts, seed, bounds_low, ftol = task
    sign = 1.0 if sense == 'min' else -1.0
    bounds = [(bounds_low, 1.0)] * d
    cons = [{'type': 'eq', 'fun': lambda p: p_A(p, lam) - a}]
    rng = np.random.default_rng(seed)
    best_val = None
    for _ in range(n_restarts):
        p0 = rng.uniform(0.05, 1.0, size=d)
        try:
            res = minimize(lambda p: sign * OTGR1(p, lam), p0, method='SLSQP',
                            bounds=bounds, constraints=cons,
                            options={'ftol': ftol, 'maxiter': 1000})
        except Exception:
            continue
        if not res.success:
            continue
        p = res.x
        if np.any(p < bounds_low - 1e-6) or np.any(p > 1.0 + 1e-6):
            continue
        if abs(p_A(p, lam) - a) > 1e-6:
            continue
        val = OTGR1(p, lam)
        if best_val is None or (sense == 'min' and val < best_val - 1e-10) or \
           (sense == 'max' and val > best_val + 1e-10):
            best_val = val
    return best_val


def frontier_over_grid(a_grid, lam, d, sense, n_restarts, seed,
                        n_workers=None):
    if n_workers is None:
        n_workers = N_WORKERS
    tasks = [(a, lam, sense, d, n_restarts, seed, BOUNDS_LOW, FTOL)
             for a in a_grid]
    if n_workers <= 1:
        results = [_solve_for_a(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            results = list(ex.map(_solve_for_a, tasks))
    return np.array([r if r is not None else np.nan for r in results])


def validate_zero_persistence(a_grid, d, n_restarts, seed):
    lam0 = np.zeros(d)
    num_min = frontier_over_grid(a_grid, lam0, d, 'min', n_restarts, seed)
    num_max = frontier_over_grid(a_grid, lam0, d, 'max', n_restarts, seed)
    errs_min, errs_max = [], []
    for a, gmin_n, gmax_n in zip(a_grid, num_min, num_max):
        gmin_e, gmax_e = zero_persistence_exact(a, d)
        errs_min.append(abs(gmin_n - gmin_e))
        errs_max.append(abs(gmax_n - gmax_e))
    return max(errs_min), max(errs_max)


# ----------------------------------------------------------------------
# Cached full-persistence reference curves
# ----------------------------------------------------------------------

def load_or_compute_full_persistence(a_grid, d, n_restarts, seed, cache_file):
    if os.path.exists(cache_file):
        data = np.load(cache_file)
        same_grid = (data['a_grid'].shape == a_grid.shape and
                     np.allclose(data['a_grid'], a_grid))
        same_settings = (int(data['n_restarts']) == n_restarts and
                          int(data['seed']) == seed and int(data['d']) == d)
        if same_grid and same_settings:
            print(f"  Loaded cached full-persistence frontiers from {cache_file}")
            return data['fp_min'], data['fp_max']
        else:
            print("  Cache found but settings differ (grid/d/restarts/seed) "
                  "-- recomputing.")

    print("  No usable cache found -- computing full-persistence frontiers "
          f"numerically ({n_restarts} restarts per point) ...")
    lam1 = np.ones(d)
    fp_min = frontier_over_grid(a_grid, lam1, d, 'min', n_restarts, seed)
    fp_max = frontier_over_grid(a_grid, lam1, d, 'max', n_restarts, seed)
    np.savez(cache_file, a_grid=a_grid, fp_min=fp_min, fp_max=fp_max,
             n_restarts=n_restarts, seed=seed, d=d)
    print(f"  Cached full-persistence frontiers to {cache_file}")
    return fp_min, fp_max


# ----------------------------------------------------------------------
# Random sample of attainable (a, G) points, under the chosen profile
# ----------------------------------------------------------------------

def random_sample(d, lam, n_samples, seed):
    rng = np.random.default_rng(seed)
    p_matrix = rng.uniform(0.01, 1.0, size=(n_samples, d))
    a_vals = np.array([p_A(p, lam) for p in p_matrix])
    g_vals = np.array([OTGR1(p, lam) for p in p_matrix])
    return a_vals, g_vals


# ----------------------------------------------------------------------
# Main experiment
# ----------------------------------------------------------------------

def main():
    lam_profile = build_lambda(D, PERSISTENCE_MODE, LAMBDA_VALUE,
                                LAMBDA_START, LAMBDA_END, LAMBDA_CUSTOM)
    print(f"Programme duration d = {D}")
    print(f"Persistence profile ({PERSISTENCE_MODE}): {np.round(lam_profile, 4)}")
    print(f"Parallel workers: {N_WORKERS} "
          f"({'sequential' if N_WORKERS <= 1 else 'multi-core'})")

    if RUN_VALIDATION:
        print("\nValidating optimizer against exact zero-persistence solution "
              f"(d={D}, {N_RESTARTS_PROFILE} restarts, seed={SEED})...")
        max_err_min, max_err_max = validate_zero_persistence(
            A_GRID, D, N_RESTARTS_PROFILE, SEED)
        print(f"  max |numeric - exact|:  G_min: {max_err_min:.3e}   "
              f"G_max: {max_err_max:.3e}")

    print("\nZero-persistence envelope (closed form) ...")
    zp_min = np.array([zero_persistence_exact(a, D)[0] for a in A_GRID])
    zp_max = np.array([zero_persistence_exact(a, D)[1] for a in A_GRID])

    print("Homogeneous benchmark curves (closed form) ...")
    hom0 = np.array([homogeneous_curve(a, D, 0.0) for a in A_GRID])
    hom1 = np.array([homogeneous_curve(a, D, 1.0) for a in A_GRID])

    print("\nFull-persistence reference frontiers (numerical, cached) ...")
    fp_min, fp_max = load_or_compute_full_persistence(
        A_GRID, D, N_RESTARTS_REF, SEED, CACHE_FILE)

    print(f"\nPersistence-specific frontier for the '{PERSISTENCE_MODE}' "
          f"profile (numerical, {N_RESTARTS_PROFILE} restarts) ...")
    prof_min = frontier_over_grid(A_GRID, lam_profile, D, 'min',
                                   N_RESTARTS_PROFILE, SEED)
    prof_max = frontier_over_grid(A_GRID, lam_profile, D, 'max',
                                   N_RESTARTS_PROFILE, SEED)

    print(f"\nGenerating random sample of {N_SAMPLES} attainable points "
          f"under the chosen profile ...")
    sample_a, sample_g = random_sample(D, lam_profile, N_SAMPLES, SAMPLE_SEED)

    # Consistency check: no sampled point should fall outside the outer
    # envelope (interpolated onto each sample's own a value), and few
    # (ideally none) should fall outside the profile-specific frontier --
    # a breach there would signal that N_RESTARTS_PROFILE/REF is too low.
    tol = 1e-3
    env_hi = np.interp(sample_a, A_GRID, fp_max)
    env_lo = np.interp(sample_a, A_GRID, zp_min)
    prof_hi = np.interp(sample_a, A_GRID, prof_max)
    prof_lo = np.interp(sample_a, A_GRID, prof_min)
    n_env_breach = np.sum((sample_g > env_hi + tol) | (sample_g < env_lo - tol))
    n_prof_breach = np.sum((sample_g > prof_hi + tol) | (sample_g < prof_lo - tol))
    print(f"  Consistency check: {n_env_breach}/{N_SAMPLES} sample points "
          f"outside the outer envelope; {n_prof_breach}/{N_SAMPLES} outside "
          f"the '{PERSISTENCE_MODE}' profile frontier "
          f"(tolerance {tol:g}).")

    # ---------------- Plot ----------------
    fig, ax = plt.subplots(figsize=(7.5, 6.5))

    # Random sample cloud (only scattered element in the figure)
    ax.scatter(sample_a, sample_g, s=4, color='0.6', alpha=0.35,
               linewidths=0, zorder=1, label='random sample (chosen profile)')

    # Outer envelope: full-persistence MAX, zero-persistence MIN (solid)
    ax.plot(A_GRID, fp_max, color='tab:red', lw=2.4, solid_capstyle='round',
            zorder=4, label=r'$\lambda=1$ maximum (envelope, numerical)')
    ax.plot(A_GRID, zp_min, color='tab:blue', lw=2.4, solid_capstyle='round',
            zorder=4, label=r'$\lambda=0$ minimum (envelope, exact)')

    # Inner reference curves: full-persistence MIN, zero-persistence MAX (dotted)
    ax.plot(A_GRID, fp_min, color='tab:red', lw=1.6, linestyle=':',
            zorder=3, label=r'$\lambda=1$ minimum (numerical)')
    ax.plot(A_GRID, zp_max, color='tab:blue', lw=1.6, linestyle=':',
            zorder=3, label=r'$\lambda=0$ maximum (exact)')

    # Homogeneous benchmark curves (solid, thinner, distinct shades)
    ax.plot(A_GRID, hom1, color='darkorange', lw=1.4,
            zorder=2, label=r'homogeneous, $\lambda=1$ (exact)')
    ax.plot(A_GRID, hom0, color='mediumpurple', lw=1.4,
            zorder=2, label=r'homogeneous, $\lambda=0$ (exact)')

    # Persistence-specific frontier for the chosen profile
    ax.plot(A_GRID, prof_max, color='tab:green', lw=1.8,
            zorder=3, label=f'{PERSISTENCE_MODE} profile maximum (numerical)')
    ax.plot(A_GRID, prof_min, color='tab:green', lw=1.8, linestyle='--',
            zorder=3, label=f'{PERSISTENCE_MODE} profile minimum (numerical)')

    ax.set_xlabel(r'Annual performance rate $a$')
    ax.set_ylabel(r'$OTGR(+1)$')
    ax.set_title(f'Attainable frontiers, d={D}, profile={PERSISTENCE_MODE} '
                 f'{np.round(lam_profile, 2)}')
    ax.legend(fontsize=7.5, loc='upper left')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    fig.tight_layout()
    fig.savefig(OUTPUT_FIGURE, dpi=200)
    print(f"\nFigure written to {OUTPUT_FIGURE}")


if __name__ == "__main__":
    main()
