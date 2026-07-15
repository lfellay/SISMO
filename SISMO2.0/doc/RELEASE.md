# SISMO 2.0 — release record (updated 2026-07-14)

**SISMO — Split Inhomogeneous Solver for Modelling Oscillations**

Adiabatic stellar oscillations with a SPLIT (inhomogeneous) Poisson treatment
on a scalar second-order Cowling core.  No coupled solver anywhere.

## Canonical run (the validated release configuration)

The normal invocation has no numerical command-line flags:

    export OMP_NUM_THREADS=<cores>
    sismo MODEL.osc.mod OUT

SISMO loads these validated settings from the plain-text `sismo.conf` file:

```text
# Mode range
l_min = 1
l_max = 2
g_min = 1
g_max = 108

# Physics and output
use_poisson = true
takata_closure = true
write_eigenfunctions = false

# Numerical controls
scan_points = 2000
eps_g = 2.5
max_iterations = 150
tolerance = 1.0e-8
inner_poisson_iterations = 1
poisson_relaxation = 0.4
poisson_inner_tolerance = 0.0
```

A custom file can be passed as the optional third positional argument:

    sismo MODEL.osc.mod OUT /path/to/custom.conf

The lookup precedence is the explicit file, `SISMO_CONFIG`, `./sismo.conf`,
then `sismo.conf` beside the executable. Each non-comment line has the form
`key = value`; `#` and `!` introduce comments. Unknown and duplicate keys are
errors.

MODEL is intSISMO output (`.osc.mod` with its `.grid.d` sidecar), or a legacy
`sta1.d` structure. `g_max` is the wanted order plus approximately eight modes
of headroom because the top of the Sturm window includes f/p modes. Set
`use_poisson = false` for a pure Cowling spectrum. `scan_points` is the
endpoint-inclusive number of initial samples per degree. Set
`write_eigenfunctions = true` only when full-grid mechanical `.eig` files are
needed; their potential-related columns are currently zero.

## Method (what this version uses)

1. OPERATOR: the syst05 mechanical pair (rows 1,2; row-1 derivative weight
   omega^2) reduced to one scalar second-order equation in y1; exponentially
   fitted staggered box scheme (y1 at nodes, y2 eliminated at midpoint faces
   with cell-local integrating factors) -> symmetric-equivalent tridiagonal
   with positive off-diagonal products.  Full model grid, no multigrid/remesh.
   Boundary rows = the block solver's own (surface delta-p = 0; centre per case).
2. STURM SCAN (per l): parallel sweep uniform-in-period + parallel bisection of
   the parity of (negative LU pivots + negative face denominators) -- the
   face-denominator term cancels the Lamb-crossing pole artifacts, so the
   parity flips exactly at eigenvalues.  Complete labeled spectrum (g + f + p);
   the asymptotic relation sets ONLY the frequency window, never a seed.
3. SPLIT REFINEMENT (per mode, OpenMP-parallel): seed = Cowling eigenvalue;
   iterate { phi <- shooting Poisson(mechanics) [frozen source, relax 0.4];
   forced response; sign-aligned incremental shape update; clamped Newton
   frequency step (finite-difference T'); best-iterate return }.
4. DIPOLE CLOSURE (`takata_closure = true`): for l=1 the shooting amplitude A
   (phi = phi_p + A*phi_h) is set by least-squares Takata J = 0 over the star
   (momentum conservation, PASJ 58, 893) instead of the vacuum surface BC.
   J is linear in A through BOTH the phi terms and the pp term.  Parameter-free;
   replaced every earlier l=1 heuristic (blends, gates, special seeds).

## Validation (1M5 intSISMO model, 25k points, l=1,2 n=1..100)

  vs OSC (same model = pure solver difference):
      l=1: median 0.0019%, 99% <0.1%, max 0.21% over all 98 modes
      l=2: median 0.0018%, 100% <0.1%, max 0.10% over all 100 modes
      (f-mode 0.07%; every mode of the problem <0.21%)
  vs GYRE (its own MESA grid; the ~0.22% floor is the mesa2SISMO re-mesh):
      l=1: 0.216% (100% <1%)   l=2: 0.217% (100% <1%)
      (OSC itself vs GYRE: 0.218% / 0.216%)
  Reference output: doc/release_reference_1M5.sismo (bit-identity checked
  after every cleanup step).  Comparison figures + scripts:
  ../1M5/compare_core2_{gyre,osc}.png, make_compare.py, make_compare_osc.py.

## Timing (18 threads, same host, l=1,2 n=1..100)

  GYRE 8.1                 9.4 s wall     71 s CPU   (1131-pt MESA profile)
  OSC  (scang, oscQI=2)   72.5 s        1095 s       (25k intSISMO)
  SISMO 1 (canonical)    72.7 s         905 s       (25k intSISMO)
  SISMO 2.0 (optimized)   1.9 s          28 s       (25k intSISMO, full grid)
  -> about 38x less CPU time than OSC and SISMO 1 on the same model; 216 modes
     including f/p.

## Optimization validation (2026-07-14)

Each step below was clean-built and timed twice on the same 24,992-point model
with 18 OpenMP threads and the canonical 216-mode command. CPU time is the more
stable comparison because wall time is short and sensitive to host scheduling.

| Stage | Wall time | CPU time | Numerical validation |
|---|---:|---:|---|
| Previous implementation | 2.86 s | 47.15 s | reference |
| Degree/operator cache | 2.20 s | 32.54 s | byte-identical |
| Reuse operator and face data | 2.10 s | 31.23 s | max relative change 1.66e-10 |
| Reusable solve/scan workspaces | 2.09 s | 31.01 s | byte-identical to preceding stage |
| Lazy eigenfunction storage | 2.10 s | 30.94 s | byte-identical; peak RSS 806 -> 372 MB |
| Default 2,000-point scan | 1.88 s | 27.99 s | same 216 ordered modes and 149/152 scan roots |

Relative to the previous implementation, the final result uses 40.6% less CPU
time and 34.3% less wall time. Default and explicit `scan_points = 2000`
configurations are byte-identical, as are 1-thread and 18-thread outputs.
Against OSC, the final medians are 0.00194% (l=1) and 0.00176% (l=2); the
maxima are 0.204% and 0.0558%, and every compared mode remains below 1%.

Compared with the previous 9,201-sample sweep, seven Cowling roots change by
at most 7.32e-10 relatively. The difficult low-order modes that do not satisfy
the convergence criterion can select a different best iterate from this tiny
seed perturbation; their largest final-frequency change is 0.291%. The direct
OSC statistics above are therefore the accuracy acceptance test for the new
default, rather than bit identity with the previous scan density.

## Source layout (release set; everything else lives in SISMO 1)

  src/sismo_precision.f90   kinds/constants
  src/sismo_types.f90       model / solution / result / config types
  src/sismo_config.f90      strict key/value configuration loading and validation
  src/sismo_io.f90          model reading (.osc.mod/sta1.d + sidecars), output writers
  src/sismo_poisson.f90     inhomogeneous Poisson: RK4 shooting + integral fallback
                             + the l=1 first-integral amplitude closure
  src/sismo_asymptotic.f90  asymptotic g-frequency (scan-window sizing only)
  src/sismo_core2.f90       the core: syst05 coefficients, scalar operator,
                             Sturm scan, split refinement, driver
  src/main.f90               configuration loading + driver call
  doc/second_order_cowling.md   derivation of the second-order core
  doc/core2_status.md           full development/validation log
  doc/release_reference_1M5.sismo   bit-identity reference output

## Known limits / next development

- Rotation: adding Coriolis to the inhomogeneous source is the project goal;
  note the l=1 closure must then use the ROTATING generalization of Takata's
  first integral.
- The tight inner loop (`inner_poisson_iterations > 1`) is present but not
  stabilized (needs Anderson + sit/search semantics); the release configuration
  uses the single-refresh loop, which is validated above.
- Mode labels are Sturm ranks from the window top (include f/p at the top);
  physical identification is by frequency.
