# SISMO 2.0

**SISMO — Split Inhomogeneous Solver for Modelling Oscillations**

Install SISMO 2.0 together with `mesa2SISMO` and `intSISMO` by running
`../install.sh` from this directory. See `../INSTALL.md` for prerequisites and
the complete conversion workflow.

Adiabatic stellar oscillation code: complete labeled spectra (g, f, p modes)
with a SPLIT (inhomogeneous) Poisson treatment of self-gravity on a scalar
second-order Cowling core.  Sturm-scan mode location (no seeding, no
asymptotics), OpenMP-parallel, full model grid.

    make
    export OMP_NUM_THREADS=<cores>
    ./sismo MODEL.osc.mod OUT

SISMO reads its numerical and physics settings from the plain-text
`sismo.conf` file. The shipped file contains the validated release
configuration:

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

To select another file explicitly, pass it as the optional third positional
argument:

    ./sismo MODEL.osc.mod OUT /path/to/custom.conf

Configuration lookup follows this order: the explicit positional file,
`SISMO_CONFIG`, `./sismo.conf`, then `sismo.conf` beside the executable. SISMO
stops with a clear error if none of these files can be read. Command-line
numerical flags are no longer accepted.

Each non-comment line is `key = value`; `#` and `!` introduce comments. Boolean
values should normally be written as `true` or `false`. Unknown or duplicate
keys are errors.

`scan_points` is the endpoint-inclusive number of initial Sturm samples for
each degree. SISMO stops with a diagnostic instead of silently skipping modes
if the requested sampling is too coarse. Set `use_poisson = false` for a pure
Cowling spectrum. Full-grid eigenfunctions are omitted by default; set
`write_eigenfunctions = true` to write one mechanical `.eig` file per mode.
Potential-related columns in these optional files are currently zero.

See doc/RELEASE.md for the validated configuration, method summary, accuracy
(max 0.21% vs OSC over every l=1,2 mode of the 1M5 benchmark) and timing
(about 1.9 s wall / 28.0 s CPU for 216 modes on a 25k-point model at 18 threads).
Design note: doc/second_order_cowling.md.  Development log: doc/core2_status.md.
SISMO 1 (the block-solver lineage and all legacy machinery) lives in ../SISMO.
