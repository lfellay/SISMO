# SISMO

**SISMO — Split Inhomogeneous Solver for Modelling Oscillations**

SISMO is research software for computing adiabatic stellar oscillation
frequencies. SISMO 2.0 locates the complete scalar Cowling spectrum with a
Sturm scan and refines each mode with an independently solved inhomogeneous
Poisson source. Asymptotic information sizes the scan window but is never used
as a mode-frequency seed.

The repository contains:

- `SISMO2.0/`: the current oscillation solver;
- `mesa2SISMO/`: conversion from a MESA profile to the intermediate model;
- `intSISMO/`: remeshing and conversion to the `.osc.mod` input format;
- `install.sh`: a source installer for all three programs.

## Build and install

Requirements are a POSIX shell, `make`, and GNU Fortran with OpenMP support.

```sh
git clone https://github.com/lfellay/SISMO.git
cd SISMO
./install.sh
export PATH="$PWD/bin:$PATH"
```

See [INSTALL.md](INSTALL.md) for the complete MESA-to-SISMO workflow and
configuration reference.

## Run

```sh
OMP_NUM_THREADS=8 sismo model.osc.mod results/model
```

Mode ranges, scan controls, Poisson choices, and convergence settings live in
`sismo.conf`. A different configuration can be supplied as the optional third
positional argument:

```sh
sismo model.osc.mod results/model /path/to/run.conf
```

## Validation

The SISMO 2.0 release configuration computes 216 ordered modes for the
24,992-point 1.5-solar-mass benchmark. The checked-in reference spectrum is
`SISMO2.0/doc/release_reference_1M5.sismo`; full accuracy and timing details
are recorded in [SISMO2.0/doc/RELEASE.md](SISMO2.0/doc/RELEASE.md).

Run the command-line tests with:

```sh
make -C SISMO2.0 test
```

Run the full numerical regression with the included benchmark model:

```sh
OMP_NUM_THREADS=2 make -C SISMO2.0 regression
```

The benchmark input is
`SISMO2.0/tests/data/release_1M5.osc.mod` (SHA-256
`fca1f0e6c00dc2624e5ea5664bc1db18f74f08b4a10660fb9498c5ec89b19316`).
Pass `MODEL=`, `REFERENCE=`, or `CONFIG=` to test other files.

## Documentation and citation

The method article is in `SISMO2.0/article/`. Citation metadata is provided in
[CITATION.cff](CITATION.cff).

The public software license is intentionally not assumed here. Until the
authors select and add a `LICENSE` file, normal copyright restrictions apply.
