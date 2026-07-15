# Installing SISMO 2.0

The installer builds the three SISMO programs and copies them, together with
the SISMO configuration, into one binary directory.

## Requirements

- A Unix-like system with a POSIX-compatible shell
- `make`
- GNU Fortran (`gfortran`) with OpenMP support
- The standard `install` utility

MESA is not required to compile SISMO. It is only needed to generate a stellar
profile for `mesa2SISMO`.

## Install

From the repository root:

```sh
./install.sh
```

By default, this installs:

```text
bin/sismo
bin/mesa2SISMO
bin/intSISMO
bin/sismo.conf
bin/sismo.conf.default
```

To install elsewhere:

```sh
./install.sh --bin-dir "$HOME/.local/bin"
```

If GNU Fortran has a different name or location:

```sh
FC=/path/to/gfortran ./install.sh
```

The destination must be writable. Existing executables with the same names
are replaced only after all three builds and command-line checks succeed.
`sismo.conf.default` is refreshed on every installation. `sismo.conf` is
created only when it does not already exist, so reinstalling SISMO preserves
your settings.

To use the default installation from any directory:

```sh
export SISMO_ROOT=/absolute/path/to/SISMO
export PATH="$SISMO_ROOT/bin:$PATH"
```

Add the `PATH` line to your shell configuration if you want it to persist.

## Verify the installation

```sh
sismo --help
mesa2SISMO --help
intSISMO --help
```

## Configure SISMO

SISMO reads its mode selection, physics choices, and numerical controls from a
plain-text `key = value` file. The installed default is:

```text
l_min = 1
l_max = 2
g_min = 1
g_max = 108

use_poisson = true
takata_closure = true
write_eigenfunctions = false

scan_points = 2000
eps_g = 2.5
max_iterations = 150
tolerance = 1.0e-8
inner_poisson_iterations = 1
poisson_relaxation = 0.4
poisson_inner_tolerance = 0.0
```

The keys have the following meanings:

- `l_min`, `l_max`: inclusive spherical-degree range.
- `g_min`, `g_max`: inclusive radial-order range selected from the scan.
- `use_poisson`: enable the split Poisson refinement; set it to `false` for
  the Cowling approximation.
- `takata_closure`: enable the dipole first-integral closure for `l = 1`.
- `write_eigenfunctions`: write one full-grid mechanical `.eig` file per mode.
- `scan_points`: endpoint-inclusive Sturm samples used for each degree.
- `eps_g`: phase used to set the scan-frequency window.
- `max_iterations`, `tolerance`: outer refinement limit and relative frequency
  convergence tolerance.
- `inner_poisson_iterations`, `poisson_relaxation`,
  `poisson_inner_tolerance`: controls for the frozen-potential refresh.

Boolean values are written as `true` and `false`. The parser also accepts
`.true.`/`.false.`, `yes`/`no`, `on`/`off`, and `1`/`0`, without regard to
case. Blank lines are allowed, and `#` or `!` begins a comment, including an
inline comment. The supplied numerical values are the validated SISMO 2.0
settings and normally do not need to be changed.

SISMO uses the first configuration it finds in this order:

1. The optional third command-line argument.
2. The file named by the `SISMO_CONFIG` environment variable.
3. `sismo.conf` in the current working directory.
4. `sismo.conf` beside the `sismo` executable.

SISMO stops with an error if it cannot open any configuration file. The
installed `sismo.conf` therefore provides the normal fallback. To create a
project-specific configuration, copy the refreshed template into the working
directory and edit the copy:

```sh
cp "$SISMO_ROOT/bin/sismo.conf.default" ./sismo.conf
```

For a configuration stored elsewhere, either pass it explicitly:

```sh
sismo model.osc.mod results/model /path/to/run.conf
```

or set the environment variable:

```sh
SISMO_CONFIG=/path/to/run.conf sismo model.osc.mod results/model
```

The explicit third argument has priority over every automatic location. To
change the fallback for every run, edit the installed `bin/sismo.conf`. To
restore the current release defaults, copy `bin/sismo.conf.default` over it.

## Complete workflow

Create a working directory first:

```sh
mkdir -p /path/to/work/results
cd /path/to/work
```

### 1. Convert a MESA profile

```sh
mesa2SISMO /path/to/profile.data model.madmod --adiabatic
```

If the output name is omitted, `mesa2SISMO` derives it from the input name.
For example, `profile.data` becomes `profile.madmod`.

### 2. Create the SISMO input model

```sh
intSISMO model.madmod 25000 16 radial
```

The arguments are:

```text
intSISMO MODEL GRID_SIZE [GRID_STEP] [GRID_MODE]
```

- `GRID_SIZE` is the requested number of output points and must be at least 8.
- `GRID_STEP` defaults to `1`.
- `GRID_MODE` defaults to `radial`; `bv` is also available.

The example requests 25,000 points with a step of 16. `intSISMO` selects the
nearest compatible size, 24,992 points. It creates:

```text
model.osc.mod
model.grid.d
```

Keep these two files together with the same basename. SISMO reads
`model.grid.d` to recover the grid-step information.

### 3. Run SISMO 2.0

The command syntax is:

```text
sismo MODEL OUT [CONFIG]
```

`MODEL` is the input model, `OUT` is the output basename, and `CONFIG` is an
optional configuration path for this run.

```sh
OMP_NUM_THREADS=8 sismo model.osc.mod results/model
```

Change `OMP_NUM_THREADS` to the number of CPU threads you want SISMO to use.
The command uses the first configuration found according to the lookup order
above. To use a particular file for this run, add it as the third positional
argument:

```sh
OMP_NUM_THREADS=8 sismo model.osc.mod results/model /path/to/run.conf
```

SISMO reports when `scan_points` is too small to isolate every crossing. Set
`write_eigenfunctions = true` only when you need one full-grid mechanical
`.eig` file per mode; potential-related columns in these optional files are
zero. The parent directory of the output base must already exist. The command
creates:

```text
results/model.sismo
results/model.comparison
```

## Generated files

- `mesa2SISMO` writes a `.madmod` stellar model.
- `intSISMO` writes an `.osc.mod` model and its `.grid.d` sidecar.
- `sismo` writes the mode spectrum to `.sismo` and comparison information to
  `.comparison`.
- With `write_eigenfunctions = true`, `sismo` also writes one mechanical
  `.eig` file per mode.

These files are written in the working directory or at the paths given on the
command line, not in the installation directory.

## Troubleshooting

### `gfortran: command not found`

Install GNU Fortran and confirm that `gfortran --version` works. If the
compiler has another name, set `FC` when running the installer.

### OpenMP or Fortran runtime library error

Build with a GNU Fortran installation that includes OpenMP. The installed
executables depend on compatible GNU Fortran and OpenMP runtime libraries.

### Permission denied for the installation directory

Use a writable location:

```sh
./install.sh --bin-dir "$HOME/.local/bin"
```

### SISMO cannot open its output

Create the output directory before running it:

```sh
mkdir -p results
```

### SISMO cannot find a configuration

Confirm that the installed active configuration exists:

```sh
ls "$SISMO_ROOT/bin/sismo.conf"
```

For a local configuration, copy the installed template:

```sh
cp "$SISMO_ROOT/bin/sismo.conf.default" ./sismo.conf
```

You can also pass the configuration as the third argument or set
`SISMO_CONFIG`. If none is available, SISMO reports all the ways to select one.

### SISMO does not recover the intended grid step

Keep the `.grid.d` file beside the corresponding `.osc.mod` file and preserve
their common basename.

### Build artifacts came from another compiler or computer

The installer cleans and rebuilds all three programs. Run it again:

```sh
./install.sh
```

## Portability

The installer builds executables for the current operating system, CPU
architecture, GNU Fortran runtime, and OpenMP runtime. When moving SISMO to
another computer or operating system, copy the source tree and run the
installer there instead of copying the existing `bin` directory.
