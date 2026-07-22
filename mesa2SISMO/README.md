# mesa2SISMO

This directory contains the MESA-side profile output setup and the
`mesa2SISMO` converter. The normal route is adiabatic: it needs only the
mechanical stellar structure and writes a compact, versioned `.madmod` file
read by `intSISMO`. The full legacy non-adiabatic format remains available as
an explicit option.

## Converter

Build the converter:

```sh
cd /path/to/SISMO/mesa2SISMO
make
```

Convert a MESA profile:

```sh
./mesa2SISMO /path/to/MESA-work/LOGS/profile1.data /path/to/models/model.madmod
```

No physics flag is needed. The required MESA columns are:

```text
mass or mass_grams
radius_cm
rho
pressure
gamma1
```

The output contains only the versioned adiabatic structure needed by
`intSISMO`. Use the matching updated `intSISMO` executable; older binaries
cannot read the `SISMOAD2` compact header. The older `--adiabatic` and
`--ad` options remain accepted as compatibility no-ops for existing scripts.
`brunt_N2` is optional and is retained only as an adiabatic grid-weighting
quantity for `intSISMO ... bv`; it defaults to zero when absent.

Then run `intSISMO` to remesh the model and write the SISMO/OSC input files:

```sh
cd /path/to/models
/path/to/SISMO/bin/intSISMO model.madmod 25000 16 radial
```

This creates `model.osc.mod` and the companion `model.grid.d`. Keep both
files together. `intSISMO` is the SISMO remesher and does not write a
`sta1.d` model.

## Adiabatic MESA profile

Use `profile_columns_sismo_adiabatic.list` as the MESA
`profile_columns_file`. It contains the small standard-column set required by
the default converter.

Install it in a MESA work directory:

```sh
cp /path/to/SISMO/mesa2SISMO/profile_columns_sismo_adiabatic.list \
   /path/to/MESA-work/profile_columns_sismo_adiabatic.list
```

Then set the profile column file in the MESA inlist:

```fortran
&star_job
   profile_columns_file = 'profile_columns_sismo_adiabatic.list'
/
```

The relative path is resolved from the MESA work directory, normally the
directory where `./rn` is executed.  An absolute path can also be used:

```fortran
&star_job
   profile_columns_file = '/path/to/SISMO/mesa2SISMO/profile_columns_sismo_adiabatic.list'
/
```

Profiles are written in `LOGS/profile*.data`; the rows are ordered from
surface to center.

## Optional full non-adiabatic output

Request the full legacy format explicitly:

```sh
./mesa2SISMO /path/to/MESA-work/LOGS/profile1.data \
  /path/to/models/model.madmod --nonad
```

For this route, use `profile_columns_mad_nonad.list` as the MESA
`profile_columns_file`. It contains the standard thermodynamic, opacity,
luminosity, nuclear, convection, composition, and atmosphere columns required
by the full conversion. The four custom EOS derivatives described below must
also be supplied through `run_star_extras`.

To use the resulting full input explicitly in `intSISMO`:

```sh
intSISMO model.madmod 25000 16 radial --nonad
```

`intSISMO` validates the full input and retains its atmosphere and zone
structure during remeshing. Its `.osc.mod` output contains only the
mechanical structure, and SISMO itself remains adiabatic.

### MESA atmosphere setup

Use an explicit MESA atmosphere boundary so the photospheric/raccord point is
well defined:

```fortran
&controls
   atm_option = 'T_tau'
   atm_T_tau_relation = 'Eddington'
   atm_T_tau_opacity = 'fixed'
/
```

For MESA's pulsation exports, these controls ask MESA to build atmospheric
points down to the requested outer optical depth:

```fortran
&controls
   atm_build_tau_outer = 1d-7
   atm_build_dlogtau = 0.025d0
   atm_build_errtol = 1d-8
/
```

In the full legacy format, the MAD radius is the first inward row with
`tau >= 2/3`. Rows above it are retained as the atmosphere; `mesa2SISMO`
never invents atmospheric layers. In MESA r26.4.1, `atm_build_*` affects
pulse-data atmospheres, not necessarily the normal `profile*.data` mesh.
Any atmosphere needed by the full non-adiabatic route must therefore be
present in the profile with all required quantities.

### Quantities requiring MESA custom output

These quantities are needed only with `--nonad`. The standard
`profile_columns_file` mechanism cannot create them, so export them through
`run_star_extras` extra profile columns:

```text
Cprho = (d ln Cp / d ln rho)_T
CpT   = (d ln Cp / d ln T)_rho
Qrho  = (d ln Q / d ln rho)_T, with Q = chiT/chiRho
QT    = (d ln Q / d ln T)_rho, with Q = chiT/chiRho
```

MESA's `how_many_extra_profile_columns` and
`data_for_extra_profile_columns` hooks append such columns to the profile
output.  They do not need to be listed in `profile_columns_mad_nonad.list`.

An example implementation is provided in `run_star_extras_mad_eos.f90`.
Inspect `src/run_star_extras.f90` before changing it. If it already contains
custom logic, do not replace it; merge the SISMO hooks and routines as
described below. If it is the unmodified MESA template and you intentionally
want to replace it, preserve a backup first:

```sh
cd /path/to/MESA-work
cp -p src/run_star_extras.f90 src/run_star_extras.f90.before_sismo
cp /path/to/SISMO/mesa2SISMO/run_star_extras_mad_eos.f90 src/run_star_extras.f90
./mk
./rn
```

To merge the MAD EOS pieces into an existing custom file:

1. Add these procedure pointers in `extras_controls`:

```fortran
s% how_many_extra_profile_columns => how_many_extra_profile_columns
s% data_for_extra_profile_columns => data_for_extra_profile_columns
```

2. Add the two routines `how_many_extra_profile_columns` and
   `data_for_extra_profile_columns` from `run_star_extras_mad_eos.f90`.
3. Add the helper `safe_log_deriv`.
4. Add the needed `use eos_def` and `use eos_lib` imports near the top of the
   module.

After running MESA, check the header of `LOGS/profile*.data`; the columns
`Cprho`, `CpT`, `Qrho`, and `QT` should appear after the standard profile
columns.

For MESA r26.4.1, AGSS09 opacity/table controls use these names:

```fortran
&star_job
   initial_zfracs = 6
/

&kap
   kap_file_prefix = 'a09'
   kap_CO_prefix = 'a09_co'
   kap_lowT_prefix = 'lowT_fa05_a09p'
/
```

For a newly created pre-main-sequence or initial model, select the desired
hydrogen fraction \(X\) indirectly through the helium fraction:

```fortran
&controls
   initial_z = 0.02d0
   initial_y = 0.28d0   ! X is then 1 - initial_y - initial_z
/
```

MESA r26.4.1 has no `initial_x` control. `initial_y` does not alter the
composition when loading a saved model or a ZAMS model; use the appropriate
MESA composition-change/relaxation controls for those workflows.

Do not add `set_uniform_initial_composition = .true.` by itself. If uniform
composition is deliberately required, first choose \(X\), \(Y\), and \(Z\)
for the active nuclear network, then explicitly provide the compatible
`initial_h1`, `initial_h2`, `initial_he3`, and `initial_he4` fractions in
`&star_job`, together with `initial_z` in `&controls` and the desired
`initial_zfracs`. Verify that all isotope fractions are non-negative and sum
to one before running MESA. For ordinary pre-main-sequence creation,
`initial_y`, `initial_z`, and `initial_zfracs` are the safer workflow.

For opacity, the standard profile list provides face derivatives:

```text
dkap_dlnrho_face = (d ln kappa / d ln rho)_T at the outer cell face
dkap_dlnT_face   = (d ln kappa / d ln T)_rho at the outer cell face
```

These can be interpolated to MAD's mesh.  Centered opacity derivatives would
be better if they are easy to expose from the MESA opacity call.

The isotope columns `h1`, `he3`, `he4`, `c12`, `n14`, and `o16` are diagnostic
composition columns.  If a chosen MESA network does not contain one of these
isotopes and MESA rejects the profile list, comment out the missing isotope;
the future converter can still use the aggregate `x`, `y`, and `z` columns.
