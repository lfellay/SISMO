# intSISMO

`intSISMO` converts a `.madmod` stellar model into the adiabatic binary
structure read by SISMO and OSC. The normal input is the compact, versioned
adiabatic format written by `mesa2SISMO`. It writes an `.osc.mod` file and a
matching `.grid.d` sidecar; the two files should remain together under the
same basename.

Compact `SISMOAD2` models require this updated reader. Existing full
`.madmod` files remain readable, but an older `intSISMO` executable cannot
read a newly generated compact model.

The compact model stores only the zone table, stellar mass, radius and
gravitational constant, plus the profiles of `r`, `m`, `rho`, `P`,
`Gamma1`, and optional `brunt_N2`. If `brunt_N2` was not exported, it is
stored as zero and radial remeshing remains available.

## Build

The root installer builds the optimized program and places it in `bin`:

```sh
cd /path/to/SISMO
./install.sh
```

For a direct debug build with bounds and floating-point checks:

```sh
cd intSISMO
make clean
make debug=yes
```

For a direct optimized build, clean first so that objects compiled with the
debug flags cannot be reused:

```sh
cd intSISMO
make clean
make debug=no
```

## Run

```sh
intSISMO MODEL.madmod GRID_SIZE [GRID_STEP] [GRID_MODE] [--nonad]
```

`GRID_SIZE` must be between 8 and 2,000,000 points. `GRID_STEP` defaults to
`1`, uses the same upper safety limit, and adjusts the requested size to the
nearest compatible multiple. `GRID_MODE` defaults to `radial` and also
accepts `bv`.

Without a physics flag, `intSISMO` reads only the adiabatic structural fields.
This is the normal route.

For example:

```sh
intSISMO model.madmod 25000 16 radial
```

This writes `model.osc.mod` and `model.grid.d`. No external tables or runtime
data files are required. See `../INSTALL.md` for the complete MESA-to-SISMO
workflow.

For a full legacy `.madmod` created by `mesa2SISMO --nonad`, opt in
explicitly:

```sh
intSISMO model.madmod 25000 16 radial --nonad
```

This option requires and validates the full non-adiabatic input and retains
its atmosphere and zone structure during remeshing. The additional
thermodynamic fields are not written to `.osc.mod`; SISMO continues to
compute adiabatic oscillations.

The converter finishes both files under unique temporary names before
publishing them. Publication is serialized with a
`model.sismo-output.lock`; if a reported file or rename operation fails, the
converter attempts to restore the previous regular-file pair. An incomplete
rollback is reported and any unrestored backup is retained for manual
recovery. Existing directories, symbolic links, and special files are rejected
rather than replaced.

Because two separate files cannot be switched in one filesystem operation,
this is rollback-safe publication, not pair-wide atomic publication for
concurrent readers. An abrupt process kill, host crash, or power loss can
leave a lock, temporary file, or explicitly named `.bak` file. Before removing
a stale lock, first confirm that no `intSISMO` process is running and inspect
any retained backup.
