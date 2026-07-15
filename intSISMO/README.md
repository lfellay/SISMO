# intSISMO

`intSISMO` converts a `.madmod` stellar model into the binary structure read by
SISMO and OSC. It writes an `.osc.mod` file and a matching `.grid.d` sidecar;
the two files should remain together under the same basename.

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

For a direct optimized build, use `make debug=no`.

## Run

```sh
intSISMO MODEL.madmod GRID_SIZE [GRID_STEP] [GRID_MODE]
```

`GRID_STEP` defaults to `1`. `GRID_MODE` defaults to `radial` and also accepts
`bv`. The requested grid size is adjusted to the nearest multiple of
`GRID_STEP`.

For example:

```sh
intSISMO model.madmod 25000 16 radial
```

This writes `model.osc.mod` and `model.grid.d`. No external tables or runtime
data files are required. See `../INSTALL.md` for the complete MESA-to-SISMO
workflow.
