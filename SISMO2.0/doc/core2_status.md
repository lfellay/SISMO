# SISMO 2.0 core — implementation status (2026-07-04)

## What is implemented (src/sismo_core2.f90, flag --core2)

1. Second-order scalar reduction of the syst05 MECHANICAL pair (rows 1,2 of
   system05_coefficients), in the code's own variables (y1, y2; y3 = scaled
   frozen potential as source).  NOTE the row-1 derivative weight: the syst05
   scheme has omega^2 * y1' = ... (diag1(1,1)=1 in difference_block), so
   a12 = d0(1,2)/omega^2 + d1(1,2) -- this term carries the g-branch
   dispersion k^2 ~ N^2/omega^2 (first implementation missed it -> sparse
   wrong spectrum; decisive symptom to remember).
2. Exponentially-fitted staggered box scheme: y1 at nodes, y2 eliminated at
   midpoint faces via cell-local integrating factors (a GLOBAL canonical
   transform is exact but overflows across ~200 pressure scale heights).
   Result: scalar tridiagonal with positive off-diagonal products
   (symmetrizable -> Sturm pivot counting).  Lamb-crossing faces (a12=0) are
   apparent singularities: entries grow ~1/a12 and enforce the exact
   degenerate constraint; LU handles them; exact zero guarded.
3. Boundary rows = the block solver's own (center_block mechanical row per
   case; surface_block row 1 = y1 - y2 + y3 = 0 i.e. delta-p = 0).  Centre
   guard skips nodes with x <= 1e-8 (1/x coefficients).
4. STURM SCAN: sweep + bisection on the parity of
   (negative LU pivots + negative face denominators) -- each face-denominator
   (Lamb) crossing flips a pivot AND the denominator count, so the parity is
   invariant under the pole artifacts and flips exactly at eigenvalues.
   Complete labeled spectrum (g + f + p) on the FULL 25k grid in ~2 s.
5. Split (inhomogeneous) refinement: seed = Cowling eigenvalue; iterate
   phi <- solve_poisson_correction(mechanics) (UNCHANGED Poisson solver);
   rhs(y3); forced-response + Newton frequency update with finite-difference
   T'; clamp 1e-3*w2/step.  Currently a SINGLE phi refresh per step (loose).

## Validation so far (1M5cmp.osc.mod, full grid)

- Cowling l=2 vs legacy block solver (its converged modes g8-g12):
  0.000-0.020% agreement.  Legacy's unconverged modes (conv=F) differ 0.1-3%;
  core2 values are exact eigenvalues of the discrete operator.
- The scan resolves the complete spectrum: Cowling g1=3.652, g2=3.168,
  g3=2.463 (~+1-2% above the full values = the expected Cowling shift),
  f=4.05, p-branch above -- structure correct, no gaps/doublets.
- Split l=2 (single-refresh): g1 3.652->3.587 vs OSC 3.593 (0.19%),
  g2->3.080 (0.44%), mid orders 0.08-0.14% vs GYRE, f-mode 4.05->3.934 vs
  OSC 3.953 (0.48%).  First pipeline that reaches the f-mode.
- Split l=1: the known dipole loose-coupling drift (~1% low orders; the
  Cowling->full travel for l=1 tops is large and the loose loop stalls).

## To do

- Port the tight-coupling machinery into core2_split_refine (inner phi
  refreshes to convergence per frequency step, poisson_relax, Anderson,
  rate-corrected inner criterion) -- the l=1 fix, same as SISMO 1.
- Convergence polish: several l=2 modes hit max_iter with the frequency
  already right (the phi-change criterion stalls at ~1e-4 with the loose
  loop); expected to resolve with the inner loop.
- Label calibration (Sturm index vs OSC/GYRE n_pg), full l=1+l=2 comparison
  figures vs OSC/GYRE, performance profile, Takata-source hook for l=1.

## Full comparison + timing (2026-07-04, l=1,2 n=1..100, 1M5cmp model, this machine)

Accuracy (nearest-frequency; figures 1M5/compare_core2_gyre.png, compare_core2_osc.png;
run file 1M5/sismo2_1M5_gmodes.sismo):
  core2 vs GYRE : l=1 0.342% (87% <1%), l=2 0.223% (97% <1%)
  core2 vs OSC  : l=1 0.089% (max 2.2%), l=2 0.0020% (95% <0.1%, max 1.27%)
  (references: SISMO1 canonical vs GYRE 0.253/0.227%; vs OSC 0.048/0.0001%.
   core2 l=2 vs OSC shows the clean discretization curve 1e-4 -> 1e-2 with |n|;
   l=1 carries the loose-refresh drift -- tight-coupling port pending.
   core2 max errors BEAT SISMO1 (l=1 2.2 vs 3.1%, l=2 1.27 vs 9.1%) because the
   f-adjacent tops and the f-mode are genuinely computed, not nearest-matched.)

Computing time (same model, same machine, except GYRE):
  OSC  (scang, oscQI=2, ~16 threads) : 76 s wall   1224 s CPU   202 modes
  SISMO 1 canonical (8 threads)     : ~95 s wall   ~560 s CPU   204 modes
  SISMO 2.0 core2 (SINGLE thread)   : 104 s wall    104 s CPU   216 modes (incl. f/p)
  GYRE (user's run, own grid/machine): ~9 s (sum of scan timers in its log; indicative)
  -> core2 uses ~5x less CPU than SISMO1 and ~12x less than OSC, before any
     parallelization (its mode loop is trivially parallel); the Sturm scan replaces
     all seeding machinery and costs ~2 s per l on the full 25k grid.

## After the tight-coupling port + parallel mode loop (2026-07-04, second full run)

The port added: inner phi-refresh loop (rate-corrected criterion, relax, factor
reuse), NaN guard, incremental shape update, and an OpenMP-parallel mode loop.
The tight inner loop is NOT yet stable for the strongly-shifted modes (dipole
low orders, f-mode): it needs the Anderson + sit/search machinery (next task).
Production configuration for now: --inner-poisson-iters 1 --poisson-relax 1.0
(single refresh + the per-step shape correction of the ported loop).

Full run l=1,2 n=1..108 (core2_v2 = 1M5/sismo2_1M5_gmodes.sismo):
  vs GYRE: l=1 0.266% (95% <1%), l=2 0.208% (100% <1%)   [OSC ref: 0.218/0.216%]
  vs OSC : l=1 0.0081% (88% <0.1%, max 2.6%), l=2 0.0076% (92% <0.1%, max 5.3%)
  g-branch tops now excellent (l=2 g1 0.03%, g2 0.14%, g3 0.08% vs OSC);
  residual outliers = the f-mode (l=2, 5.3%) and the two f-adjacent dipole tops
  (2-3%) -- exactly the strong-self-gravity modes awaiting the stabilized tight loop.

Timing (same model/machine):
  OSC  (~16 threads)      : 76 s wall   1224 s CPU   202 modes
  SISMO 1 (8 threads)    : ~95 s wall   ~560 s CPU   204 modes
  SISMO 2.0 (8 threads)  : 39.5 s wall   172 s CPU   216 modes (incl. f/p)
  GYRE (own log, indicative): ~9 s
  -> SISMO 2.0 is now the fastest of the three on this model in wall time
     (2.4x vs SISMO 1, 1.9x vs OSC) and CPU (3.3x vs SISMO 1, 7x vs OSC).

## Takata dipole source adopted for l=1 (2026-07-04, core2_v3 = canonical)

Canonical command adds --takata-source --takata-blend 0.8 (l=1-gated in the
Poisson module; l=2 bit-unchanged).  Effect: moves the dipole lambda-update's
zero toward the truth, so the Cowling->full travel LANDS instead of stalling:
the catastrophic f-adjacent mode went 10% -> 2.7%, and the low branch now
matches or beats SISMO 1 mode-by-mode.
  v3 vs GYRE: l=1 0.222% (98% <1%)  l=2 0.208% (100% <1%)
              [SISMO 1: 0.253/0.227%; OSC-vs-GYRE floor: 0.218/0.216%]
  v3 vs OSC : l=1 0.0081% (max 2.3%)  l=2 0.0076% (max 5.3%)
  36 s wall / 133 s CPU @ 8 threads.
SISMO 2.0 now beats SISMO 1 on BOTH degrees and is ~2.5x faster in wall time.
Remaining >1% vs OSC: the f-mode and two f-adjacent dipole tops (2-3%) -- the
stabilized tight loop (Anderson + sit/search) or the full Takata reduction.

## Sign-alignment fix (2026-07-04, core2_v4 = canonical): l=2 residuals SOLVED

The l=2 errors near n~-10 (and the f-mode 5.3%) were NOT physics: the refine's
incremental update could land with negative overlap, and after normalization the
iterate FLIPPED SIGN every pass (phi-change oscillating ~2.0 = the flip
signature), so (phi, shape) never converged; the mode stalled at its seed or in
a limit cycle.  Converged modes were exact (0.00%); stuck ones showed 0.1-5%.
FIX: sign-align the updated eigenvector with the previous iterate (the sign is
arbitrary) in both the inner refresh and the lambda step + best-iterate return
+ poisson_relax 0.4.  Result (core2_v4, canonical command now
--inner-poisson-iters 1 --poisson-relax 0.4 --takata-source --takata-blend 0.8):
  vs OSC : l=2 0.0018% median, 99% <0.1%, MAX 0.32% over all 100 modes
           (f-mode 3.9502 vs 3.9529 = 0.07%);  l=1 0.0099%, max 3.0%
  vs GYRE: l=2 0.217% (100% <1%) = the re-mesh floor (OSC itself: 0.216%);
           l=1 0.221% (97% <1%) vs OSC's 0.218%.
  24 s wall / 57 s CPU @ 8 threads (most modes now converge in ~20 iters).
Remaining >1% anywhere: only the two f-adjacent dipole tops (l=1, 2-3%).

## The l=1 inverse power law vs OSC (2026-07-04, core2_v5 = canonical; best-gate 1e-1)

Measured: l=1 error vs OSC ~ n^-2.33 ; the PHYSICAL dipole self-gravity
correction (Cowling seed vs full value) ~ n^-2.39 ; ratio err/shift = 1.1-1.45,
nearly constant over n=1..30.  Interpretation: the power law is NORMAL -- it is
the textbook decay of the dipole self-gravity coupling with radial order
(rapidly oscillating rho' self-cancels in the Poisson integral, ~n^-2).  The
l=1 error TRACKS it at ratio ~1.3 because the dipole refine still fails to
apply most of the correction (no stationary point for the frozen-phi dipole
update; the takata blend stabilizes but the travel does not complete) -- the
quantitative signature of the known l=1 frontier, NOT a numerical bug.
Proof by contrast: l=2, where the split iteration converges, sits at 0.0018%
median = THREE orders below its own Cowling shift.
v5 numbers: vs OSC l=2 max 0.10% (100% <0.1%!), l=1 max 2.3%;
vs GYRE l=1 0.221%/97%, l=2 0.217%/100%.  25 s wall @ 8 threads.
When the dipole frontier work lands (stabilized tight loop or the Takata
reduced system), this power-law envelope should collapse to the l=2-like floor.

## FIRST-INTEGRAL AMPLITUDE CLOSURE: the dipole SOLVED (2026-07-04, core2_v6 = FINAL canonical)

--takata-closure: for l=1, the shooting Poisson amplitude A (phi = phi_p + A*phi_h)
is set by least-squares Takata J=0 over the star (momentum conservation) instead
of the vacuum surface BC.  J is linear in A through BOTH the phi terms and the
pp term (pack relation -P*pp = -rho*sp*s2*(y2-y3)); A = -<H,M>/<H,H>, dx-weighted,
computed inside solve_poisson_shooting_correction from ymech/particular/homogeneous.
Parameter-free (replaces --takata-source/--takata-blend), fully split.

RESULT: the l=1 power-law envelope COLLAPSED to the discretization floor:
  vs OSC : l=1 0.0019% median, 99% <0.1%, MAX 0.21% over all 98 modes
           (the f-adjacent dipole: 2.7378 vs 2.7325 = 0.20%, was 2.7-10%);
           l=2 0.0018%, 100% <0.1%, max 0.10%.
  vs GYRE: l=1 0.216% (100% <1%) -- BELOW OSC's own 0.218%; l=2 0.217% (100% <1%).
  EVERY mode of the l=1,2 n=1..100 problem is <0.21% vs OSC.  35 s wall @ 8 thr.
The validated settings are now selected through the canonical configuration:

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
tolerance = 1e-8
inner_poisson_iterations = 1
poisson_relaxation = 0.4
poisson_inner_tolerance = 0.0
```

Current invocation: `sismo MODEL.osc.mod OUT [CONFIG]`. If CONFIG is omitted,
lookup proceeds through `SISMO_CONFIG`, `./sismo.conf`, then the configuration
beside the executable.
The blend (--takata-source) and all l=1 special-casing are now obsolete.

## Timing at 18 threads (full machine; same host, 2026-07-04)

  code                     wall      CPU     modes   grid
  GYRE 8.1 (MESA build)    9.4 s     71 s    200     1131-pt MESA profile + internal remesh
  OSC (scang, oscQI=2)    72.5 s   1095 s    202     25k intSISMO + per-mode WKB
  SISMO 1 (canonical)    72.7 s    905 s    204     25k intSISMO + multigrid/WKB
  SISMO 2.0 (v6)         29.4 s    159 s    216     25k intSISMO, full grid, no remesh

  - On the SAME 25k model, SISMO 2.0 is 2.5x faster than OSC and SISMO 1 in
    wall time and 6-7x cheaper in CPU, while also computing the f/p modes.
  - GYRE's wall time is on its own (much smaller) input grid with internal
    adaptive remeshing -- not directly comparable, listed for context.
  - SISMO 2.0 scaling 8->18 threads: 35->29 s; the serial Sturm scan (~10 s)
    now dominates (Amdahl); parallelizing the scan is the next speed lever.

## Parallel Sturm scan (2026-07-04 final): 9.1 s wall at 18 threads

The sweep points and the per-crossing bisections are independent factorizations
-> both OpenMP-parallel via a thread-safe sturm_count_at (local workspace).
Spectra verified BIT-IDENTICAL to the serial scan.  Final 18-thread table
(same host; l=1,2 n=1..100):

  code                     wall      CPU     modes   grid
  GYRE 8.1                 9.4 s     71 s    200     1131-pt MESA profile
  OSC (scang, oscQI=2)    72.5 s   1095 s    202     25k intSISMO
  SISMO 1 (canonical)    72.7 s    905 s    204     25k intSISMO
  SISMO 2.0 (final)       9.1 s    147 s    216     25k intSISMO (full grid)

  SISMO 2.0 now matches GYRE's wall time on a 22x larger input grid, and is
  8x faster than OSC / SISMO 1 on the same model -- while agreeing with OSC
  to max 0.21% (l=1) / 0.10% (l=2) over every mode, fully split/inhomogeneous.
