# The Cowling oscillation equations as a single second-order ODE

Basis for the new core matrix inversion: reduce the adiabatic Cowling mechanics
(radial momentum + horizontal momentum + continuity + adiabatic relation) to ONE
exact second-order differential equation, whose discretization is a scalar
tridiagonal (3-point) operator instead of the current 4-variable block system.
The inhomogeneous terms (phi' now, Coriolis later) enter as a right-hand side.

## 1. Starting point: the two first-order Cowling equations

Variables: xi_r (radial displacement), p' (Eulerian pressure perturbation).
The horizontal momentum equation is algebraic in the Cowling approximation and
defines xi_h; continuity + adiabatic eliminate rho'.  Standard form (Unno et al.
1989, eqs. 15.18-15.19 with Phi' = 0):

    d xi_r/dr = ( g/c^2 - 2/r ) xi_r  +  (1/(rho c^2)) ( L_l^2/sigma^2 - 1 ) p'      (1)
    d p'/dr   = rho ( sigma^2 - N^2 ) xi_r  -  ( g/c^2 ) p'                          (2)

with
    c^2    = Gamma1 P / rho                    (adiabatic sound speed)
    L_l^2  = l(l+1) c^2 / r^2                  (Lamb frequency)
    N^2    = g ( (1/Gamma1) dlnP/dr - dlnrho/dr )   (Brunt-Vaisala)
    xi_h   = p' / (rho sigma^2 r)              (horizontal displacement, Cowling)

Shorthand for the reduction:

    d xi_r/dr = A xi_r + B p'        A = g/c^2 - 2/r
                                     B = (L_l^2/sigma^2 - 1) / (rho c^2)
    d p'/dr   = C xi_r + D p'        C = rho (sigma^2 - N^2)
                                     D = -g/c^2

Note A + D = -2/r (the g/c^2 terms cancel) — used below.

## 2. Exact second-order equation for p'

From (2): xi_r = ( p'' - D p' ) / C   [ ' = d/dr ].  Substituting into (1):

    p''' - (A + D + C'/C) p'' + (AD - BC + (C'/C) D - D') p' = 0

i.e. with the shorthands evaluated:

    p''  +  ( 2/r - C'/C ) p'  +  Q(r) p  = 0                                        (3)

    Q(r) = K(r) + 2g/(r c^2) - g^2/c^4 + (g/c^2)' - (g/c^2) (ln|C|)'

    K(r) = (sigma^2 - N^2)(sigma^2 - L_l^2) / (sigma^2 c^2)

(K is the classical dispersion term: propagation where K > 0, i.e. sigma^2 above
both N^2 and L_l^2 -> p modes, or below both -> g modes.  The remaining terms in
Q are the exact acoustic-cutoff / stratification corrections; no approximation
has been made.)

## 3. Self-adjoint (Sturm–Liouville) form

The integrating factor mu with mu'/mu = 2/r - C'/C is mu = r^2 / C, giving the
exact self-adjoint form

    d/dr [  r^2/(rho (sigma^2 - N^2))  dp'/dr  ]
        +  r^2/(rho (sigma^2 - N^2)) * Q(r) * p'  =  0                               (4)

with the algebraic back-substitutions

    xi_r = [ dp'/dr + (g/c^2) p' ] / ( rho (sigma^2 - N^2) )
    xi_h = p' / (rho sigma^2 r)

This is a lambda-nonlinear Sturm–Liouville problem (sigma^2 appears in the
coefficients, not merely as a linear eigenvalue).

Turning points: the p'-form is singular where sigma^2 = N^2 (C = 0), i.e. at the
g-cavity boundaries.  The twin reduction in xi_r (eliminate p' via (1)) is

    xi_r'' + (2/r - B'/B) xi_r' + [ K + 2g/(rc^2) - g^2/c^4
             + (B'/B)(g/c^2 - 2/r) - (g/c^2 - 2/r)' ] xi_r = 0                       (5)

singular instead where sigma^2 = L_l^2.  Every classical single-variable form
carries one of the two turning-point families; the choice of working variable
(or a hybrid/canonical variable) is an implementation decision for the
discretization step.

## 4. Dimensionless (code) variables

With x = r/R and the intSISMO/OSC structure arrays (GM = R = 1 units):

    s2 = qx3 = m(x)/x^3      ->  g_hat      = s2 * x
    s3 = P/rho, s5 = Gamma1  ->  c^2_hat    = s5 * s3
    s4 = rho
    s6 = aosc                ->  N^2_hat    = s2 * s6 * x^2
                                 (from brunt_integral: sqrt(s2*|s6|) = N/x)
    L_l^2_hat = l(l+1) * s5 * s3 / x^2
    sigma^2   = omega^2 (the code's dimensionless lambda)

The exact syst05 scalings of the dependent variables (a, pp) must be pinned to
these when the operator is assembled, so results remain comparable with the
existing block solver and with OSC.

## 5. The inhomogeneous terms (for the split architecture)

Keeping Phi' (and later Coriolis) as a FROZEN source, the first-order pair gains
right-hand sides (Unno 15.18-15.19):

    d xi_r/dr = A xi_r + B p' + R1        R1 = l(l+1) Phi' / (sigma^2 r^2)
    d p'/dr   = C xi_r + D p' + R2        R2 = -rho dPhi'/dr
    xi_h      = ( p'/rho + Phi' ) / (sigma^2 r)

Carrying R1, R2 through the same elimination, the second-order equation (3)
becomes inhomogeneous with the exact right-hand side

    RHS(r) = R2' - ( C'/C + A ) R2 + C R1                                            (6)

and in self-adjoint form:  d/dr[ mu p'' ] + mu Q p' = mu * RHS.

So the split structure is preserved verbatim: the homogeneous operator (4) is
the Cowling mechanics; the frozen potential (and any added physics) drives it
through (6).  The core inversion becomes a scalar tridiagonal solve.

## 6. What this buys

- 3-point scalar stencil (tridiagonal) instead of 4x4 block-tridiagonal:
  cheaper factorization, trivially parallel scans.
- Self-adjoint discretization of (4) -> symmetric matrix -> the pivot signs of
  the LDL^T factorization give a rigorous Sturm count of the COWLING spectrum
  (mode counting/labeling without eigenfunctions), while the physics stays
  split: the inhomogeneous sources never enter the operator.
- The single-variable form is the natural place to inject Coriolis terms as
  sources on the RHS, mirroring the Phi' treatment.

Open implementation choices (next steps): working variable (p' vs xi_r vs a
canonical/hybrid variable) w.r.t. the turning-point singularities; boundary
conditions (regular center: p' ~ r^l; surface: delta p = 0 or atmospheric
match); discretization of (4) on the existing staggered midpoint grid.
