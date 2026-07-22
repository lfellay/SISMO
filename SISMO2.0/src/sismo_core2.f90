! =============================================================================
! SISMO 2.0 core: the Cowling mechanics as ONE second-order scalar equation.
!
! The syst05 mechanical pair (rows 1,2 of system05_coefficients; y3 = the scaled
! frozen potential as inhomogeneous source) is
!     y1' = a11 y1 + a12 y2 + a13 y3
!     y2' = a21 y1 + a22 y2 + a23 y3
! with a11 = w2*d1(1,1), a12 = d0(1,2)+w2*d1(1,2), a13 = w2*d1(1,3),
!      a21 = d0(2,1)+w2*d1(2,1), a22 = d0(2,2), a23 = d0(2,3)   (all real).
!
! y2 is eliminated on the staggered grid (y1 at the model nodes, y2 at the
! midpoint faces) with CELL-LOCAL exponential integrating factors (a global
! integrating-factor transform is exact analytically but spans the star's ~200
! pressure scale heights and overflows).  Face relation across [x_i, x_{i+1}]:
!     y1_{i+1} = e^p y1_i + h*phi1(p) * ( a12_f y2_f + a13_f y3_f ),  p = a11_f h
! Node relation across the dual cell [f-, f+]:
!     y2_{f+}  = e^q y2_{f-} + V*phi1(q) * ( a21_i y1_i + a23_i y3_i ), q = a22_i V
! (phi1(z) = (e^z - 1)/z).  Substituting the face relations into the node
! relation gives a SCALAR TRIDIAGONAL system in y1 whose off-diagonal products
! are strictly positive (both couplings across a face share the same factor
! 1/(h phi1 a12), so their product is a positive square times exponentials).
! A tridiagonal with positive off-diagonal products is symmetrizable, so the
! LU pivot signs form a rigorous Sturm sequence: the negative-pivot count
! changes by one at every eigenvalue -> complete, labeled Cowling spectrum by
! sweep + bisection, with NO eigenfunctions involved.
!
! The turning point (a12 = 0, sigma^2 = Lamb) is an apparent singularity: the
! face coefficients ~ 1/a12 grow there, which enforces the exact degenerate
! constraint y1' = a11 y1 + a13 y3; the physical solution is smooth through it
! and the LU handles the stiff row (guarded against exact zero).
!
! The split (inhomogeneous) architecture is preserved verbatim: the operator is
! pure Cowling mechanics; the frozen potential y3 = phi_code/(x^l qx3) from the
! existing Poisson solver drives the right-hand side, exactly as in the block
! solver.  Boundary rows are the block solver's own (center_block row 1 /
! surface_block row 1: y1 - y2 + y3 = 0, i.e. delta-p = 0).
! =============================================================================
module sismo_core2
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use sismo_precision, only : dp, tiny, czero
  use sismo_types, only : stellar_model, mode_solution, mode_result, mode_frequency, &
       source_terms, iteration_config, allocate_solution, allocate_sources, zero_sources
  use sismo_poisson, only : solve_poisson_correction, prepare_poisson_cache, poisson_degree_cache
  use sismo_asymptotic, only : initial_asymptotic_g_frequency
  use sismo_io, only : open_result_file, write_result_row, write_comparison_file, write_eigenfunction_file
  implicit none
  private

  public :: run_core2

  ! assembled tridiagonal T(sigma^2): row i:  lo(i)*y_{i-1} + di(i)*y_i + up(i)*y_{i+1}
  type :: tri_matrix
     integer :: n = 0
     integer :: i0 = 1                        ! first active node (centre guard)
     real(dp), allocatable :: lo(:), di(:), up(:)
  end type tri_matrix

  ! Degree-dependent part of the mechanical operator.  Only a12 and a21
  ! depend on frequency; all remaining coefficients and integrating-factor
  ! geometry are prepared once and shared read-only by the OpenMP workers.
  type :: core2_degree_cache
     integer :: n = 0
     integer :: l = 0
     integer :: i0 = 1
     real(dp) :: bc_l1 = 0.0_dp
     real(dp), allocatable :: a11(:), a12_0(:), a12_1(:), a13(:)
     real(dp), allocatable :: a21_0(:), a21_1(:), a22(:), a23(:)
     real(dp), allocatable :: face_hphi(:), face_exp(:), face_a13(:)
     real(dp), allocatable :: node_vphi(:), node_exp(:)
  end type core2_degree_cache

  ! Frequency-dependent face data shared by the operator, source-vector, and
  ! y2 reconstruction at one value of sigma^2.
  type :: core2_operator_faces
     integer :: n = 0
     real(dp), allocatable :: a12(:), a21(:), a12_face(:)
     real(dp), allocatable :: den(:), cR(:), cL(:)
  end type core2_operator_faces

contains

  ! ---------------------------------------------------------------- utilities
  real(dp) function phi1(z) result(f)
    real(dp), intent(in) :: z
    if (abs(z) < 1.0d-6) then
       f = 1.0_dp + 0.5_dp*z + z*z/6.0_dp
    else
       f = (exp(z) - 1.0_dp)/z
    end if
  end function phi1

  real(dp) function cap_exp(z) result(e)
    real(dp), intent(in) :: z
    e = exp(max(-60.0_dp, min(60.0_dp, z)))
  end function cap_exp

  ! The syst05 coefficient matrices (identical to OSC SYSTEM 5; moved here
  ! from the legacy block-solver module for the release build).
  subroutine system05_coefficients(model, l, i, d0, d1)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l, i
    complex(dp), intent(out) :: d0(4,4), d1(4,4)

    real(dp) :: ell2, x, s2, s3, s4, s5, s6
    real(dp) :: c33, c21, c42, d13

    d0 = czero
    d1 = czero
    ell2 = real(l*(l+1), dp)
    x = max(model%x(i), tiny)
    s2 = max(tiny, model%qx3(i))
    s3 = max(tiny, model%pressure(i))/max(tiny, model%rho(i))
    s4 = model%rho(i)
    s5 = max(tiny, model%gamma1(i))
    s6 = model%aosc(i)

    c33 = (3.0_dp - s4/s2 - real(l, dp))/x
    c21 = x*s6
    c42 = x*s4/(s3*s5)
    d13 = x*s2/(s3*s5)

    d0(1,2) = cmplx(ell2*s2/x, 0.0_dp, kind=dp)
    d0(2,1) = cmplx(c21, 0.0_dp, kind=dp)
    d0(2,2) = cmplx(c33 - c21, 0.0_dp, kind=dp)
    d0(2,3) = cmplx(c21, 0.0_dp, kind=dp)
    d0(3,3) = cmplx(c33, 0.0_dp, kind=dp)
    d0(3,4) = cmplx(1.0_dp/x, 0.0_dp, kind=dp)
    d0(4,1) = -cmplx(c21*s4/s2, 0.0_dp, kind=dp)
    d0(4,2) = cmplx(c42, 0.0_dp, kind=dp)
    d0(4,3) = cmplx(ell2/x - c42, 0.0_dp, kind=dp)
    d0(4,4) = cmplx(c33 - 1.0_dp/x, 0.0_dp, kind=dp)

    d1(1,1) = cmplx(-real(l+1, dp)/x + d13, 0.0_dp, kind=dp)
    d1(1,2) = -cmplx(d13, 0.0_dp, kind=dp)
    d1(1,3) = cmplx(d13, 0.0_dp, kind=dp)
    d1(2,1) = cmplx(1.0_dp/(x*s2), 0.0_dp, kind=dp)
  end subroutine system05_coefficients

  subroutine prepare_core2_cache(model, l, cache)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    type(core2_degree_cache), intent(inout) :: cache

    complex(dp) :: d0(4,4), d1(4,4)
    real(dp) :: h, p, q, vol, s2, rho1
    integer :: i, n

    n = model%n
    if (cache%n /= n) then
       if (allocated(cache%a11)) then
          deallocate(cache%a11, cache%a12_0, cache%a12_1, cache%a13)
          deallocate(cache%a21_0, cache%a21_1, cache%a22, cache%a23)
          deallocate(cache%face_hphi, cache%face_exp, cache%face_a13)
          deallocate(cache%node_vphi, cache%node_exp)
       end if
       allocate(cache%a11(n), cache%a12_0(n), cache%a12_1(n), cache%a13(n))
       allocate(cache%a21_0(n), cache%a21_1(n), cache%a22(n), cache%a23(n))
       allocate(cache%face_hphi(n), cache%face_exp(n), cache%face_a13(n))
       allocate(cache%node_vphi(n), cache%node_exp(n))
       cache%n = n
    end if
    cache%l = l

    do i = 1, n
       call system05_coefficients(model, l, i, d0, d1)
       cache%a11(i) = real(d1(1,1), dp)
       cache%a12_0(i) = real(d0(1,2), dp)
       cache%a12_1(i) = real(d1(1,2), dp)
       cache%a13(i) = real(d1(1,3), dp)
       cache%a21_0(i) = real(d0(2,1), dp)
       cache%a21_1(i) = real(d1(2,1), dp)
       cache%a22(i) = real(d0(2,2), dp)
       cache%a23(i) = real(d0(2,3), dp)
    end do

    cache%i0 = 1
    do while (cache%i0 < n-2 .and. model%x(cache%i0) <= 1.0d-8)
       cache%i0 = cache%i0 + 1
    end do

    cache%face_hphi = 0.0_dp
    cache%face_exp = 0.0_dp
    cache%face_a13 = 0.0_dp
    do i = 1, n-1
       h = model%x(i+1) - model%x(i)
       p = 0.5_dp*(cache%a11(i) + cache%a11(i+1))*h
       cache%face_hphi(i) = h*phi1(max(-60.0_dp, min(60.0_dp, p)))
       cache%face_exp(i) = cap_exp(p)
       cache%face_a13(i) = 0.5_dp*(cache%a13(i) + cache%a13(i+1))
    end do

    cache%node_vphi = 0.0_dp
    cache%node_exp = 0.0_dp
    do i = cache%i0, n
       if (i == cache%i0) then
          vol = 0.5_dp*(model%x(i+1) - model%x(i))
       else if (i == n) then
          vol = 0.5_dp*(model%x(n) - model%x(n-1))
       else
          vol = 0.5_dp*(model%x(i+1) - model%x(i-1))
       end if
       q = cache%a22(i)*vol
       cache%node_vphi(i) = vol*phi1(q)
       cache%node_exp(i) = cap_exp(q)
    end do

    if (l == 1) then
       s2 = max(tiny, model%qx3(cache%i0))
       rho1 = model%rho(cache%i0)
       cache%bc_l1 = (rho1/3.0_dp - s2)/max(tiny, rho1/3.0_dp)
    else
       cache%bc_l1 = 0.0_dp
    end if
  end subroutine prepare_core2_cache

  ! The two frequency-dependent node coefficients of the mechanical pair.
  subroutine core2_frequency_coefficients(cache, w2, a12, a21)
    type(core2_degree_cache), intent(in) :: cache
    real(dp), intent(in) :: w2
    real(dp), intent(out) :: a12(:), a21(:)

    integer :: i

    do i = 1, cache%n
       ! Preserve the original arithmetic ordering exactly.
       a12(i) = cache%a12_0(i)/w2 + cache%a12_1(i)
       a21(i) = cache%a21_0(i) + w2*cache%a21_1(i)
    end do
  end subroutine core2_frequency_coefficients

  ! ------------------------------------------------------------ assembly
  ! Builds T(w2) and retains its frequency-dependent face data for reuse.
  ! Face f between nodes (i,i+1):  y2_f = (y1_{i+1} - e^{p_f} y1_i)/(h phi1(p_f) a12_f)
  !                                       - (a13_f/a12_f) y3_f
  ! Node i (dual cell V_i): y2_{f+} - e^{q_i} y2_{f-} - V phi1(q_i) a21_i y1_i
  !                          = V phi1(q_i) a23_i y3_i
  subroutine core2_assemble(cache, w2, T, faces, nnegden)
    type(core2_degree_cache), intent(in) :: cache
    real(dp), intent(in) :: w2
    type(tri_matrix), intent(inout) :: T
    type(core2_operator_faces), intent(inout) :: faces
    integer, intent(out), optional :: nnegden

    real(dp) :: den, eq
    integer :: i, n, i0

    n = cache%n
    if (T%n /= n) then
       if (allocated(T%lo)) deallocate(T%lo, T%di, T%up)
       allocate(T%lo(n), T%di(n), T%up(n))
       T%n = n
    end if
    T%lo = 0.0_dp; T%di = 0.0_dp; T%up = 0.0_dp
    if (present(nnegden)) nnegden = 0

    if (faces%n /= n) then
       if (allocated(faces%a12)) then
          deallocate(faces%a12, faces%a21, faces%a12_face)
          deallocate(faces%den, faces%cR, faces%cL)
       end if
       allocate(faces%a12(n), faces%a21(n), faces%a12_face(n))
       allocate(faces%den(n), faces%cR(n), faces%cL(n))
       faces%n = n
    end if
    call core2_frequency_coefficients(cache, w2, faces%a12, faces%a21)

    i0 = cache%i0
    T%i0 = i0

    ! face coefficients (face f between nodes i and i+1, stored at index i)
    faces%a12_face = 0.0_dp
    faces%den = 0.0_dp
    faces%cR = 0.0_dp
    faces%cL = 0.0_dp
    do i = 1, n-1
       faces%a12_face(i) = 0.5_dp*(faces%a12(i) + faces%a12(i+1))
       den = cache%face_hphi(i)*faces%a12_face(i)
       if (abs(den) < 1.0d-290) den = sign(1.0d-290, den + tiny)
       faces%den(i) = den
       faces%cR(i) = 1.0_dp/den
       faces%cL(i) = -cache%face_exp(i)/den
       if (present(nnegden) .and. i >= i0) then
          if (den < 0.0_dp) nnegden = nnegden + 1
       end if
    end do

    ! interior node equations i = i0..n-1 (node n carries the surface BC)
    do i = i0, n-1
       eq = cache%node_exp(i)
       ! y2_{f+} - e^q y2_{f-} - V phi1(q) a21 y1_i = V phi1(q) a23 y3_i
       ! f+ = face(i); f- = face(i-1) for interior, or the centre BC at i = i0.
       T%up(i) = faces%cR(i)
       T%di(i) = faces%cL(i) - cache%node_vphi(i)*faces%a21(i)
       if (i == i0) then
          if (cache%l == 1) then
             T%di(i) = T%di(i) - eq*cache%bc_l1
          end if
       else
          T%di(i) = T%di(i) - eq*faces%cR(i-1)
          T%lo(i) = -eq*faces%cL(i-1)
       end if
    end do

    ! surface row at node n:  y1 - y2 + y3 = 0  with y2(x_n) from the last face
    ! half-cell:  y2(x_n) = e^q y2_{f-} + V phi1(q) (a21 y1_n + a23 y3_n)
    eq = cache%node_exp(n)
    ! y1_n + y3_n = e^q [cR(n-1) y1_n + cL(n-1) y1_{n-1} + sF(n-1)]
    !               + V phi1(q) (a21 y1_n + a23 y3_n)
    T%di(n) = 1.0_dp - eq*faces%cR(n-1) - cache%node_vphi(n)*faces%a21(n)
    T%lo(n) = -eq*faces%cL(n-1)
    T%up(n) = 0.0_dp

    ! Dirichlet rows: inactive nodes below i0, and the l>=2 centre rule
    do i = 1, i0-1
       T%di(i) = 1.0_dp; T%lo(i) = 0.0_dp; T%up(i) = 0.0_dp
    end do
    if (cache%l >= 2) then
       T%di(i0) = 1.0_dp; T%lo(i0) = 0.0_dp; T%up(i0) = 0.0_dp
    end if
  end subroutine core2_assemble

  ! Builds only the source vector for an already assembled T(w2).
  subroutine core2_build_rhs(cache, faces, y3, rhs)
    type(core2_degree_cache), intent(in) :: cache
    type(core2_operator_faces), intent(in) :: faces
    real(dp), intent(in) :: y3(:)
    real(dp), intent(out) :: rhs(:)

    real(dp), allocatable :: sF(:)
    real(dp) :: eq
    integer :: i, n, i0

    n = cache%n
    i0 = cache%i0
    allocate(sF(n))
    sF = 0.0_dp
    rhs(1:n) = 0.0_dp
    do i = i0, n-1
       ! Use the same guarded denominator as the assembled operator.  Written
       ! this way the face_hphi factor cancels analytically away from an exact
       ! Lamb crossing, while the guard remains effective at the crossing.
       sF(i) = -cache%face_a13(i)*cache%face_hphi(i)/faces%den(i) &
            *0.5_dp*(y3(i) + y3(i+1))
    end do

    do i = i0, n-1
       eq = cache%node_exp(i)
       if (i == i0) then
          if (cache%l == 1) then
             rhs(i) = cache%node_vphi(i)*cache%a23(i)*y3(i) - sF(i) + &
                  eq*cache%bc_l1*y3(i)
          else
             rhs(i) = cache%node_vphi(i)*cache%a23(i)*y3(i) - sF(i)
          end if
       else
          rhs(i) = cache%node_vphi(i)*cache%a23(i)*y3(i) - sF(i) + eq*sF(i-1)
       end if
    end do

    eq = cache%node_exp(n)
    rhs(n) = -y3(n) + eq*sF(n-1) + cache%node_vphi(n)*cache%a23(n)*y3(n)
    rhs(1:i0-1) = 0.0_dp
    if (cache%l >= 2) rhs(i0) = 0.0_dp
  end subroutine core2_build_rhs

  ! ----------------------------------------------------- LU / Sturm machinery
  ! LU without pivoting; pivots d_i = di_i - lo_i*up_{i-1}/d_{i-1}.
  ! For a tridiagonal with positive off-diagonal products (ours, by construction)
  ! the negative-pivot count is a rigorous Sturm count.
  subroutine tri_pivots(T, piv, nneg, info)
    type(tri_matrix), intent(in) :: T
    real(dp), intent(out) :: piv(:)
    integer, intent(out) :: nneg, info

    integer :: i
    info = 0
    nneg = 0
    piv(1) = T%di(1)
    if (piv(1) < 0.0_dp) nneg = nneg + 1
    do i = 2, T%n
       if (abs(piv(i-1)) < 1.0d-300) then
          piv(i-1) = sign(1.0d-300, piv(i-1) + tiny)
       end if
       piv(i) = T%di(i) - T%lo(i)*T%up(i-1)/piv(i-1)
       if (piv(i) < 0.0_dp) nneg = nneg + 1
    end do
  end subroutine tri_pivots

  subroutine tri_solve(T, piv, b, x, work)
    type(tri_matrix), intent(in) :: T
    real(dp), intent(in) :: piv(:), b(:)
    real(dp), intent(out) :: x(:)
    real(dp), intent(out) :: work(:)

    integer :: i
    work(1) = b(1)
    do i = 2, T%n
       work(i) = b(i) - T%lo(i)*work(i-1)/piv(i-1)
    end do
    x(T%n) = work(T%n)/piv(T%n)
    do i = T%n-1, 1, -1
       x(i) = (work(i) - T%up(i)*x(i+1))/piv(i)
    end do
  end subroutine tri_solve

  ! Keep identity rows homogeneous in every inverse solve.  In particular, a
  ! random inverse-iteration seed must not leak into inactive centre nodes.
  subroutine zero_dirichlet_rhs(cache, b)
    type(core2_degree_cache), intent(in) :: cache
    real(dp), intent(inout) :: b(:)

    if (cache%i0 > 1) b(1:cache%i0-1) = 0.0_dp
    if (cache%l >= 2) b(cache%i0) = 0.0_dp
  end subroutine zero_dirichlet_rhs

  ! Allocation-free Sturm count using workspace owned by one OpenMP thread.
  integer function sturm_count_with_workspace(cache, w2, Tl, facesl, pivl) result(indicator)
    type(core2_degree_cache), intent(in) :: cache
    real(dp), intent(in) :: w2
    type(tri_matrix), intent(inout) :: Tl
    type(core2_operator_faces), intent(inout) :: facesl
    real(dp), intent(out) :: pivl(:)
    integer :: info, nneg, nnegden
    call core2_assemble(cache, w2, Tl, facesl, nnegden=nnegden)
    call tri_pivots(Tl, pivl, nneg, info)
    indicator = nneg + nnegden
  end function sturm_count_with_workspace



  ! Sweep [om_lo, om_hi] uniformly in period, bisect every parity change.
  ! PARALLEL: the sweep points are independent factorizations, and so is each
  ! crossing's bisection -- both loops are OpenMP-parallel with one reusable
  ! matrix/factorization workspace per thread.
  subroutine core2_sturm_scan(cache, om_lo, om_hi, nscan, eigs, neig)
    type(core2_degree_cache), intent(in) :: cache
    integer, intent(in) :: nscan
    real(dp), intent(in) :: om_lo, om_hi
    real(dp), allocatable, intent(out) :: eigs(:)
    integer, intent(out) :: neig

    real(dp), allocatable :: tgrid(:), found(:)
    real(dp), allocatable :: scan_piv(:)
    integer, allocatable :: counts(:), cross(:)
    type(tri_matrix) :: scan_T
    type(core2_operator_faces) :: scan_faces
    real(dp) :: t_lo, t_hi, tL, tR, tm, om
    integer :: k, ncross, j, cL, cm, itb

    t_lo = 1.0_dp/max(tiny, om_hi)
    t_hi = 1.0_dp/max(tiny, om_lo)
    allocate(tgrid(nscan), counts(nscan))
    do k = 1, nscan
       tgrid(k) = t_lo + (t_hi - t_lo)*real(k-1, dp)/real(nscan-1, dp)
    end do
    !$omp parallel private(scan_T, scan_faces, scan_piv)
    scan_T%n = 0
    scan_faces%n = 0
    allocate(scan_piv(cache%n))
    !$omp do schedule(static)
    do k = 1, nscan
       counts(k) = sturm_count_with_workspace(cache, (1.0_dp/tgrid(k))**2, &
            scan_T, scan_faces, scan_piv)
    end do
    !$omp end do
    deallocate(scan_piv)
    !$omp end parallel

    ! Collect intervals containing one crossing.  A larger count jump means
    ! the requested initial scan is too coarse and is rejected below.
    allocate(cross(nscan-1), found(nscan-1))
    ncross = 0
    do k = 2, nscan
       if (abs(counts(k) - counts(k-1)) > 1) then
          write(*,'(A,I0,A,I0,A)') ' SISMO: scan_points=', nscan, &
               ' is too small near scan sample ', k, '; increase it in sismo.conf.'
          error stop 2
       end if
       if (abs(counts(k) - counts(k-1)) == 1) then
          ncross = ncross + 1
          cross(ncross) = k
       end if
    end do

    ! bisect each crossing in parallel
    !$omp parallel private(scan_T, scan_faces, scan_piv, tL, tR, tm, om, cL, cm, itb)
    scan_T%n = 0
    scan_faces%n = 0
    allocate(scan_piv(cache%n))
    !$omp do schedule(dynamic)
    do j = 1, ncross
       tL = tgrid(cross(j)-1)
       tR = tgrid(cross(j))
       cL = counts(cross(j)-1)
       do itb = 1, 60
          tm = 0.5_dp*(tL + tR)
          om = 1.0_dp/tm
          cm = sturm_count_with_workspace(cache, om*om, scan_T, scan_faces, scan_piv)
          if (mod(abs(cm - cL), 2) == 1) then
             tR = tm
          else
             tL = tm
          end if
          if ((tR - tL) < 1.0d-12*tR) exit
       end do
       found(j) = 2.0_dp/(tL + tR)
    end do
    !$omp end do
    deallocate(scan_piv)
    !$omp end parallel

    neig = ncross
    allocate(eigs(max(1, neig)))
    do j = 1, neig
       eigs(j) = found(j)
    end do
    call sort_desc(eigs, neig)
  end subroutine core2_sturm_scan

  subroutine sort_desc(a, n)
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: n
    integer :: i, j
    real(dp) :: t
    do i = 2, n
       t = a(i); j = i - 1
       do while (j >= 1)
          if (a(j) >= t) exit
          a(j+1) = a(j); j = j - 1
       end do
       a(j+1) = t
    end do
  end subroutine sort_desc

  ! Cowling eigenvector by inverse iteration at a slightly shifted frequency.
  subroutine core2_eigenvector(cache, w2, y1)
    type(core2_degree_cache), intent(in) :: cache
    real(dp), intent(in) :: w2
    real(dp), allocatable, intent(out) :: y1(:)

    type(tri_matrix) :: T
    type(core2_operator_faces) :: faces
    real(dp), allocatable :: piv(:), b(:), solve_work(:)
    real(dp) :: nrm
    integer :: i, it, nneg, info

    allocate(y1(cache%n), piv(cache%n), b(cache%n), solve_work(cache%n))
    call core2_assemble(cache, w2*(1.0_dp + 1.0d-9), T, faces)
    call tri_pivots(T, piv, nneg, info)
    do i = 1, cache%n
       y1(i) = 1.0_dp - 2.0_dp*mod(i*2654435761_8/65536_8, 1000_8)/1000.0_dp
    end do
    do it = 1, 3
       b = y1
       call zero_dirichlet_rhs(cache, b)
       call tri_solve(T, piv, b, y1, solve_work)
       nrm = maxval(abs(y1))
       if (nrm > tiny) y1 = y1/nrm
    end do
  end subroutine core2_eigenvector

  ! y2 at the nodes from the converged y1 (+ optional source), via face relations
  ! averaged onto nodes (used to build xi_h for the Poisson source).
  subroutine core2_y2_nodes(cache, w2, y1, y3, y2, faces)
    type(core2_degree_cache), intent(in) :: cache
    real(dp), intent(in) :: w2
    real(dp), intent(in) :: y1(:)
    real(dp), intent(in), optional :: y3(:)
    real(dp), intent(out) :: y2(:)
    type(core2_operator_faces), intent(in), optional :: faces

    real(dp), allocatable :: a12(:), y2f(:)
    real(dp) :: den, src
    integer :: i, n

    n = cache%n
    allocate(y2f(n))
    if (.not. present(faces)) then
       allocate(a12(n))
       do i = 1, n
          a12(i) = cache%a12_0(i)/w2 + cache%a12_1(i)
       end do
    end if
    y2f = 0.0_dp
    do i = 1, n-1
       if (present(faces)) then
          den = faces%den(i)
       else
          den = cache%face_hphi(i)*0.5_dp*(a12(i) + a12(i+1))
          if (abs(den) < 1.0d-290) den = sign(1.0d-290, den + tiny)
       end if
       src = 0.0_dp
       if (present(y3)) then
          ! As in core2_build_rhs, retain face_hphi explicitly so that den is
          ! the sole (regularized) division at a Lamb crossing.
          src = -cache%face_a13(i)*cache%face_hphi(i)/den &
               *0.5_dp*(y3(i) + y3(i+1))
       end if
       y2f(i) = (y1(i+1) - cache%face_exp(i)*y1(i))/den + src
    end do
    ! Lamb-crossing guard: at the face where a12 -> 0 the reconstruction divides
    ! rounding noise by ~zero; the physical y2 is smooth there, so replace clear
    ! outliers by the neighbour average (they would otherwise contaminate the
    ! Poisson source and destabilise the phi iteration).
    do i = 2, n-2
       if (abs(y2f(i)) > 10.0_dp*(abs(y2f(i-1)) + abs(y2f(i+1)))) then
          y2f(i) = 0.5_dp*(y2f(i-1) + y2f(i+1))
       end if
    end do
    y2(1) = y2f(1)
    do i = 2, n-1
       y2(i) = 0.5_dp*(y2f(i-1) + y2f(i))
    end do
    y2(n) = y2f(n-1)
  end subroutine core2_y2_nodes

  ! Fill a mode_solution (the Poisson solver's interface) from (y1, y2):
  ! sol%xi_r = x^{l-1} y1 ; sol%xi_h = y2 x^{l-1} qx3 / sigma^2  (pack convention)
  subroutine core2_fill_solution(model, l, w2, y1, y2, sol, y3)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    real(dp), intent(in) :: w2, y1(:), y2(:)
    type(mode_solution), intent(inout) :: sol
    real(dp), intent(in), optional :: y3(:)

    real(dp) :: sd, sp, s2, s3, y3i, dx
    integer :: i, ilo, ihi

    if (sol%n /= model%n) call allocate_solution(sol, model%n)
    sol%l = l
    sol%omega = cmplx(sqrt(max(tiny, w2)), 0.0_dp, kind=dp)
    do i = 1, model%n
       sd = max(model%x(i), 0.0_dp)**max(0, l-1)
       if (max(0, l-1) == 0) sd = 1.0_dp
       sp = max(model%x(i), 0.0_dp)**max(0, l)
       if (max(0, l) == 0) sp = 1.0_dp
       s2 = max(tiny, model%qx3(i))
       s3 = max(tiny, model%pressure(i))/max(tiny, model%rho(i))
       y3i = 0.0_dp
       if (present(y3)) y3i = y3(i)
       sol%xi_r(i) = cmplx(sd*y1(i), 0.0_dp, kind=dp)
       sol%xi_h(i) = cmplx(y2(i)*sd*s2/max(tiny, w2), 0.0_dp, kind=dp)
       ! Horizontal momentum gives P'/P = x^l*S2/S3*(y2-y3).
       sol%dpp(i) = cmplx(sp*s2/s3*(y2(i) - y3i), 0.0_dp, kind=dp)
       if (present(y3)) then
          ! y3 is the regularized potential phi/(x^l q/x^3); mode_solution
          ! stores the physical potential used by the eigenfunction writer.
          sol%phi(i) = cmplx(sp*s2*y3(i), 0.0_dp, kind=dp)
       else
          sol%phi(i) = cmplx(0.0_dp, 0.0_dp, kind=dp)
       end if
       sol%theta(i) = cmplx(0.0_dp, 0.0_dp, kind=dp)
    end do
    if (present(y3) .and. model%n >= 2) then
       ! Store the physical d(Phi')/dx corresponding to the physical potential
       ! above.  The split iteration carries only y3, so recover its derivative
       ! on the same full grid for eigenfunction output.
       do i = 1, model%n
          ilo = max(1, i-1)
          ihi = min(model%n, i+1)
          dx = model%x(ihi) - model%x(ilo)
          if (dx > tiny) sol%theta(i) = (sol%phi(ihi) - sol%phi(ilo)) &
               /cmplx(dx, 0.0_dp, kind=dp)
       end do
    end if
  end subroutine core2_fill_solution

  ! ------------------------------------------------- split (inhomogeneous) solve
  ! Refine one mode with the frozen-phi source: iterate
  !   phi <- Poisson(mechanics) ; solve T(w2) w = rhs(phi) ;
  !   dw2 = <y,z>/<y,t>, z = T^{-1}(rhs - T y), t = T^{-1} T' y
  subroutine core2_split_refine(model, l, mechanical_cache, w2_seed, config, poisson_cache, &
       w2_out, y1, y3_out, niter, conv)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    type(core2_degree_cache), intent(in) :: mechanical_cache
    real(dp), intent(in) :: w2_seed
    type(iteration_config), intent(in) :: config
    type(poisson_degree_cache), intent(in) :: poisson_cache
    real(dp), intent(out) :: w2_out
    real(dp), allocatable, intent(inout) :: y1(:)
    real(dp), allocatable, intent(out) :: y3_out(:)
    integer, intent(out) :: niter
    logical, intent(out) :: conv

    ! TIGHT-COUPLING port of the SISMO-1 split iteration:
    !   outer loop: one clamped frequency (Newton) step per pass;
    !   inner loop: at FIXED w2, co-iterate (shape, phi) to self-consistency --
    !     phi <- relax*Poisson(mechanics) + (1-relax)*phi_old,
    !     y1  <- normalized forced response T^{-1} rhs(phi)
    !     stop on the rate-corrected error estimate change*q/(1-q) < ptol
    !     (capped at n_inner refreshes; T factored ONCE per outer pass).
    ! This is what stops the dipole loose-coupling drift: the lambda step is
    ! taken only with phi converged, so the update's zero is the split fixed
    ! point, not the lagged-phi drift direction.
    type(tri_matrix) :: T, Tp
    type(core2_operator_faces) :: faces, facesp
    type(mode_solution) :: sol
    type(source_terms) :: src
    real(dp), allocatable :: piv(:), rhs(:), z(:), tvec(:), y2(:), y3(:), y3new(:), b(:), solve_work(:)
    real(dp) :: w2, num, den, dw2, dw2_step, nrm, relch, eps_fd
    real(dp) :: pchange, pchange_prev, qrate, perr
    real(dp) :: relax, ptol, clamp
    integer, parameter :: l1_stagnation_min_iter = 30
    integer, parameter :: l1_stable_best_required = 20
    integer, parameter :: l1_unclamped_required = 10
    integer, parameter :: l1_alternating_clamps_required = 10
    real(dp), parameter :: l1_small_correction_limit = 1.0d-2
    real(dp) :: w2_best, best_res, w2_candidate, best_w2_marker, previous_dw2
    real(dp), allocatable :: y1_best(:), y3_best(:)
    integer :: it, inner, n_inner, i, n, nneg, info
    integer :: stable_best_passes, unclamped_passes, clamp_alternations
    logical :: inner_ok, best_valid, clamp_active, previous_clamped

    n = model%n
    allocate(piv(n), rhs(n), z(n), tvec(n), y2(n), y3(n), y3new(n), b(n), solve_work(n))
    call allocate_sources(src, n)
    w2 = w2_seed
    conv = .false.
    niter = 0
    y3 = 0.0_dp
    allocate(y1_best(n), y3_best(n))
    w2_best = w2_seed
    y1_best = y1
    y3_best = y3
    best_res = huge(1.0_dp)
    best_valid = .false.
    best_w2_marker = w2_seed
    stable_best_passes = 0
    unclamped_passes = 0
    clamp_alternations = 0
    previous_dw2 = 0.0_dp
    previous_clamped = .false.
    relax = max(0.05_dp, min(1.0_dp, config%poisson_relax))
    n_inner = max(1, config%inner_poisson_iters)
    ptol = config%poisson_inner_tol            ! <=0: fixed n_inner refreshes
    clamp = 1.0d-3

    do it = 1, max(10, config%max_iter)
       niter = it
       ! ---- factor T(w2) once per outer pass
       call core2_assemble(mechanical_cache, w2, T, faces)
       call tri_pivots(T, piv, nneg, info)
       ! ---- inner loop: co-converge (y1, phi) at fixed w2 -- INCREMENTAL shape
       ! update y1 <- y1 + z (SISMO-1 semantics; keeps sign/shape continuity),
       ! z = T^{-1}(rhs - T y1) the residual correction toward the forced response.
       pchange_prev = huge(1.0_dp)
       inner_ok = .false.
       do inner = 1, n_inner
          call core2_y2_nodes(mechanical_cache, w2, y1, y3=y3, y2=y2, faces=faces)
          call core2_fill_solution(model, l, w2, y1, y2, sol)
          call zero_sources(src)
          call solve_poisson_correction(model, l, sol, src, poisson_cache)
          do i = 1, n
             y3new(i) = real(sol%phi(i), dp)/(max(tiny, max(model%x(i), 0.0_dp)**max(0,l)) &
                  *max(tiny, model%qx3(i)))
          end do
          y3new = relax*y3new + (1.0_dp - relax)*y3
          nrm = maxval(abs(y3new))
          pchange = 0.0_dp
          if (nrm > tiny) pchange = maxval(abs(y3new - y3))/nrm
          y3 = y3new
          call core2_build_rhs(mechanical_cache, faces, y3, rhs)
          b = rhs
          b(1) = b(1) - (T%di(1)*y1(1) + T%up(1)*y1(2))
          do i = 2, n-1
             b(i) = b(i) - (T%lo(i)*y1(i-1) + T%di(i)*y1(i) + T%up(i)*y1(i+1))
          end do
          b(n) = b(n) - (T%lo(n)*y1(n-1) + T%di(n)*y1(n))
          call zero_dirichlet_rhs(mechanical_cache, b)
          call tri_solve(T, piv, b, z, solve_work)
          ! sign-align: the eigenvector sign is arbitrary; if the update lands
          ! with negative overlap the normalized iterate flips sign every pass
          ! and phi never converges (pchange ~ 2).  Align with the previous y1.
          if (dot_product(y1, y1 + z) < 0.0_dp) then
             y1 = -(y1 + z)
          else
             y1 = y1 + z
          end if
          nrm = maxval(abs(y1))
          if (nrm > tiny) y1 = y1/nrm
          if (ptol > 0.0_dp .and. inner >= 2 .and. pchange_prev > tiny) then
             qrate = pchange/pchange_prev
             if (qrate < 0.999_dp) then
                perr = pchange*qrate/(1.0_dp - qrate)
                if (perr <= ptol) then
                   inner_ok = .true.
                   exit
                end if
             end if
          end if
          pchange_prev = pchange
       end do
       if (.not. inner_ok) inner_ok = (pchange <= max(1.0d-4, 10.0_dp*max(0.0_dp, ptol)))
       ! ---- one frequency step with converged phi
       ! z = T^{-1}(rhs - T y1)
       b = rhs
       b(1) = b(1) - (T%di(1)*y1(1) + T%up(1)*y1(2))
       do i = 2, n-1
          b(i) = b(i) - (T%lo(i)*y1(i-1) + T%di(i)*y1(i) + T%up(i)*y1(i+1))
       end do
       b(n) = b(n) - (T%lo(n)*y1(n-1) + T%di(n)*y1(n))
       call zero_dirichlet_rhs(mechanical_cache, b)
       call tri_solve(T, piv, b, z, solve_work)
       ! t = T^{-1} T' y1  (finite-difference T')
       ! Scale the perturbation to max(1,|w2|), retaining the established
       ! low-frequency step while adding a multiple-ULP floor so that
       ! w2 + eps_fd is always numerically distinct from w2.
       eps_fd = max(1.0d-7*max(1.0_dp, abs(w2)), 32.0_dp*spacing(w2))
       call core2_assemble(mechanical_cache, w2 + eps_fd, Tp, facesp)
       b(1) = (Tp%di(1)-T%di(1))*y1(1) + (Tp%up(1)-T%up(1))*y1(2)
       do i = 2, n-1
          b(i) = (Tp%lo(i)-T%lo(i))*y1(i-1) + (Tp%di(i)-T%di(i))*y1(i) &
               + (Tp%up(i)-T%up(i))*y1(i+1)
       end do
       b(n) = (Tp%lo(n)-T%lo(n))*y1(n-1) + (Tp%di(n)-T%di(n))*y1(n)
       b = b/eps_fd
       call zero_dirichlet_rhs(mechanical_cache, b)
       call tri_solve(T, piv, b, tvec, solve_work)
       num = dot_product(y1, z)
       den = dot_product(y1, tvec)
       if (abs(den) <= tiny) exit
       dw2 = num/den
       if (ieee_is_nan(dw2)) exit
       relch = abs(dw2)/w2
       clamp_active = relch > clamp
       if (clamp_active) then
          unclamped_passes = 0
          if (previous_clamped .and. ((dw2 < 0.0_dp .and. previous_dw2 > 0.0_dp) .or. &
               (dw2 > 0.0_dp .and. previous_dw2 < 0.0_dp))) then
             clamp_alternations = clamp_alternations + 1
          else
             clamp_alternations = 1
          end if
       else
          unclamped_passes = unclamped_passes + 1
          clamp_alternations = 0
       end if
       previous_dw2 = dw2
       previous_clamped = clamp_active
       if (best_valid) stable_best_passes = stable_best_passes + 1
       dw2_step = dw2
       if (abs(dw2_step) > clamp*w2) dw2_step = sign(clamp*w2, dw2_step)
       b = y1 + z - dw2_step*tvec
       if (dot_product(y1, b) < 0.0_dp) b = -b
       nrm = maxval(abs(b))
       if (nrm > tiny) b = b/nrm
       w2_candidate = max(tiny, w2 + dw2_step)
       ! remember the best-visited iterate: the (phi,shape,lambda) loop can enter
       ! a small limit cycle instead of converging; the point of smallest raw
       ! Newton step is the closest pass to the resonance.  Snapshot the whole
       ! updated split state so frequency, mechanics, and frozen potential all
       ! come from the same algorithmic iterate.
       if (it >= 3 .and. pchange < 1.0d-1 .and. relch < best_res) then
          if (.not. best_valid .or. &
               abs(w2_candidate - best_w2_marker)/max(w2_seed, tiny) > &
               max(tiny, 0.1_dp*config%tol)) then
             best_w2_marker = w2_candidate
             stable_best_passes = 0
          end if
          best_valid = .true.
          best_res = relch
          w2_best = w2_candidate
          y1_best = b
          if (config%write_eigenfunctions) then
             ! The Newton candidate (b,w2_candidate) is one step newer than
             ! y3.  Refresh its potential before snapshotting so an early
             ! stagnation exit cannot return a mixed-iterate eigenfunction.
             call core2_y2_nodes(mechanical_cache, w2_candidate, b, y3=y3, y2=y2)
             call core2_fill_solution(model, l, w2_candidate, b, y2, sol)
             call zero_sources(src)
             call solve_poisson_correction(model, l, sol, src, poisson_cache)
             do i = 1, n
                y3_best(i) = real(sol%phi(i), dp) &
                     /(max(tiny, max(model%x(i), 0.0_dp)**max(0,l)) &
                     *max(tiny, model%qx3(i)))
             end do
          else
             y3_best = y3
          end if
       end if
       y1 = b
       w2 = w2_candidate
       ! converge on the RAW (unclamped) Newton step vanishing, with phi converged
       if (it >= 3 .and. relch < max(1.0d-7, config%tol) .and. inner_ok) then
          conv = .true.
          exit
       end if
       ! Dipole iterations often enter a bounded fixed/period-2 cycle after
       ! reaching their best frequency.  Stop only when that best frequency has
       ! been stable for many passes, the non-Cowling correction is small, and
       ! either raw Newton steps are no longer clamped or a persistent alternating
       ! clamped cycle has been observed.  Large low-order corrections retain the
       ! full iteration budget.  This is a stagnation exit: conv remains false.
       if (l == 1 .and. it >= l1_stagnation_min_iter .and. best_valid .and. &
            stable_best_passes >= l1_stable_best_required .and. &
            abs(w2_best - w2_seed)/max(w2_seed, tiny) <= l1_small_correction_limit .and. &
            (unclamped_passes >= l1_unclamped_required .or. &
             clamp_alternations >= l1_alternating_clamps_required)) exit
    end do
    if (conv) then
       w2_out = w2
       allocate(y3_out(n))
       y3_out = y3
    else
       ! unconverged (limit cycle / cap): return the best pass, not the last
       w2_out = w2_best
       y1 = y1_best
       allocate(y3_out(n))
       y3_out = y3_best
    end if
    if (config%write_eigenfunctions) then
       ! One final Poisson refresh makes the optional exported potential belong
       ! to the returned mechanical shape even when refinement stopped on the
       ! bounded best iterate rather than on the last outer pass.
       call core2_y2_nodes(mechanical_cache, w2_out, y1, y3=y3_out, y2=y2)
       call core2_fill_solution(model, l, w2_out, y1, y2, sol)
       call zero_sources(src)
       call solve_poisson_correction(model, l, sol, src, poisson_cache)
       do i = 1, n
          y3_out(i) = real(sol%phi(i), dp) &
               /(max(tiny, max(model%x(i), 0.0_dp)**max(0,l)) &
               *max(tiny, model%qx3(i)))
       end do
    end if
  end subroutine core2_split_refine

  ! ------------------------------------------------------------------ driver
  subroutine run_core2(model, config, l_min, l_max, g_min, g_max, output_base)
    type(stellar_model), intent(in) :: model
    type(iteration_config), intent(in) :: config
    integer, intent(in) :: l_min, l_max, g_min, g_max
    character(len=*), intent(in) :: output_base

    type(mode_result), allocatable :: results(:)
    type(poisson_degree_cache) :: poisson_cache
    type(core2_degree_cache) :: mechanical_cache
    real(dp), allocatable :: eigs(:), y1(:), y2(:), y3(:)
    real(dp) :: om_lo, om_hi, w2, w2_out, om_asym_top, om_asym_bot
    integer :: l, k, neig, nres, unit, nsel, niter, i, nstart
    logical :: conv

    allocate(results((l_max-l_min+1)*(g_max-g_min+1)))
    nres = 0

    do l = max(1, l_min), l_max
       call prepare_core2_cache(model, l, mechanical_cache)
       ! frequency window from the asymptotic relation (range only, never seeds)
       ! Always start the ordered scan at global radial order one.  Starting
       ! the window at g_min makes the first eigenvalue found look like n=-1
       ! and therefore mislabels every result when g_min > 1.
       om_asym_top = initial_asymptotic_g_frequency(model, l, 1, config%eps_g)
       om_asym_bot = initial_asymptotic_g_frequency(model, l, g_max + 3, config%eps_g)
       om_hi = 3.0_dp*om_asym_top
       om_lo = 0.75_dp*om_asym_bot
       call core2_sturm_scan(mechanical_cache, om_lo, om_hi, config%sturm_scan_points, eigs, neig)
       write(*,'(A,I3,A,I5,A,1P,2E12.4)') ' core2: l=', l, '  Sturm eigenvalues found=', neig, &
            '  window=', om_lo, om_hi
       nsel = 0
       nstart = nres
       do k = 1, neig
          nsel = nsel + 1
          if (nsel < g_min .or. nsel > g_max) cycle
          nres = nres + 1
          results(nres)%mode%l = l
          results(nres)%mode%m = 0
          results(nres)%mode%k = 1
          results(nres)%mode%n = -nsel
          results(nres)%mode%ltilde = real(l, dp)
          results(nres)%mode%omega_ad = eigs(k)
          results(nres)%mode%freq_com_cd = eigs(k)*model%fdy*0.0864_dp
          results(nres)%mode%freq_in_cd = results(nres)%mode%freq_com_cd
          results(nres)%num = nsel
          results(nres)%err_takata = 0.0_dp
          results(nres)%final_stride = 1
          results(nres)%inertia = 0.0_dp
          results(nres)%work = 0.0_dp
       end do
       if (config%use_poisson) call prepare_poisson_cache(model, l, poisson_cache)
       ! refine the selected modes in parallel (independent eigenproblems)
       !$omp parallel do private(w2, w2_out, y1, y2, y3, niter, conv) schedule(dynamic)
       do k = nstart+1, nres
          w2 = results(k)%mode%omega_ad**2
          call core2_eigenvector(mechanical_cache, w2, y1)
          if (config%use_poisson) then
             call core2_split_refine(model, l, mechanical_cache, w2, config, poisson_cache, &
                  w2_out, y1, y3, niter, conv)
             results(k)%iterations = niter
             results(k)%converged = conv
          else
             w2_out = w2
             allocate(y3(model%n))
             y3 = 0.0_dp
             results(k)%iterations = 0
             results(k)%converged = .true.
          end if
          if (config%write_eigenfunctions) then
             call allocate_solution(results(k)%solution, model%n)
             allocate(y2(model%n))
             call core2_y2_nodes(mechanical_cache, w2_out, y1, y3=y3, y2=y2)
             call core2_fill_solution(model, l, w2_out, y1, y2, results(k)%solution, y3=y3)
             deallocate(y2)
          else
             results(k)%solution%n = model%n
             results(k)%solution%l = l
             results(k)%solution%omega = cmplx(sqrt(max(tiny, w2_out)), 0.0_dp, kind=dp)
          end if
          if (allocated(y1)) deallocate(y1)
          if (allocated(y3)) deallocate(y3)
       end do
       !$omp end parallel do
    end do

    call open_result_file(output_base, unit)
    do i = 1, nres
       call write_result_row(unit, results(i))
    end do
    close(unit)
    call write_comparison_file(output_base, results, nres)
    if (config%write_eigenfunctions) then
       do i = 1, nres
          call write_eigenfunction_file(model, output_base, results(i))
       end do
    end if
    write(*,*) 'SISMO core2: wrote ', trim(output_base)//'.sismo'
  end subroutine run_core2

end module sismo_core2
