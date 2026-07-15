module sismo_poisson
  use sismo_precision, only : dp, tiny, czero
  use sismo_types, only : stellar_model, mode_solution, source_terms, &
       allocate_sources, zero_sources
  implicit none
  private

  public :: solve_poisson_correction, set_takata_closure, prepare_poisson_cache

  type, public :: poisson_degree_cache
     integer :: n = 0
     integer :: l = -1
     logical :: prepared = .false.
     complex(dp), allocatable :: homogeneous(:,:)
     ! (x, s2, s4, c33, c21, c42) at left/mid/right for each cell.  Only the
     ! frequency-independent raw RHS values are cached; the RK4 arithmetic and
     ! expression grouping remain identical to the uncached path.
     real(dp), allocatable :: rk4_coeff(:,:,:)
     complex(dp) :: bc_homogeneous = czero
  end type poisson_degree_cache

  ! Release configuration: the split Poisson is solved by RK4 shooting with the
  ! vacuum surface match (with an automatic Green-function integral fallback when
  ! the homogeneous match degenerates).  For l=1 the amplitude can instead be
  ! closed by Takata's first integral (momentum conservation) -- the release
  ! default for dipole runs (takata_closure = true in sismo.conf).
  logical, save :: use_takata_closure = .false.

contains

  subroutine set_takata_closure(flag)
    logical, intent(in) :: flag
    use_takata_closure = flag
  end subroutine set_takata_closure

  subroutine prepare_poisson_cache(model, l, cache)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    type(poisson_degree_cache), intent(inout) :: cache

    complex(dp), allocatable :: zero_mechanics(:,:)
    integer :: i, n, start_i

    n = model%n
    cache%prepared = .false.
    cache%n = n
    cache%l = l
    cache%bc_homogeneous = czero
    if (allocated(cache%homogeneous)) deallocate(cache%homogeneous)
    if (allocated(cache%rk4_coeff)) deallocate(cache%rk4_coeff)
    allocate(cache%homogeneous(2,n), zero_mechanics(2,n))
    cache%homogeneous = czero
    zero_mechanics = czero
    if (n <= 1) then
       allocate(cache%rk4_coeff(6,3,0))
       cache%prepared = .true.
       return
    end if

    cache%homogeneous(1,1) = cone()
    cache%homogeneous(2,1) = cmplx(real(max(l, 0), dp), 0.0_dp, kind=dp)
    start_i = 1
    if (model%x(1) <= tiny .and. n >= 2) then
       cache%homogeneous(:,2) = cache%homogeneous(:,1)
       start_i = 2
    end if
    do i = start_i, n-1
       call rk4_poisson_step(model, l, zero_mechanics, i, cache%homogeneous(:,i), .false., &
            cache%homogeneous(:,i+1))
    end do
    cache%bc_homogeneous = cmplx(real(l+1, dp), 0.0_dp, kind=dp)*cache%homogeneous(1,n) &
         + cache%homogeneous(2,n)
    allocate(cache%rk4_coeff(6,3,n-1))
    cache%rk4_coeff = 0.0_dp
    do i = start_i, n-1
       if (model%x(i+1) - model%x(i) <= tiny) cycle
       call poisson_ode_coefficients(model, l, i, 0.0_dp, cache%rk4_coeff(:,1,i))
       call poisson_ode_coefficients(model, l, i, 0.5_dp, cache%rk4_coeff(:,2,i))
       call poisson_ode_coefficients(model, l, i, 1.0_dp, cache%rk4_coeff(:,3,i))
    end do
    cache%prepared = .true.
  end subroutine prepare_poisson_cache


  subroutine solve_poisson_correction(model, l, sol, sources, cache)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    type(mode_solution), intent(inout) :: sol
    type(source_terms), intent(inout) :: sources
    type(poisson_degree_cache), intent(in) :: cache

    if (model%n <= 1) return
    call solve_poisson_shooting_correction(model, l, sol, sources, cache)
  end subroutine solve_poisson_correction


  subroutine solve_poisson_shooting_correction(model, l, sol, sources, cache)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    type(mode_solution), intent(inout) :: sol
    type(source_terms), intent(inout) :: sources
    type(poisson_degree_cache), intent(in) :: cache

    complex(dp), allocatable :: ymech(:,:), particular(:,:), raw_phi(:,:)
    complex(dp) :: bc_particular, bc_homogeneous, amplitude
    real(dp) :: sd, sp, s2
    integer :: i, n, start_i, m

    n = model%n
    if (n <= 1) return
    if (.not. cache%prepared .or. cache%n /= n .or. cache%l /= l) &
         error stop 'SISMO: invalid Poisson degree cache'

    call pack_mechanical_shape(model, l, sol, ymech)
    allocate(particular(2,n), raw_phi(2,n))
    particular = czero
    start_i = 1
    if (model%x(1) <= tiny .and. n >= 2) then
       particular(:,2) = particular(:,1)
       start_i = 2
    end if

    do i = start_i, n-1
       call rk4_poisson_step_cached(model, l, i, cache%rk4_coeff(:,:,i), ymech(:,i), ymech(:,i+1), &
            particular(:,i), particular(:,i+1))
    end do

    ! Match the vacuum exterior condition at the photosphere (interior-atmosphere
    ! raccord) rather than the top of the extended atmosphere: the homogeneous
    ! amplitude this fixes scales phi everywhere, so distorting it by integrating
    ! the growing (~x^l) solution through the thick low-density atmosphere corrupts
    ! phi back in the envelope.  The vacuum match belongs where the star ends.
    m = n
    bc_particular = cmplx(model%rho(m)/max(tiny, model%qx3(m)), 0.0_dp, kind=dp)*ymech(1,m) &
         + cmplx(real(l+1, dp), 0.0_dp, kind=dp)*particular(1,m) + particular(2,m)
    bc_homogeneous = cache%bc_homogeneous
    if (abs(bc_homogeneous) <= tiny) then
       call solve_poisson_integral_correction(model, l, sol, sources)
       return
    end if
    amplitude = -bc_particular/bc_homogeneous
    if (use_takata_closure .and. l == 1) then
       ! Takata J(x) = w2*x*rho*a - P*pp - (w2+2*s2)*phi + x*(w2-s2)*theta = 0.
       ! With pp expressed through the mechanics and phi (pack relation
       ! y2 = y3 + pp*s3/(sp*s2)  =>  -P*pp = -rho*sp*s2*(y2 - y3)), J is linear
       ! in the amplitude A through BOTH the phi and the pp terms:
       !   J_i = M_i + A*H_i
       !   M_i = w2*x*rho*y1 - rho*x*s2*y2 + (rho - w2 - 2*s2)*x*s2*part1
       !         + x*(w2 - s2)*s2*part2
       !   H_i = (rho - w2 - 2*s2)*x*s2*homo1 + x*(w2 - s2)*s2*homo2
       ! (l=1: sd=1, sp=x).  A = -<H,M>/<H,H> (dx-weighted least squares).
       block
         complex(dp) :: mm, hh, num_c, den_c, sig2
         real(dp) :: xw, rhow, s2w, wdx
         integer :: iw
         sig2 = sol%omega*sol%omega
         num_c = czero
         den_c = czero
         do iw = 2, n-1
            xw = model%x(iw)
            if (xw <= tiny) cycle
            rhow = model%rho(iw)
            s2w = max(tiny, model%qx3(iw))
            wdx = 0.5_dp*(model%x(iw+1) - model%x(iw-1))
            mm = sig2*cmplx(xw*rhow, 0.0_dp, kind=dp)*ymech(1,iw) &
                 - cmplx(rhow*xw*s2w, 0.0_dp, kind=dp)*ymech(2,iw) &
                 + (cmplx(rhow, 0.0_dp, kind=dp) - sig2 - cmplx(2.0_dp*s2w, 0.0_dp, kind=dp)) &
                   *cmplx(xw*s2w, 0.0_dp, kind=dp)*particular(1,iw) &
                 + (sig2 - cmplx(s2w, 0.0_dp, kind=dp))*cmplx(xw*s2w, 0.0_dp, kind=dp)*particular(2,iw)
            hh = (cmplx(rhow, 0.0_dp, kind=dp) - sig2 - cmplx(2.0_dp*s2w, 0.0_dp, kind=dp)) &
                   *cmplx(xw*s2w, 0.0_dp, kind=dp)*cache%homogeneous(1,iw) &
                 + (sig2 - cmplx(s2w, 0.0_dp, kind=dp))*cmplx(xw*s2w, 0.0_dp, kind=dp)*cache%homogeneous(2,iw)
            num_c = num_c + cmplx(wdx, 0.0_dp, kind=dp)*conjg(hh)*mm
            den_c = den_c + cmplx(wdx, 0.0_dp, kind=dp)*conjg(hh)*hh
         end do
         if (abs(den_c) > tiny) amplitude = -num_c/den_c
       end block
    end if
    raw_phi = particular + amplitude*cache%homogeneous

    do i = 1, n
       sd = regularization_scale(model%x(i), max(0, l-1))
       sp = regularization_scale(model%x(i), max(0, l))
       s2 = max(tiny, model%qx3(i))
       sol%phi(i) = cmplx(sp*s2, 0.0_dp, kind=dp)*raw_phi(1,i)
       sol%theta(i) = cmplx(sd*s2, 0.0_dp, kind=dp)*raw_phi(2,i)
    end do

    ! NOTE: a vacuum continuation phi ~ x^-(l+1) beyond the photosphere was tried
    ! here and REVERTED -- replacing the atmosphere theta perturbs the near-surface
    ! radial source (sources%radial -= theta below), which destabilised the coupled
    ! iteration for some modes (l=1 n=-6,-7 lost their mechanical eigenfunction).
    ! The interior-ODE continuation through the thin (x=1..1.003) atmosphere is
    ! benign; the photosphere MATCH POINT (m above) is what matters.

    do i = 2, n-1
       sources%radial(i) = sources%radial(i) - sol%theta(i)
       sources%horizontal(i) = sources%horizontal(i) + sol%phi(i)
    end do
    sources%horizontal(n) = sources%horizontal(n) + sol%phi(n)
  end subroutine solve_poisson_shooting_correction

  complex(dp) function cone() result(value)
    value = cmplx(1.0_dp, 0.0_dp, kind=dp)
  end function cone

  subroutine rk4_poisson_step(model, l, ymech, i, y0, include_forcing, y1)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l, i
    complex(dp), intent(in) :: ymech(:,:)
    complex(dp), intent(in) :: y0(2)
    logical, intent(in) :: include_forcing
    complex(dp), intent(out) :: y1(2)

    complex(dp) :: k1(2), k2(2), k3(2), k4(2), ym(2)
    real(dp) :: h

    h = model%x(i+1) - model%x(i)
    if (h <= tiny) then
       y1 = y0
       return
    end if

    ym = 0.5_dp*(ymech(:,i) + ymech(:,i+1))
    call poisson_ode_rhs(model, l, i, 0.0_dp, ymech(:,i), y0, include_forcing, k1)
    call poisson_ode_rhs(model, l, i, 0.5_dp, ym, y0 + cmplx(0.5_dp*h, 0.0_dp, kind=dp)*k1, &
         include_forcing, k2)
    call poisson_ode_rhs(model, l, i, 0.5_dp, ym, y0 + cmplx(0.5_dp*h, 0.0_dp, kind=dp)*k2, &
         include_forcing, k3)
    call poisson_ode_rhs(model, l, i, 1.0_dp, ymech(:,i+1), y0 + cmplx(h, 0.0_dp, kind=dp)*k3, &
         include_forcing, k4)
    y1 = y0 + cmplx(h/6.0_dp, 0.0_dp, kind=dp)*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
  end subroutine rk4_poisson_step

  subroutine rk4_poisson_step_cached(model, l, i, coeff, ymech_left, ymech_right, y0, y1)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l, i
    real(dp), intent(in) :: coeff(6,3)
    complex(dp), intent(in) :: ymech_left(2), ymech_right(2), y0(2)
    complex(dp), intent(out) :: y1(2)

    complex(dp) :: k1(2), k2(2), k3(2), k4(2), ym(2)
    real(dp) :: h

    h = model%x(i+1) - model%x(i)
    if (h <= tiny) then
       y1 = y0
       return
    end if

    ym = 0.5_dp*(ymech_left + ymech_right)
    call poisson_ode_rhs_from_coefficients(coeff(:,1), l, ymech_left, y0, .true., k1)
    call poisson_ode_rhs_from_coefficients(coeff(:,2), l, ym, &
         y0 + cmplx(0.5_dp*h, 0.0_dp, kind=dp)*k1, .true., k2)
    call poisson_ode_rhs_from_coefficients(coeff(:,2), l, ym, &
         y0 + cmplx(0.5_dp*h, 0.0_dp, kind=dp)*k2, .true., k3)
    call poisson_ode_rhs_from_coefficients(coeff(:,3), l, ymech_right, &
         y0 + cmplx(h, 0.0_dp, kind=dp)*k3, .true., k4)
    y1 = y0 + cmplx(h/6.0_dp, 0.0_dp, kind=dp)*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
  end subroutine rk4_poisson_step_cached

  subroutine poisson_ode_rhs(model, l, i, weight, ymech_at, y, include_forcing, dydx)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l, i
    real(dp), intent(in) :: weight
    complex(dp), intent(in) :: ymech_at(2), y(2)
    logical, intent(in) :: include_forcing
    complex(dp), intent(out) :: dydx(2)

    real(dp) :: coeff(6)

    call poisson_ode_coefficients(model, l, i, weight, coeff)
    call poisson_ode_rhs_from_coefficients(coeff, l, ymech_at, y, include_forcing, dydx)
  end subroutine poisson_ode_rhs

  subroutine poisson_ode_coefficients(model, l, i, weight, coeff)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l, i
    real(dp), intent(in) :: weight
    real(dp), intent(out) :: coeff(6)

    real(dp) :: x, s2, s3, s4, s5, s6, c33, c21, c42

    coeff = 0.0_dp
    x = interpolated_scalar(model%x(i), model%x(i+1), weight)
    if (x <= tiny) return
    s2 = max(tiny, interpolated_scalar(model%qx3(i), model%qx3(i+1), weight))
    s4 = interpolated_scalar(model%rho(i), model%rho(i+1), weight)
    s5 = max(tiny, interpolated_scalar(model%gamma1(i), model%gamma1(i+1), weight))
    s3 = max(tiny, interpolated_scalar(model%pressure(i), model%pressure(i+1), weight))/max(tiny, s4)
    s6 = interpolated_scalar(model%aosc(i), model%aosc(i+1), weight)

    c33 = (3.0_dp - s4/s2 - real(l, dp))/x
    c21 = x*s6
    c42 = x*s4/(s3*s5)

    coeff(1) = x
    coeff(2) = s2
    coeff(3) = s4
    coeff(4) = c33
    coeff(5) = c21
    coeff(6) = c42
  end subroutine poisson_ode_coefficients

  subroutine poisson_ode_rhs_from_coefficients(coeff, l, ymech_at, y, include_forcing, dydx)
    real(dp), intent(in) :: coeff(6)
    integer, intent(in) :: l
    complex(dp), intent(in) :: ymech_at(2), y(2)
    logical, intent(in) :: include_forcing
    complex(dp), intent(out) :: dydx(2)

    real(dp) :: ell2, x, s2, s4, c33, c21, c42

    x = coeff(1)
    if (x <= tiny) then
       dydx = czero
       return
    end if

    ell2 = real(l*(l+1), dp)
    s2 = coeff(2)
    s4 = coeff(3)
    c33 = coeff(4)
    c21 = coeff(5)
    c42 = coeff(6)
    dydx(1) = cmplx(c33, 0.0_dp, kind=dp)*y(1) + y(2)/cmplx(x, 0.0_dp, kind=dp)
    dydx(2) = cmplx(ell2/x - c42, 0.0_dp, kind=dp)*y(1) &
         + cmplx(c33 - 1.0_dp/x, 0.0_dp, kind=dp)*y(2)
    if (include_forcing) then
       dydx(2) = dydx(2) - cmplx(c21*s4/s2, 0.0_dp, kind=dp)*ymech_at(1) &
            + cmplx(c42, 0.0_dp, kind=dp)*ymech_at(2)
    end if
  end subroutine poisson_ode_rhs_from_coefficients

  real(dp) function interpolated_scalar(left, right, weight) result(value)
    real(dp), intent(in) :: left, right, weight

    value = (1.0_dp - weight)*left + weight*right
  end function interpolated_scalar

  subroutine pack_mechanical_shape(model, l, sol, ymech)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    type(mode_solution), intent(in) :: sol
    complex(dp), allocatable, intent(out) :: ymech(:,:)

    complex(dp) :: sigma2
    integer :: i
    real(dp) :: sd, s2

    allocate(ymech(2, model%n))
    ymech = czero
    sigma2 = sol%omega*sol%omega
    if (abs(sigma2) <= tiny) sigma2 = cmplx(tiny, 0.0_dp, kind=dp)
    do i = 1, model%n
       sd = regularization_scale(model%x(i), max(0, l-1))
       s2 = max(tiny, model%qx3(i))
       ymech(1,i) = sol%xi_r(i)/cmplx(max(sd, tiny), 0.0_dp, kind=dp)
       ymech(2,i) = sigma2*sol%xi_h(i)/cmplx(max(sd*s2, tiny), 0.0_dp, kind=dp)
    end do
  end subroutine pack_mechanical_shape

  real(dp) function regularization_scale(x, power) result(scale)
    real(dp), intent(in) :: x
    integer, intent(in) :: power

    if (power <= 0) then
       scale = 1.0_dp
    else
       scale = max(x, 0.0_dp)**power
    end if
  end function regularization_scale

  subroutine solve_poisson_integral_correction(model, l, sol, sources)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    type(mode_solution), intent(inout) :: sol
    type(source_terms), intent(inout) :: sources

    integer :: n, i
    complex(dp), allocatable :: rho_source(:), phi(:), dphi(:)

    n = model%n

    allocate(rho_source(n), phi(n), dphi(n))
    rho_source = czero
    phi = czero
    dphi = czero

    do i = 2, n-1
       rho_source(i) = cmplx(model%rho(i), 0.0_dp, kind=dp) * density_contrast(model, sol, l, i)
    end do
    if (n >= 2) then
       if (l == 0) then
          rho_source(1) = rho_source(2)
       else
          rho_source(1) = czero
       end if
       rho_source(n) = rho_source(n-1)
    end if

    call solve_poisson_integral(model, l, rho_source, phi, dphi)
    sol%phi = phi
    sol%theta = dphi

    do i = 2, n-1
       sources%radial(i) = sources%radial(i) - dphi(i)
       sources%horizontal(i) = sources%horizontal(i) + phi(i)
    end do
    if (n >= 2) sources%horizontal(n) = sources%horizontal(n) + phi(n)
  end subroutine solve_poisson_integral_correction

  subroutine solve_poisson_integral(model, l, rho_source, phi, dphi)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l
    complex(dp), intent(in) :: rho_source(:)
    complex(dp), intent(out) :: phi(:), dphi(:)

    complex(dp), allocatable :: inner(:), outer(:)
    real(dp) :: x, dx, factor
    integer :: n, i

    n = model%n
    allocate(inner(n), outer(n))
    inner = czero
    outer = czero

    ! Green-function solution of
    !   (1/x^2) d/dx(x^2 dphi/dx) - l(l+1) phi/x^2 = rho_source,
    ! with the regular center solution and the vacuum exterior solution.
    do i = 2, n
       dx = max(tiny, model%x(i) - model%x(i-1))
       inner(i) = inner(i-1) + cmplx(0.5_dp*dx, 0.0_dp, kind=dp) &
            * (inner_kernel(model%x(i-1), l, rho_source(i-1)) + inner_kernel(model%x(i), l, rho_source(i)))
    end do

    do i = n-1, 1, -1
       dx = max(tiny, model%x(i+1) - model%x(i))
       outer(i) = outer(i+1) + cmplx(0.5_dp*dx, 0.0_dp, kind=dp) &
            * (outer_kernel(model%x(i), l, rho_source(i)) + outer_kernel(model%x(i+1), l, rho_source(i+1)))
    end do

    factor = -1.0_dp/real(2*l + 1, dp)
    do i = 1, n
       x = model%x(i)
       if (x <= tiny) then
          if (l == 0) then
             phi(i) = cmplx(factor, 0.0_dp, kind=dp)*outer(i)
             dphi(i) = czero
          else
             phi(i) = czero
             if (l == 1) then
                dphi(i) = cmplx(factor, 0.0_dp, kind=dp)*outer(i)
             else
                dphi(i) = czero
             end if
          end if
       else
          phi(i) = cmplx(factor, 0.0_dp, kind=dp) &
               * (inner(i)/cmplx(x**(l+1), 0.0_dp, kind=dp) &
               + cmplx(x**l, 0.0_dp, kind=dp)*outer(i))
          if (l == 0) then
             dphi(i) = -cmplx(factor, 0.0_dp, kind=dp)*inner(i)/cmplx(x*x, 0.0_dp, kind=dp)
          else
             dphi(i) = cmplx(factor, 0.0_dp, kind=dp) &
                  * (-cmplx(real(l+1, dp), 0.0_dp, kind=dp)*inner(i)/cmplx(x**(l+2), 0.0_dp, kind=dp) &
                  + cmplx(real(l, dp)*x**(l-1), 0.0_dp, kind=dp)*outer(i))
          end if
       end if
    end do
  end subroutine solve_poisson_integral

  complex(dp) function inner_kernel(x, l, rho_source) result(value)
    real(dp), intent(in) :: x
    integer, intent(in) :: l
    complex(dp), intent(in) :: rho_source

    if (x <= tiny) then
       value = czero
    else
       value = cmplx(x**(l+2), 0.0_dp, kind=dp)*rho_source
    end if
  end function inner_kernel

  complex(dp) function outer_kernel(x, l, rho_source) result(value)
    real(dp), intent(in) :: x
    integer, intent(in) :: l
    complex(dp), intent(in) :: rho_source

    if (x <= tiny) then
       if (l <= 1) then
          value = rho_source
       else
          value = czero
       end if
    else
       value = cmplx(x**(1-l), 0.0_dp, kind=dp)*rho_source
    end if
  end function outer_kernel

  complex(dp) function density_contrast(model, sol, l, i) result(drho)
    type(stellar_model), intent(in) :: model
    type(mode_solution), intent(in) :: sol
    integer, intent(in) :: l, i

    real(dp) :: ell2, x, dx, dlnrho_dx
    complex(dp) :: divxi

    ell2 = real(l*(l+1), dp)
    x = max(model%x(i), tiny)
    dx = max(tiny, model%x(i+1)-model%x(i-1))
    dlnrho_dx = model%drhdxx(i)/max(tiny, model%rho(i))
    divxi = (sol%xi_r(i+1)-sol%xi_r(i-1))/cmplx(dx, 0.0_dp, kind=dp) &
         + cmplx(2.0_dp/x, 0.0_dp, kind=dp)*sol%xi_r(i) &
         - cmplx(ell2/x, 0.0_dp, kind=dp)*sol%xi_h(i)
    ! Poisson is sourced by the Eulerian density perturbation.  The previous
    ! version used only the Lagrangian -div(xi) term, and damped it by an
    ! arbitrary 0.05 factor, which displaced the non-Cowling correction.
    drho = -divxi - cmplx(dlnrho_dx, 0.0_dp, kind=dp)*sol%xi_r(i)
  end function density_contrast

end module sismo_poisson
