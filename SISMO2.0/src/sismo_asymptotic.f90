module sismo_asymptotic
  use sismo_precision, only : dp, pi, tiny
  use sismo_types, only : stellar_model
  implicit none
  private

  public :: initial_asymptotic_g_frequency

contains



  real(dp) function initial_asymptotic_g_frequency(model, l, g_order, eps_g) result(omega)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l, g_order
    real(dp), intent(in) :: eps_g

    real(dp) :: integral, denom

    omega = 0.0_dp
    if (l < 1 .or. g_order < 1) return

    integral = brunt_integral(model, l)
    denom = pi*(real(g_order, dp) + eps_g)
    if (integral > tiny .and. denom > tiny) then
       omega = sqrt(real(l*(l+1), dp))*integral/denom
    end if
  end function initial_asymptotic_g_frequency

  real(dp) function brunt_integral(model, l) result(integral)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: l

    integer :: i, imax
    logical :: use_imported_aosc
    real(dp) :: xmid, dx, n2l, n2r, n2mid, nmid, kl, kr

    integral = 0.0_dp
    if (model%n < 2 .or. l < 1) return

    imax = model%n
    if (model%nrac > 2) imax = min(imax, model%nrac)
    use_imported_aosc = allocated(model%aosc) .and. any(abs(model%aosc(1:imax)) > tiny)

    do i = 2, imax
       dx = model%x(i) - model%x(i-1)
       if (dx <= tiny) cycle

       if (use_imported_aosc) then
          kl = sqrt(max(0.0_dp, model%qx3(i-1)*abs(model%aosc(i-1))))
          kr = sqrt(max(0.0_dp, model%qx3(i)*abs(model%aosc(i))))
          if (kl <= tiny .or. kr <= tiny) cycle
          integral = integral + 0.5_dp*(kl+kr)*dx
       else
          xmid = 0.5_dp*(model%x(i) + model%x(i-1))
          if (xmid <= tiny) cycle

          n2l = brunt_squared(model, i-1)
          n2r = brunt_squared(model, i)
          if (n2l <= tiny .or. n2r <= tiny) cycle
          n2mid = 0.5_dp*(n2l+n2r)
          nmid = sqrt(max(tiny, n2mid))

          integral = integral + nmid*dx/xmid
       end if
    end do
  end function brunt_integral

  real(dp) function brunt_squared(model, i) result(n2)
    type(stellar_model), intent(in) :: model
    integer, intent(in) :: i

    real(dp) :: dlnrho_dx, dlnp_dx, a_x, x, rho, pressure, gamma1, qx3

    n2 = 0.0_dp
    if (i < 1 .or. i > model%n) return
    x = max(tiny, model%x(i))
    rho = max(tiny, model%rho(i))
    pressure = max(tiny, model%pressure(i))
    gamma1 = max(tiny, model%gamma1(i))
    qx3 = max(tiny, model%qx3(i))

    dlnrho_dx = log_derivative(model%x, model%rho, i)
    dlnp_dx = log_derivative(model%x, model%pressure, i)
    if (abs(dlnp_dx) <= tiny) dlnp_dx = -rho*qx3*x/pressure

    a_x = dlnrho_dx - dlnp_dx/gamma1
    n2 = -a_x*qx3*x
    if (n2 < 0.0_dp) n2 = 0.0_dp
  end function brunt_squared

  real(dp) function log_derivative(x, y, i) result(d)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: i

    integer :: ilo, ihi, n
    real(dp) :: dx

    n = size(x)
    d = 0.0_dp
    if (n < 2 .or. i < 1 .or. i > n) return

    ilo = max(1, i-1)
    ihi = min(n, i+1)
    if (ihi == ilo) return
    if (y(ilo) <= tiny .or. y(ihi) <= tiny) return

    dx = x(ihi) - x(ilo)
    if (abs(dx) <= tiny) return
    d = (log(y(ihi)) - log(y(ilo)))/dx
  end function log_derivative







end module sismo_asymptotic
