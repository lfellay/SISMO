module sismo_precision
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 300)
  real(dp), parameter :: tiny = 1.0d-99
  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: czero = cmplx(0.0_dp, 0.0_dp, kind=dp)
  complex(dp), parameter :: cone = cmplx(1.0_dp, 0.0_dp, kind=dp)
  complex(dp), parameter :: ci = cmplx(0.0_dp, 1.0_dp, kind=dp)

end module sismo_precision
