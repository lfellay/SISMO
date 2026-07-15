
!   physctes.f95   2019-11-08

!===============================================================================

! The module physctes provides mathematical, physical and astronomical
! constants. The physical constants are the CODATA Internationally recommended
! values of the Fundamental Physical Constants, 2018. The astronomical constants
! come from various sources. CGS units are used unless otherwise specified.

!===============================================================================

module physctes
  implicit none
  ! 0x7FF8000000000000 is a quiet NaN in IEEE 754 double precision.
  ! transfer() from a named integer constant is a valid Fortran constant expression.
  integer(8), parameter :: NaN_bits = 9221120237041090560_8   ! = z'7FF8000000000000'
  real(8),    parameter :: NaN = transfer(NaN_bits, 0d0)
  real(8), parameter :: &
  euler=exp(1d0),                  & ! constante d'Euler = 2.7182818285E+00
pi=acos(-1d0),                   & ! pi = 3.1415926536E+00
cLight=2.99792458d10,            & ! speed of light
electrSI=1.602176634d-19,        & ! elementary charge in SI units
electr=electrSI*cLight/10d0,     & ! elementary charge = 4.8032047126E-10
electrV=1d7*electrSI,            & ! electron-Volt = 1.6021766340E-12
Grav=6.67430d-8,                 & ! Newtonian constant of gravitation
hPlanck=6.62607015d-27,          & ! Planck constant
hbar=hPlanck/2d0/pi,             & ! h/(2pi) = 1.0545718176E-27
kBoltz=1.380649d-16,             & ! Boltzmann constant
umass=1.66053906660d-24,         & ! atomic mass unit (m_u)
melectru=5.48579909065d-4,       & ! electron mass / atomic mass unit
melectr=melectru*umass,          & ! electron mass = 9.1093837015E-28
sigmaStefBoltz=2d0*pi**5/15d0*kBoltz**4/cLight**2/hPlanck**3, &
                                  ! Stefan-Boltzmann constant = 5.6703744192E-05
aRad=4*sigmaStefBoltz/cLight,    & ! radiation constant = 7.5657332503E-15
NAvog=6.02214076d23,             & ! Avogadro number
Rgas=NAvog*kBoltz,               & ! gas constant = 8.3144626182E+07
year=3.15569259747d7,            & ! tropical year 1900, January 0 at 12 h
                                   ! (CGPM 1960, resol. 9)
au=1.495978707d13,               & ! astronomical unit (IAU 2012, resol. B2)
GMsun=1.3271244d26,              & ! solar mass parameter (IAU 2015, resol. B3)
Msun=GMsun/Grav,                 & ! solar mass = 1.9884098707E+33
Lsun=3.828d33,                   & ! solar luminosity (IAU 2015, resol. B3)
Lbol0=3.0128d35,                 & ! luminosity of a Mbol=0 source
                                   ! (IAU 2015, resol. B2)
Mbolsun=-2.5d0*log10(Lsun/Lbol0),& ! solar bolometric magnitude
                                   ! = 4.7399959339E+00
Rsun=6.957d10,                   & ! solar radius (IAU 2015, resol. B3)
Teffsun=5772d0,                  & ! solar effective temperature
                                   ! (IAU 2015, resol. B3)
AgeSunyr=4.572d9,                & ! solar age since birth, in years
                                   ! (Mamajek's Star Notes)
AgeSun=AgeSunyr*year               ! solar age since birth = 1.4427826556E+17

end module physctes
