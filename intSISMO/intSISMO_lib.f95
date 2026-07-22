module intSISMO_lib
  use physctes
  use madmodlib
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

  implicit none

  private

  integer, parameter :: UNIT_STDOUT = 6
  real(8), parameter :: TINY_DR = 1d-30
  real(8), parameter :: BV_DENSITY_WEIGHT = 2d0

  public :: write_osc_model

contains

  !---------------------------------------------------------------------
  ! Write the adiabatic OSC input model.
  !
  ! The OSC-42 container is a little-endian stream file with columns:
  !   r, m, P, rho, Gamma1, A
  ! where A = dlnrho/dr - dlnP/(Gamma1 dr).
  !---------------------------------------------------------------------
  subroutine write_osc_model(filename, star, NN1, NN1a, N6, zones_osc, nz_osc, &
       grid_mode, output_ok, points_written)
    implicit none

    character(len=*), intent(in) :: filename
    type(tMadmod), intent(in) :: star
    integer, intent(in) :: NN1, NN1a, N6, nz_osc
    integer, intent(in) :: zones_osc(:,:)
    character(len=*), intent(in), optional :: grid_mode
    logical, intent(out) :: output_ok
    integer, intent(out) :: points_written

    integer :: i, io, ios, ncap, N3, N5_osc, n_osc
    real(8) :: gconst
    real(8) :: osc_globals(3)
    real(8), allocatable :: r21(:), r22(:), d1(:), d2(:)
    real(8), allocatable :: rr1(:), rr21(:), rr22(:)
    real(8), allocatable :: rho21(:), rho22(:), p22(:), gam1(:)
    real(8), allocatable :: mr3(:), a_osc(:), osc_s(:,:)
    integer, allocatable :: osc_z(:,:)
    character(len=32) :: mode
    character(len=512) :: iomsg
    logical :: stage_ok

    output_ok = .false.
    points_written = 0

    if (NN1a < 8 .or. NN1 < 8) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - source model is too small for OSC.'
       return
    end if
    if (N6 < 8) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - invalid OSC target grid size: ', N6
       return
    end if
    if (nz_osc < 1) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - empty OSC source zone table.'
       return
    end if
    if (NN1 /= star%npi .or. NN1a /= star%np) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - source dimensions disagree with the MAD model.'
       write(UNIT_STDOUT,*) 'write_osc_model: NN1/star%npi=', NN1, star%npi, &
            ' NN1a/star%np=', NN1a, star%np
       return
    end if
    if (.not. allocated(star%r) .or. .not. allocated(star%m) .or. &
         .not. allocated(star%rho) .or. .not. allocated(star%P) .or. &
         .not. allocated(star%Gam1) .or. .not. allocated(star%n2)) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - required source arrays are not allocated.'
       return
    end if
    if (size(star%r) < NN1a .or. size(star%m) < NN1a .or. &
         size(star%rho) < NN1a .or. size(star%P) < NN1a .or. &
         size(star%Gam1) < NN1a .or. size(star%n2) < NN1a) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - a required source array is shorter than NN1a.'
       return
    end if
    if (NN1 > NN1a .or. NN1 > size(star%r) .or. NN1-5 > NN1a-3) then
       write(UNIT_STDOUT,*) &
            'write_osc_model: ERROR - interior grid extent is inconsistent with the source arrays.'
       write(UNIT_STDOUT,*) 'write_osc_model: NN1/NN1a/radius-size=', &
            NN1, NN1a, size(star%r)
       return
    end if
    if (size(zones_osc,1) < 2 .or. size(zones_osc,2) < nz_osc) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - zone table dimensions are too small.'
       write(UNIT_STDOUT,*) 'write_osc_model: shape=', size(zones_osc,1), &
            size(zones_osc,2), ' required second dimension=', nz_osc
       return
    end if
    do i = 1, nz_osc
       if (zones_osc(1,i) < 1 .or. zones_osc(2,i) < zones_osc(1,i) .or. &
            zones_osc(2,i) > NN1a) then
          write(UNIT_STDOUT,*) 'write_osc_model: ERROR - invalid zone bounds at zone ', i, &
               ': ', zones_osc(1,i), zones_osc(2,i)
          return
       end if
    end do
    if (.not. ieee_is_finite(star%mass) .or. star%mass <= 0d0 .or. &
         .not. ieee_is_finite(star%radius) .or. star%radius <= 0d0) then
       write(UNIT_STDOUT,*) &
            'write_osc_model: ERROR - stellar mass and radius must be finite and positive.'
       return
    end if
    if (.not. ieee_is_finite(star%Grav)) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - non-finite gravitational constant.'
       return
    end if
    gconst = star%Grav
    if (gconst <= 0d0) gconst = Grav
    if (.not. ieee_is_finite(gconst) .or. gconst <= 0d0) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - invalid effective gravitational constant.'
       return
    end if
    call validate_source_radii(star%r, zones_osc, nz_osc, NN1a, stage_ok)
    if (.not. stage_ok) return
    call validate_source_values(star%m, star%P, star%rho, star%Gam1, star%n2, &
         NN1a, stage_ok)
    if (.not. stage_ok) return

    N3 = NN1a - 3
    ncap = max(N6 + 64, NN1a + 64)
    iomsg = ''
    allocate(r21(ncap), r22(ncap), d1(ncap), d2(ncap), stat=ios, errmsg=iomsg)
    if (ios /= 0) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - cannot allocate grid work arrays: ', &
            trim(iomsg)
       go to 900
    end if

    mode = 'radial'
    if (present(grid_mode)) mode = lowercase_ascii(adjustl(grid_mode))
    select case (trim(mode))
    case ('radial', 'bv')
       continue
    case default
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - grid_mode must be radial or bv: ', &
            trim(mode)
       go to 900
    end select

    call build_osc_grid(star%r, star%n2, star%npi, NN1, N3, N6, mode, &
         r21, r22, d1, d2, N5_osc, stage_ok)
    if (.not. stage_ok) go to 900
    n_osc = N5_osc - 1
    if (n_osc < 3) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - generated OSC grid is too small.'
       go to 900
    end if

    iomsg = ''
    allocate(rr1(NN1a), rr21(N5_osc), rr22(N5_osc), stat=ios, errmsg=iomsg)
    if (ios /= 0) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - cannot allocate coordinate arrays: ', &
            trim(iomsg)
       go to 900
    end if
    iomsg = ''
    allocate(rho21(N5_osc), rho22(0:n_osc), p22(0:n_osc), gam1(0:n_osc), &
         stat=ios, errmsg=iomsg)
    if (ios /= 0) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - cannot allocate thermodynamic arrays: ', &
            trim(iomsg)
       go to 900
    end if
    iomsg = ''
    allocate(mr3(0:n_osc), a_osc(n_osc), osc_s(6,n_osc), osc_z(2,1), &
         stat=ios, errmsg=iomsg)
    if (ios /= 0) then
       write(UNIT_STDOUT,*) 'write_osc_model: ERROR - cannot allocate OSC output arrays: ', &
            trim(iomsg)
       go to 900
    end if

    do i = 1, NN1a
       rr1(i) = star%r(i)**2
    end do
    do i = 1, N5_osc
       rr21(i) = r21(i)**2
       rr22(i) = r22(i)**2
    end do

    call SPLINE1(zones_osc, nz_osc, star%rho, rho21(1:N5_osc), rr1, &
         rr21(1:N5_osc), stage_ok)
    if (.not. stage_ok) go to 900
    do i = 1, max(0, N5_osc-2)
       rho22(i) = 0.5d0 * (rho21(i) + rho21(i+1))
    end do
    rho22(n_osc) = rho21(N5_osc)

    call SPLINE12(zones_osc, nz_osc, star%P, p22, rr1, rr22(1:n_osc), stage_ok)
    if (.not. stage_ok) go to 900
    call SPLINE12(zones_osc, nz_osc, star%Gam1, gam1, rr1, &
         rr22(1:n_osc), stage_ok)
    if (.not. stage_ok) go to 900
    call interpolate_mass_on_osc_grid(star%r, star%m, zones_osc, nz_osc, &
         NN1a, N3, rr1, rr22(1:n_osc), mr3, stage_ok)
    if (.not. stage_ok) go to 900

    call compute_osc_a_profile(n_osc, r22(1:n_osc), p22(1:n_osc), &
         rho22(1:n_osc), gam1(1:n_osc), a_osc)

    do i = 1, n_osc
       osc_s(1,i) = r22(i)
       osc_s(2,i) = mr3(i) * r22(i)**3
       osc_s(3,i) = p22(i)
       osc_s(4,i) = rho22(i)
       osc_s(5,i) = gam1(i)
       osc_s(6,i) = a_osc(i)
    end do

    osc_z(1,1) = 1
    osc_z(2,1) = n_osc

    osc_globals = (/ star%radius, star%mass, gconst /)

    iomsg = ''
    open(newunit=io, file=trim(filename), status='new', action='write', &
         form='unformatted', access='stream', convert='little_endian', &
         iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       write(UNIT_STDOUT,'(A,A,A,I0,A,A)') 'write_osc_model: ERROR - cannot open ', &
            trim(filename), ' (iostat=', ios, '): ', trim(iomsg)
       go to 900
    end if

    write(io, iostat=ios, iomsg=iomsg) n_osc, 1, osc_globals, osc_z, osc_s
    if (ios /= 0) then
       call discard_output_file(filename, io)
       write(UNIT_STDOUT,'(A,A,A,I0,A,A)') 'write_osc_model: ERROR - write failed for ', &
            trim(filename), ' (iostat=', ios, '): ', trim(iomsg)
       go to 900
    end if

    close(io, iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       call discard_output_file(filename, io)
       write(UNIT_STDOUT,'(A,A,A,I0,A,A)') 'write_osc_model: ERROR - close failed for ', &
            trim(filename), ' (iostat=', ios, '): ', trim(iomsg)
       go to 900
    end if

    points_written = n_osc
    output_ok = .true.

900 continue
    if (allocated(r21)) deallocate(r21)
    if (allocated(r22)) deallocate(r22)
    if (allocated(d1)) deallocate(d1)
    if (allocated(d2)) deallocate(d2)
    if (allocated(rr1)) deallocate(rr1)
    if (allocated(rr21)) deallocate(rr21)
    if (allocated(rr22)) deallocate(rr22)
    if (allocated(rho21)) deallocate(rho21)
    if (allocated(rho22)) deallocate(rho22)
    if (allocated(p22)) deallocate(p22)
    if (allocated(gam1)) deallocate(gam1)
    if (allocated(mr3)) deallocate(mr3)
    if (allocated(a_osc)) deallocate(a_osc)
    if (allocated(osc_s)) deallocate(osc_s)
    if (allocated(osc_z)) deallocate(osc_z)
  end subroutine write_osc_model


  subroutine validate_source_radii(radius, zones, nz, n, values_ok)
    implicit none

    real(8), intent(in) :: radius(:)
    integer, intent(in) :: zones(:,:), nz, n
    logical, intent(out) :: values_ok

    integer :: i, k, first, last

    values_ok = .false.
    do i = 1, n
       if (.not. ieee_is_finite(radius(i)) .or. radius(i) < 0d0) then
          write(UNIT_STDOUT,*) 'write_osc_model: ERROR - invalid radius at source point ', i, ': ', radius(i)
          return
       end if
    end do
    ! A duplicated interface radius is valid when it separates two zones.
    ! DERIVEE only receives triples within one zone, where the abscissae must
    ! be strictly ordered to keep every denominator non-zero.
    do k = 1, nz
       first = zones(1,k)
       last = zones(2,k)
       do i = first + 1, last
          if (radius(i) <= radius(i-1)) then
             write(UNIT_STDOUT,*) &
                  'write_osc_model: ERROR - radii must increase strictly within each source zone.'
             write(UNIT_STDOUT,*) 'write_osc_model: zone/points/radii=', k, &
                  i-1, radius(i-1), i, radius(i)
             return
          end if
       end do
    end do
    values_ok = .true.
  end subroutine validate_source_radii


  subroutine validate_source_values(mass, pressure, density, gamma1, n2, n, &
       values_ok)
    implicit none

    real(8), intent(in) :: mass(:), pressure(:), density(:), gamma1(:), n2(:)
    integer, intent(in) :: n
    logical, intent(out) :: values_ok

    integer :: i

    values_ok = .false.
    do i = 1, n
       if (.not. ieee_is_finite(mass(i)) .or. mass(i) < 0d0) then
          write(UNIT_STDOUT,*) &
               'write_osc_model: ERROR - invalid enclosed mass at source point ', i, ': ', mass(i)
          return
       end if
       if (i > 1) then
          if (mass(i) < mass(i-1)) then
             write(UNIT_STDOUT,*) &
                  'write_osc_model: ERROR - enclosed mass decreases at source point ', i, &
                  ': ', mass(i-1), mass(i)
             return
          end if
       end if
       if (.not. ieee_is_finite(pressure(i)) .or. pressure(i) <= 0d0) then
          write(UNIT_STDOUT,*) &
               'write_osc_model: ERROR - invalid pressure at source point ', i, ': ', pressure(i)
          return
       end if
       if (.not. ieee_is_finite(density(i)) .or. density(i) <= 0d0) then
          write(UNIT_STDOUT,*) &
               'write_osc_model: ERROR - invalid density at source point ', i, ': ', density(i)
          return
       end if
       if (.not. ieee_is_finite(gamma1(i)) .or. gamma1(i) <= 0d0) then
          write(UNIT_STDOUT,*) &
               'write_osc_model: ERROR - invalid Gamma1 at source point ', i, ': ', gamma1(i)
          return
       end if
       if (.not. ieee_is_finite(n2(i))) then
          write(UNIT_STDOUT,*) &
               'write_osc_model: ERROR - non-finite Brunt-Vaisala value at source point ', i, ': ', n2(i)
          return
       end if
    end do
    values_ok = .true.
  end subroutine validate_source_values


  function lowercase_ascii(text) result(lower)
    implicit none

    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code

    lower = text
    do i = 1, len(text)
       code = iachar(lower(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) then
          lower(i:i) = achar(code + iachar('a') - iachar('A'))
       end if
    end do
  end function lowercase_ascii


  subroutine discard_output_file(filename, io)
    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(in) :: io

    integer :: cleanup_ios
    logical :: exists, opened

    inquire(unit=io, opened=opened)
    if (opened) close(io, status='delete', iostat=cleanup_ios)

    inquire(file=trim(filename), exist=exists)
    if (.not. exists) return
    call unlink(trim(filename), cleanup_ios)
  end subroutine discard_output_file


  subroutine build_osc_grid(r1, n2, npi, NN1, N3, N6, grid_mode, r21, r22, &
       d1, d2, N5, grid_ok)
    implicit none

    real(8), intent(in) :: r1(:), n2(:)
    integer, intent(in) :: npi, NN1, N3, N6
    character(len=*), intent(in) :: grid_mode
    real(8), intent(out) :: r21(:), r22(:), d1(:), d2(:)
    integer, intent(out) :: N5
    logical, intent(out) :: grid_ok

    integer :: i, j, alloc_stat
    real(8) :: total_density, target
    real(8), allocatable :: rrep(:), dens(:), cdf(:)
    logical :: use_bv, cdf_ok, density_ok
    character(len=512) :: alloc_msg

    grid_ok = .false.
    if (N3 < 8 .or. NN1 < 8) then
       write(UNIT_STDOUT,*) 'build_osc_grid: ERROR - source grid is too small.'
       return
    end if
    if (N3 > size(r1) .or. NN1 > size(r1) .or. NN1-5 > N3) then
       write(UNIT_STDOUT,*) 'build_osc_grid: ERROR - source-grid extents are inconsistent.'
       write(UNIT_STDOUT,*) 'build_osc_grid: N3/NN1/radius-size=', N3, NN1, size(r1)
       return
    end if
    if (size(r21) < N6 + 2 .or. size(r22) < N6 + 1 .or. &
         size(d1) < N6 + 1 .or. size(d2) < N6 + 1) then
       write(UNIT_STDOUT,*) 'build_osc_grid: ERROR - work arrays are too small.'
       return
    end if

    alloc_msg = ''
    allocate(rrep(N3), dens(N3-1), cdf(N3), stat=alloc_stat, errmsg=alloc_msg)
    if (alloc_stat /= 0) then
       write(UNIT_STDOUT,*) 'build_osc_grid: ERROR - cannot allocate grid-density arrays: ', &
            trim(alloc_msg)
       if (allocated(rrep)) deallocate(rrep)
       if (allocated(dens)) deallocate(dens)
       if (allocated(cdf)) deallocate(cdf)
       return
    end if
    r21 = 0d0
    r22 = 0d0
    d1 = 0d0
    d2 = 0d0

    do i = 1, NN1-5
       if (i < 5) then
          rrep(i) = dble(i+4) / r1(i+5)
       else
          rrep(i) = 9d0 / (r1(i+5) - r1(i-4))
       end if
    end do
    do i = NN1-4, N3
       rrep(i) = rrep(NN1-5)
    end do
    do i = 1, N3-1
       dens(i) = max(TINY_DR, rrep(i))
    end do

    use_bv = trim(grid_mode) == 'bv'
    if (use_bv) then
       call apply_bv_density(n2, npi, N3, dens, density_ok)
       if (.not. density_ok) then
          deallocate(rrep, dens, cdf)
          return
       end if
    end if

    call build_density_cdf(r1, dens, N3, cdf, total_density, cdf_ok)
    if (.not. cdf_ok) then
       deallocate(rrep, dens, cdf)
       return
    end if
    if (total_density <= TINY_DR) then
       write(UNIT_STDOUT,*) 'build_osc_grid: ERROR - non-positive integrated grid density.'
       deallocate(rrep, dens, cdf)
       return
    end if

    N5 = N6 + 1
    r21(1) = 0d0
    do j = 2, N5
       target = total_density * dble(j - 1) / dble(N6)
       r21(j) = invert_density_cdf(r1, cdf, N3, target)
    end do
    r21(N5) = r1(N3)
    r21(N5+1) = 2d0 * r21(N5)

    do i = 1, N5-1
       d1(i) = max(TINY_DR, r21(i+1) - r21(i))
    end do
    d1(N5) = 1d50

    d2(1) = 0d0
    do i = 1, N5-2
       r22(i) = 0.5d0 * (r21(i) + r21(i+1))
       d2(i+1) = 0.5d0 * (d1(i) + d1(i+1))
    end do
    d2(N5-1) = d2(N5-1) + 0.5d0 * d1(N5-1)
    r22(N5-1) = r21(N5)
    r22(N5) = r21(N5+1)

    deallocate(rrep, dens, cdf)
    grid_ok = .true.
  end subroutine build_osc_grid


  subroutine apply_bv_density(n2, npi, N3, dens, density_ok)
    implicit none

    real(8), intent(in) :: n2(:)
    integer, intent(in) :: npi, N3
    real(8), intent(inout) :: dens(:)
    logical, intent(out) :: density_ok

    integer :: i, imax, alloc_stat
    real(8) :: bv, bvmax
    real(8), allocatable :: bv_interval(:)
    character(len=512) :: alloc_msg

    density_ok = .false.
    alloc_msg = ''
    allocate(bv_interval(size(dens)), stat=alloc_stat, errmsg=alloc_msg)
    if (alloc_stat /= 0) then
       write(UNIT_STDOUT,*) 'apply_bv_density: ERROR - cannot allocate BV work array: ', &
            trim(alloc_msg)
       return
    end if
    bv_interval = 0d0
    bvmax = 0d0
    imax = min(N3-1, size(n2)-1, max(0, npi-1))

    do i = 1, imax
       bv = sqrt(0.5d0 * (max(0d0, n2(i)) + max(0d0, n2(i+1))))
       bv_interval(i) = bv
       bvmax = max(bvmax, bv)
    end do

    if (bvmax > TINY_DR) then
       do i = 1, imax
          dens(i) = dens(i) * (1d0 + BV_DENSITY_WEIGHT * bv_interval(i) / bvmax)
       end do
    end if

    deallocate(bv_interval)
    density_ok = .true.
  end subroutine apply_bv_density


  subroutine build_density_cdf(r1, dens, N3, cdf, total_density, cdf_ok)
    implicit none

    real(8), intent(in) :: r1(:), dens(:)
    integer, intent(in) :: N3
    real(8), intent(out) :: cdf(:), total_density
    logical, intent(out) :: cdf_ok

    integer :: i
    real(8) :: dr

    cdf = 0d0
    total_density = 0d0
    cdf_ok = .false.
    do i = 1, N3-1
       dr = r1(i+1) - r1(i)
       if (dr < -TINY_DR) then
          write(UNIT_STDOUT,*) 'build_density_cdf: ERROR - source radius grid decreases at i=', i
          return
       else if (dr <= TINY_DR) then
          cdf(i+1) = cdf(i)
          cycle
       end if
       cdf(i+1) = cdf(i) + max(TINY_DR, dens(i)) * dr
    end do
    total_density = cdf(N3)
    cdf_ok = .true.
  end subroutine build_density_cdf


  real(8) function invert_density_cdf(r1, cdf, N3, target) result(radius)
    implicit none

    real(8), intent(in) :: r1(:), cdf(:), target
    integer, intent(in) :: N3

    integer :: i
    real(8) :: frac, span

    if (target <= 0d0) then
       radius = r1(1)
       return
    end if

    do i = 1, N3-1
       if (target <= cdf(i+1) .and. cdf(i+1) > cdf(i) + TINY_DR) then
          span = cdf(i+1) - cdf(i)
          if (span <= TINY_DR) then
             radius = r1(i)
          else
             frac = (target - cdf(i)) / span
             radius = r1(i) + frac * (r1(i+1) - r1(i))
          end if
          return
       end if
    end do

    radius = r1(N3)
  end function invert_density_cdf


  subroutine interpolate_mass_on_osc_grid(r1, m1, zones, nz, N1, N3, rr1, &
       rr22, mr3, interpolation_ok)
    implicit none

    real(8), intent(in) :: r1(:), m1(:), rr1(:), rr22(:)
    integer, intent(in) :: zones(:,:), nz, N1, N3
    real(8), intent(out) :: mr3(0:)
    logical, intent(out) :: interpolation_ok

    integer :: i, nout, alloc_stat
    real(8), allocatable :: mr31(:)
    character(len=512) :: alloc_msg

    interpolation_ok = .false.
    if (N1 < 2 .or. N1 > size(r1) .or. N1 > size(m1) .or. &
         N1 > size(rr1) .or. N3+1 > N1) then
       write(UNIT_STDOUT,*) 'interpolate_mass_on_osc_grid: ERROR - inconsistent source extents.'
       return
    end if
    alloc_msg = ''
    allocate(mr31(N1), stat=alloc_stat, errmsg=alloc_msg)
    if (alloc_stat /= 0) then
       write(UNIT_STDOUT,*) &
            'interpolate_mass_on_osc_grid: ERROR - cannot allocate mass work array: ', &
            trim(alloc_msg)
       return
    end if
    interpolation_ok = .true.
    mr31 = 0d0
    mr3 = 0d0

    do i = 2, min(N3+1, N1)
       if (abs(r1(i)) > TINY_DR) then
          mr31(i) = m1(i) / r1(i)**3
       else
          mr31(i) = 0d0
       end if
    end do
    if (N1 >= 2) mr31(1) = mr31(2)

    nout = min(size(rr22), ubound(mr3, 1))
    if (nout > 0) then
       call SPLINE12(zones, nz, mr31, mr3(0:nout), rr1, rr22(1:nout), &
            interpolation_ok)
    end if

    deallocate(mr31)
  end subroutine interpolate_mass_on_osc_grid


  subroutine compute_osc_a_profile(n, r, p, rho, gam1, a)
    implicit none

    integer, intent(in) :: n
    real(8), intent(in) :: r(:), p(:), rho(:), gam1(:)
    real(8), intent(out) :: a(:)

    integer :: i, ilo, ihi
    real(8) :: dr, dlnrho_dr, dlnp_dr

    a = 0d0
    if (n < 2) return

    do i = 1, n
       ilo = max(1, i-1)
       ihi = min(n, i+1)
       dr = r(ihi) - r(ilo)
       if (dr <= TINY_DR .or. p(ilo) <= TINY_DR .or. p(ihi) <= TINY_DR .or. &
            rho(ilo) <= TINY_DR .or. rho(ihi) <= TINY_DR .or. gam1(i) <= TINY_DR) then
          a(i) = 0d0
       else
          dlnrho_dr = (log(rho(ihi)) - log(rho(ilo))) / dr
          dlnp_dr = (log(p(ihi)) - log(p(ilo))) / dr
          a(i) = dlnrho_dr - dlnp_dr / gam1(i)
       end if
    end do
  end subroutine compute_osc_a_profile


  function DERIVEE(x, x1, x2, x3, f1, f2, f3)
    implicit none
    real(8), intent(in) :: x, x1, x2, x3, f1, f2, f3
    real(8) :: DERIVEE

    DERIVEE = (2d0*x - x2 - x3) / ((x1 - x2) * (x1 - x3)) * f1 &
         + (2d0*x - x1 - x3) / ((x2 - x1) * (x2 - x3)) * f2 &
         + (2d0*x - x1 - x2) / ((x3 - x1) * (x3 - x2)) * f3
  end function DERIVEE


  subroutine SPLINE12(zones, nz, f, fi, x1, x2, interpolation_ok)
    implicit none

    integer, intent(in) :: nz
    integer, intent(in) :: zones(:,:)
    real(8), intent(in) :: f(:), x1(:), x2(:)
    real(8), intent(out) :: fi(0:)
    logical, intent(out) :: interpolation_ok

    integer :: N2
    real(8) :: df(size(f))

    fi = 0d0
    df = 0d0
    interpolation_ok = .true.
    N2 = min(size(x2), ubound(fi, 1))
    if (N2 < 1) then
       if (ubound(fi, 1) >= 0) fi(0) = f(1)
       return
    end if

    call spline_compute_df(zones, nz, f, x1, df)
    fi(0) = f(1)
    call spline_hermite_fill(zones, nz, f, x1, x2, df, fi(1:N2), N2, &
         'SPLINE12', interpolation_ok)
  end subroutine SPLINE12


  subroutine SPLINE1(zones, nz, f, fi, x1, x2, interpolation_ok)
    implicit none

    integer, intent(in) :: nz
    integer, intent(in) :: zones(:,:)
    real(8), intent(in) :: f(:), x1(:), x2(:)
    real(8), intent(out) :: fi(:)
    logical, intent(out) :: interpolation_ok

    integer :: N2
    real(8) :: df(size(f))

    fi = 0d0
    df = 0d0
    interpolation_ok = .true.
    N2 = min(size(x2), size(fi))
    if (N2 < 1) return

    call spline_compute_df(zones, nz, f, x1, df)
    call spline_hermite_fill(zones, nz, f, x1, x2, df, fi(1:N2), N2, &
         'SPLINE1', interpolation_ok)
  end subroutine SPLINE1


  subroutine spline_compute_df(zones, nz, f, x1, df)
    implicit none

    integer, intent(in) :: nz
    integer, intent(in) :: zones(:,:)
    real(8), intent(in) :: f(:), x1(:)
    real(8), intent(out) :: df(:)

    integer :: i, k, i_start, i_end

    df = 0d0
    do k = 1, nz
       i_start = max(1, zones(1,k))
       i_end = min(size(f), zones(2,k))
       if (i_end - i_start < 2) cycle

       i = i_start
       df(i) = DERIVEE(x1(i), x1(i), x1(i+1), x1(i+2), f(i), f(i+1), f(i+2))

       do i = i_start+1, i_end-1
          df(i) = DERIVEE(x1(i), x1(i-1), x1(i), x1(i+1), f(i-1), f(i), f(i+1))
       end do

       i = i_end
       df(i) = DERIVEE(x1(i), x1(i-2), x1(i-1), x1(i), f(i-2), f(i-1), f(i))
    end do
  end subroutine spline_compute_df


  subroutine spline_hermite_fill(zones, nz, f, x1, x2, df, fi, N2, name, &
       interpolation_ok)
    implicit none

    integer, intent(in) :: nz, N2
    integer, intent(in) :: zones(:,:)
    real(8), intent(in) :: f(:), x1(:), x2(:), df(:)
    real(8), intent(out) :: fi(:)
    character(len=*), intent(in) :: name
    logical, intent(out) :: interpolation_ok

    integer :: i, j, k, i_start, i_end, fill_start
    real(8) :: dx, dx2, dx3, epsu

    epsu = 1d-12 * max(1d0, abs(x1(size(x1))))
    interpolation_ok = .false.

    j = 1
    do k = 1, nz
       i_start = max(1, zones(1,k))
       fill_start = i_start
       if (k > 1) fill_start = max(1, i_start-1)
       i_end = min(size(f)-1, zones(2,k)-1)
       do i = fill_start, i_end
          do
             if (j > N2) exit
             if (x2(j) > x1(i+1)) exit
             if (x2(j) < x1(i) - epsu) then
                write(UNIT_STDOUT,*) trim(name)//': ERROR - target not covered by zone intervals.'
                write(UNIT_STDOUT,'(A,I0,3(1X,1PE16.8))') &
                     trim(name)//':        j/x2(j)/interval= ', j, x2(j), x1(i), x1(i+1)
                return
             end if
             dx = x1(i+1) - x1(i)
             if (abs(dx) < TINY_DR) then
                fi(j) = 0.5d0 * (f(i) + f(i+1))
             else
                dx2 = dx * dx
                dx3 = dx2 * dx
                fi(j) = (x2(j)-x1(i)) * (x2(j)-x1(i+1))**2 * df(i) / dx2 &
                     + (x2(j)-x1(i))**2 * (x2(j)-x1(i+1)) * df(i+1) / dx2 &
                     + (x2(j)-x1(i+1))**2 * f(i) * (2d0*x2(j)-3d0*x1(i)+x1(i+1)) / dx3 &
                     - (x2(j)-x1(i))**2 * f(i+1) * (2d0*x2(j)-3d0*x1(i+1)+x1(i)) / dx3
             end if
             j = j + 1
          end do
          if (j > N2) exit
       end do
       if (j > N2) exit
    end do

    if (j <= N2) then
       write(UNIT_STDOUT,*) trim(name)//': ERROR - target grid extends beyond zone coverage.'
       write(UNIT_STDOUT,'(A,I0,1X,1PE16.8)') &
            trim(name)//':        first uncovered j/x2(j)= ', j, x2(j)
       return
    end if
    interpolation_ok = .true.
  end subroutine spline_hermite_fill

end module intSISMO_lib
