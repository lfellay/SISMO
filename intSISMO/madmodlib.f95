
!   madmodlib.f95   2018-07-26  (modernised 2026-06)

!===============================================================================

! Le module madmodlib lit un fichier *.madmod produit par cles2mad.
! Il offre le type tMadmod, la procédure readmadmod et DeleteMadmod.

!===============================================================================

module madmodlib
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter :: MAX_MODEL_POINTS = 2000000
  character(len=8), parameter :: COMPACT_MAGIC = 'SISMOAD2'

  public :: tMadmod, readmadmod, DeleteMadmod

  type tMadmod
    integer :: npi = 0, npa = 0, np = 0, nz = 0, step = 0
    logical :: compact_adiabatic = .false.
    logical :: nonadiabatic_available = .false.
    integer,  allocatable :: zones(:,:)
    real(8) :: mass = 0d0, radius = 0d0, Teff = 0d0, Lum = 0d0
    real(8) :: age = 0d0, logg = 0d0, LgOverL = 0d0, alphaConv = 0d0
    real(8) :: alphaOver = 0d0, X0 = 0d0, Z0 = 0d0, diff = 0d0
    real(8) :: Grav = 0d0, aRad = 0d0, sigmaSB = 0d0
    ! Interior + atmosphere fields (1:np)
    real(8), allocatable, dimension(:) :: r, m, &
      rho, T, L, P, Cv, &
      Gam1, Gam31, Prho, PT, &
      Cp, Cprho, CpT, Q, Qrho, &
      QT, kappa, krho, kT, en, &
      X, Z, gradP, gradT, &
      gradRad, gradAd, grad, GamC, &
      AlphaC, dm1, dm2, yH1, &
      yH2, yHe3, yHe4, yLi7, &
      yC12, yC13, yN14, yN15, &
      yO16, yO17, yO18, yNe20, &
      yOthers, ZOthers, ZZOthers, &
      MOthers, eg, mix, n2
    ! Atmosphere-only fields (npi+1:np); only allocated when npa > 0
    real(8), allocatable, dimension(:) :: tau, rr, &
      mm, Pg, Cvg, Gam1g, &
      Gam31g, Prhog, PTg, Cpg, &
      Cprhog, CpTg, Qg, Qrhog, &
      QTg, accrad, dTdTe, dTdg, &
      PR, gradPR
  end type tMadmod

contains

  subroutine readmadmod(fn, star, require_nonadiabatic)
    character(len=*), intent(in) :: fn
    type(tMadmod),    intent(out) :: star
    logical, intent(in), optional :: require_nonadiabatic

    integer :: io, np, npi, npa, nz, np1, ios, alloc_stat
    integer(int64) :: file_size, expected_size, integer_bytes, real_bytes
    character(len=512) :: iomsg, alloc_msg
    character(len=8) :: magic
    logical :: require_full

    require_full = .false.
    if (present(require_nonadiabatic)) require_full = require_nonadiabatic

    open(newunit=io, file=fn, action='read', status='old', &
         form='unformatted', access='stream', convert='little_endian', &
         iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - cannot open ', trim(fn), &
            ' (iostat=', ios, '): ', trim(iomsg)
       stop 1
    end if

    inquire(unit=io, size=file_size, iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - cannot determine size of ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       stop 1
    end if
    if (file_size < len(COMPACT_MAGIC)) then
       write(*,'(A,A)') 'readmadmod: ERROR - model file is too short: ', trim(fn)
       close(io)
       stop 1
    end if

    read(io, pos=1, iostat=ios, iomsg=iomsg) magic
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - format probe failed in ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       stop 1
    end if
    if (magic == COMPACT_MAGIC) then
       call read_compact_model(io, fn, file_size, star)
       if (require_full) then
          write(*,'(A)') 'readmadmod: ERROR - --nonad requires a full legacy non-adiabatic .madmod.'
          write(*,'(A,A)') 'readmadmod: compact adiabatic input: ', trim(fn)
          call DeleteMadmod(star)
          stop 1
       end if
       return
    end if

    read(io, pos=1, iostat=ios, iomsg=iomsg) star%npi, star%npa, star%np, star%nz, star%step
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - header read failed in ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       stop 1
    end if

    npi = star%npi
    npa = star%npa
    np  = star%np
    nz  = star%nz

    if (npi < 1 .or. npa < 0 .or. np < 1 .or. nz < 1 .or. nz > np .or. &
         int(npi, int64) + int(npa, int64) /= int(np, int64)) then
       write(*,'(A,A)') 'readmadmod: ERROR - invalid dimensions in ', trim(fn)
       write(*,'(A,I0,A,I0,A,I0,A,I0)') 'readmadmod: npi=', npi, &
            ' npa=', npa, ' np=', np, ' nz=', nz
       write(*,'(A)') 'readmadmod: expected npi >= 1, npa >= 0, np = npi + npa, and 1 <= nz <= np.'
       close(io)
       stop 1
    end if
    if (np > MAX_MODEL_POINTS) then
       write(*,'(A,A,A,I0,A,I0)') 'readmadmod: ERROR - model in ', trim(fn), &
            ' has ', np, ' points; safety limit is ', MAX_MODEL_POINTS
       close(io)
       stop 1
    end if

    integer_bytes = int(storage_size(star%npi)/8, int64)
    real_bytes = int(storage_size(star%mass)/8, int64)
    expected_size = 5_int64*integer_bytes + 15_int64*real_bytes + &
         3_int64*int(nz,int64)*integer_bytes + &
         52_int64*int(np,int64)*real_bytes + &
         20_int64*int(npa,int64)*real_bytes
    if (file_size < expected_size) then
       write(*,'(A,A)') 'readmadmod: ERROR - file is too short for its declared dimensions: ', trim(fn)
       write(*,'(A,I0,A,I0)') 'readmadmod: file bytes=', file_size, &
            ' minimum expected bytes=', expected_size
       close(io)
       stop 1
    end if

    read(io, iostat=ios, iomsg=iomsg) star%mass, star%radius, star%Teff, &
         star%Lum, star%age, star%logg, star%LgOverL, star%alphaConv, &
         star%alphaOver, star%X0, star%Z0, star%diff, star%Grav, star%aRad, &
         star%sigmaSB
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - global scalar read failed in ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       stop 1
    end if

    alloc_msg = ''
    allocate( &
         star%zones(3,nz), &
         star%r(np),       star%m(np),       star%rho(np),     star%T(np),   &
         star%L(np),       star%P(np),        star%Cv(np),      star%Gam1(np), &
         star%Gam31(np),   star%Prho(np),     star%PT(np),      star%Cp(np),  &
         star%Cprho(np),   star%CpT(np),      star%Q(np),       star%Qrho(np), &
         star%QT(np),      star%kappa(np),    star%krho(np),    star%kT(np),  &
         star%en(np),      star%X(np),        star%Z(np),       star%gradP(np), &
         star%gradT(np),   star%gradRad(np),  star%gradAd(np),  star%grad(np), &
         star%GamC(np),    star%AlphaC(np),   star%dm1(np),     star%dm2(np), &
         star%yH1(np),     star%yH2(np),      star%yHe3(np),    star%yHe4(np), &
         star%yLi7(np),    star%yC12(np),     star%yC13(np),    star%yN14(np), &
         star%yN15(np),    star%yO16(np),     star%yO17(np),    star%yO18(np), &
         star%yNe20(np),   star%yOthers(np),  star%ZOthers(np), star%ZZOthers(np), &
         star%MOthers(np), star%eg(np),       star%mix(np),     star%n2(np), &
         stat=alloc_stat, errmsg=alloc_msg)
    if (alloc_stat /= 0) then
       write(*,'(A,A,A,A)') 'readmadmod: ERROR - allocation failed for ', &
            trim(fn), ': ', trim(alloc_msg)
       close(io)
       call DeleteMadmod(star)
       stop 1
    end if

    read(io, iostat=ios, iomsg=iomsg) &
         star%zones, star%r, star%m, star%rho, star%T, star%L, star%P, &
         star%Cv, star%Gam1, star%Gam31, star%Prho, star%PT, star%Cp, star%Cprho, &
         star%CpT, star%Q, star%Qrho, star%QT, star%kappa, star%krho, star%kT, star%en, &
         star%X, star%Z, star%gradP, star%gradT, star%gradRad, star%gradAd, star%grad, &
         star%GamC, star%AlphaC, star%dm1, star%dm2, star%yH1, star%yH2, star%yHe3, &
         star%yHe4, star%yLi7, star%yC12, star%yC13, star%yN14, star%yN15, star%yO16, &
         star%yO17, star%yO18, star%yNe20, star%yOthers, star%ZOthers, star%ZZOthers, &
         star%MOthers, star%eg, star%mix, star%n2
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - array data read failed in ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       call DeleteMadmod(star)
       stop 1
    end if

    if (npa == 0) then
       star%compact_adiabatic = .false.
       star%nonadiabatic_available = valid_nonadiabatic_block(star)
       if (require_full .and. .not. star%nonadiabatic_available) then
          write(*,'(A)') 'readmadmod: ERROR - --nonad requires populated non-adiabatic fields.'
          write(*,'(A,A)') 'readmadmod: neutral or incomplete legacy input: ', trim(fn)
          close(io)
          call DeleteMadmod(star)
          stop 1
       end if
       close(io, iostat=ios, iomsg=iomsg)
       if (ios /= 0) then
          write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - close failed for ', &
               trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
          call DeleteMadmod(star)
          stop 1
       end if
       return
    end if

    np1 = npi + 1
    alloc_msg = ''
    allocate( &
         star%tau(np1:np),    star%rr(np1:np),    star%mm(np1:np),  &
         star%Pg(np1:np),     star%Cvg(np1:np),   star%Gam1g(np1:np), &
         star%Gam31g(np1:np), star%Prhog(np1:np), star%PTg(np1:np), &
         star%Cpg(np1:np),    star%Cprhog(np1:np),star%CpTg(np1:np), &
         star%Qg(np1:np),     star%Qrhog(np1:np), star%QTg(np1:np), &
         star%accrad(np1:np), star%dTdTe(np1:np), star%dTdg(np1:np), &
         star%PR(np1:np),     star%gradPR(np1:np), &
         stat=alloc_stat, errmsg=alloc_msg)
    if (alloc_stat /= 0) then
       write(*,'(A,A,A,A)') 'readmadmod: ERROR - atmosphere allocation failed for ', &
            trim(fn), ': ', trim(alloc_msg)
       close(io)
       call DeleteMadmod(star)
       stop 1
    end if

    read(io, iostat=ios, iomsg=iomsg) &
         star%tau, star%rr, star%mm, star%Pg, star%Cvg, star%Gam1g, &
         star%Gam31g, star%Prhog, star%PTg, star%Cpg, star%Cprhog, star%CpTg, &
         star%Qg, star%Qrhog, star%QTg, star%accrad, star%dTdTe, star%dTdg, &
         star%PR, star%gradPR
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - atmosphere data read failed in ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       call DeleteMadmod(star)
       stop 1
    end if

    close(io, iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - close failed for ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       call DeleteMadmod(star)
       stop 1
    end if
    star%compact_adiabatic = .false.
    star%nonadiabatic_available = valid_nonadiabatic_block(star)
    if (require_full .and. .not. star%nonadiabatic_available) then
       write(*,'(A)') 'readmadmod: ERROR - --nonad requires populated non-adiabatic fields.'
       write(*,'(A,A)') 'readmadmod: neutral or incomplete legacy input: ', trim(fn)
       call DeleteMadmod(star)
       stop 1
    end if
  end subroutine readmadmod


  subroutine read_compact_model(io, fn, file_size, star)
    integer, intent(in) :: io
    character(len=*), intent(in) :: fn
    integer(int64), intent(in) :: file_size
    type(tMadmod), intent(inout) :: star

    integer :: ios, alloc_stat, np, nz
    integer(int64) :: expected_size, integer_bytes, real_bytes
    character(len=512) :: iomsg, alloc_msg

    read(io, pos=len(COMPACT_MAGIC)+1, iostat=ios, iomsg=iomsg) &
         star%npi, star%npa, star%np, star%nz, star%step
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - compact header read failed in ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       stop 1
    end if

    np = star%np
    nz = star%nz
    if (np < 1 .or. np > MAX_MODEL_POINTS .or. star%npi /= np .or. &
         star%npa /= 0 .or. nz < 1 .or. nz > np) then
       write(*,'(A,A)') 'readmadmod: ERROR - invalid compact adiabatic dimensions in ', trim(fn)
       write(*,'(A,I0,A,I0,A,I0,A,I0)') 'readmadmod: npi=', star%npi, &
            ' npa=', star%npa, ' np=', np, ' nz=', nz
       close(io)
       stop 1
    end if

    integer_bytes = int(storage_size(star%npi)/8, int64)
    real_bytes = int(storage_size(star%mass)/8, int64)
    expected_size = int(len(COMPACT_MAGIC), int64) + 5_int64*integer_bytes + &
         3_int64*real_bytes + 3_int64*int(nz,int64)*integer_bytes + &
         6_int64*int(np,int64)*real_bytes
    if (file_size /= expected_size) then
       write(*,'(A,A)') 'readmadmod: ERROR - compact model has an unexpected size: ', trim(fn)
       write(*,'(A,I0,A,I0)') 'readmadmod: file bytes=', file_size, &
            ' expected bytes=', expected_size
       close(io)
       stop 1
    end if

    read(io, iostat=ios, iomsg=iomsg) star%mass, star%radius, star%Grav
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - compact globals read failed in ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       stop 1
    end if

    alloc_msg = ''
    allocate(star%zones(3,nz), star%r(np), star%m(np), star%rho(np), &
         star%P(np), star%Gam1(np), star%n2(np), stat=alloc_stat, errmsg=alloc_msg)
    if (alloc_stat /= 0) then
       write(*,'(A,A,A,A)') 'readmadmod: ERROR - compact allocation failed for ', &
            trim(fn), ': ', trim(alloc_msg)
       close(io)
       call DeleteMadmod(star)
       stop 1
    end if

    read(io, iostat=ios, iomsg=iomsg) star%zones, star%r, star%m, star%rho, &
         star%P, star%Gam1, star%n2
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - compact data read failed in ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       close(io)
       call DeleteMadmod(star)
       stop 1
    end if

    close(io, iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       write(*,'(A,A,A,I0,A,A)') 'readmadmod: ERROR - compact close failed for ', &
            trim(fn), ' (iostat=', ios, '): ', trim(iomsg)
       call DeleteMadmod(star)
       stop 1
    end if
    star%compact_adiabatic = .true.
    star%nonadiabatic_available = .false.
  end subroutine read_compact_model


  logical function valid_nonadiabatic_block(star)
    type(tMadmod), intent(in) :: star
    integer :: n

    valid_nonadiabatic_block = .false.
    ! Some historical full files use neutral placeholders in atmosphere-only
    ! records. Validate the populated interior block without rejecting those
    ! otherwise readable legacy models.
    n = star%npi
    if (n < 1) return
    if (.not. allocated(star%T) .or. .not. allocated(star%Cv) .or. &
         .not. allocated(star%Gam31) .or. .not. allocated(star%Prho) .or. &
         .not. allocated(star%PT) .or. .not. allocated(star%Cp) .or. &
         .not. allocated(star%Cprho) .or. .not. allocated(star%CpT) .or. &
         .not. allocated(star%Q) .or. .not. allocated(star%Qrho) .or. &
         .not. allocated(star%QT) .or. .not. allocated(star%kappa) .or. &
         .not. allocated(star%L) .or. .not. allocated(star%gradAd) .or. &
         .not. allocated(star%grad) .or. .not. allocated(star%gradRad) .or. &
         .not. allocated(star%X) .or. .not. allocated(star%Z)) return
    if (.not. ieee_is_finite(star%Lum) .or. star%Lum <= 0d0) return
    if (any(.not. ieee_is_finite(star%T(1:n))) .or. any(star%T(1:n) <= 0d0)) return
    if (any(.not. ieee_is_finite(star%Cv(1:n))) .or. any(star%Cv(1:n) <= 0d0)) return
    if (any(.not. ieee_is_finite(star%Gam31(1:n))) .or. any(star%Gam31(1:n) <= 0d0)) return
    if (any(.not. ieee_is_finite(star%Prho(1:n))) .or. any(star%Prho(1:n) <= 0d0)) return
    if (any(.not. ieee_is_finite(star%PT(1:n))) .or. any(star%PT(1:n) <= 0d0)) return
    if (any(.not. ieee_is_finite(star%Cp(1:n))) .or. any(star%Cp(1:n) <= 0d0)) return
    if (any(.not. ieee_is_finite(star%kappa(1:n))) .or. any(star%kappa(1:n) <= 0d0)) return
    if (any(.not. ieee_is_finite(star%gradAd(1:n))) .or. any(star%gradAd(1:n) <= 0d0)) return
    if (any(.not. ieee_is_finite(star%L(1:n))) .or. &
         any(.not. ieee_is_finite(star%grad(1:n))) .or. &
         any(.not. ieee_is_finite(star%gradRad(1:n))) .or. &
         any(.not. ieee_is_finite(star%X(1:n))) .or. &
         any(.not. ieee_is_finite(star%Z(1:n)))) return
    if (any(.not. ieee_is_finite(star%Cprho(1:n))) .or. &
         any(.not. ieee_is_finite(star%CpT(1:n))) .or. &
         any(.not. ieee_is_finite(star%Q(1:n))) .or. &
         any(.not. ieee_is_finite(star%Qrho(1:n))) .or. &
         any(.not. ieee_is_finite(star%QT(1:n)))) return
    valid_nonadiabatic_block = .true.
  end function valid_nonadiabatic_block


  subroutine DeleteMadmod(star)
    type(tMadmod), intent(inout) :: star

    if (allocated(star%zones))    deallocate(star%zones)
    if (allocated(star%r))        deallocate(star%r)
    if (allocated(star%m))        deallocate(star%m)
    if (allocated(star%rho))      deallocate(star%rho)
    if (allocated(star%T))        deallocate(star%T)
    if (allocated(star%L))        deallocate(star%L)
    if (allocated(star%P))        deallocate(star%P)
    if (allocated(star%Cv))       deallocate(star%Cv)
    if (allocated(star%Gam1))     deallocate(star%Gam1)
    if (allocated(star%Gam31))    deallocate(star%Gam31)
    if (allocated(star%Prho))     deallocate(star%Prho)
    if (allocated(star%PT))       deallocate(star%PT)
    if (allocated(star%Cp))       deallocate(star%Cp)
    if (allocated(star%Cprho))    deallocate(star%Cprho)
    if (allocated(star%CpT))      deallocate(star%CpT)
    if (allocated(star%Q))        deallocate(star%Q)
    if (allocated(star%Qrho))     deallocate(star%Qrho)
    if (allocated(star%QT))       deallocate(star%QT)
    if (allocated(star%kappa))    deallocate(star%kappa)
    if (allocated(star%krho))     deallocate(star%krho)
    if (allocated(star%kT))       deallocate(star%kT)
    if (allocated(star%en))       deallocate(star%en)
    if (allocated(star%X))        deallocate(star%X)
    if (allocated(star%Z))        deallocate(star%Z)
    if (allocated(star%gradP))    deallocate(star%gradP)
    if (allocated(star%gradT))    deallocate(star%gradT)
    if (allocated(star%gradRad))  deallocate(star%gradRad)
    if (allocated(star%gradAd))   deallocate(star%gradAd)
    if (allocated(star%grad))     deallocate(star%grad)
    if (allocated(star%GamC))     deallocate(star%GamC)
    if (allocated(star%AlphaC))   deallocate(star%AlphaC)
    if (allocated(star%dm1))      deallocate(star%dm1)
    if (allocated(star%dm2))      deallocate(star%dm2)
    if (allocated(star%yH1))      deallocate(star%yH1)
    if (allocated(star%yH2))      deallocate(star%yH2)
    if (allocated(star%yHe3))     deallocate(star%yHe3)
    if (allocated(star%yHe4))     deallocate(star%yHe4)
    if (allocated(star%yLi7))     deallocate(star%yLi7)
    if (allocated(star%yC12))     deallocate(star%yC12)
    if (allocated(star%yC13))     deallocate(star%yC13)
    if (allocated(star%yN14))     deallocate(star%yN14)
    if (allocated(star%yN15))     deallocate(star%yN15)
    if (allocated(star%yO16))     deallocate(star%yO16)
    if (allocated(star%yO17))     deallocate(star%yO17)
    if (allocated(star%yO18))     deallocate(star%yO18)
    if (allocated(star%yNe20))    deallocate(star%yNe20)
    if (allocated(star%yOthers))  deallocate(star%yOthers)
    if (allocated(star%ZOthers))  deallocate(star%ZOthers)
    if (allocated(star%ZZOthers)) deallocate(star%ZZOthers)
    if (allocated(star%MOthers))  deallocate(star%MOthers)
    if (allocated(star%eg))       deallocate(star%eg)
    if (allocated(star%mix))      deallocate(star%mix)
    if (allocated(star%n2))       deallocate(star%n2)
    ! Atmosphere arrays: only allocated when npa > 0
    if (allocated(star%tau))      deallocate(star%tau)
    if (allocated(star%rr))       deallocate(star%rr)
    if (allocated(star%mm))       deallocate(star%mm)
    if (allocated(star%Pg))       deallocate(star%Pg)
    if (allocated(star%Cvg))      deallocate(star%Cvg)
    if (allocated(star%Gam1g))    deallocate(star%Gam1g)
    if (allocated(star%Gam31g))   deallocate(star%Gam31g)
    if (allocated(star%Prhog))    deallocate(star%Prhog)
    if (allocated(star%PTg))      deallocate(star%PTg)
    if (allocated(star%Cpg))      deallocate(star%Cpg)
    if (allocated(star%Cprhog))   deallocate(star%Cprhog)
    if (allocated(star%CpTg))     deallocate(star%CpTg)
    if (allocated(star%Qg))       deallocate(star%Qg)
    if (allocated(star%Qrhog))    deallocate(star%Qrhog)
    if (allocated(star%QTg))      deallocate(star%QTg)
    if (allocated(star%accrad))   deallocate(star%accrad)
    if (allocated(star%dTdTe))    deallocate(star%dTdTe)
    if (allocated(star%dTdg))     deallocate(star%dTdg)
    if (allocated(star%PR))       deallocate(star%PR)
    if (allocated(star%gradPR))   deallocate(star%gradPR)

    star%npi = 0
    star%npa = 0
    star%np  = 0
    star%nz  = 0
    star%compact_adiabatic = .false.
    star%nonadiabatic_available = .false.
    star%mass = 0d0
    star%radius = 0d0
    star%Grav = 0d0
  end subroutine DeleteMadmod

end module madmodlib
