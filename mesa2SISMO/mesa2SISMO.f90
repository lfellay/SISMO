program mesa2SISMO
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 300)
  integer, parameter :: max_words = 512
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: clight = 2.99792458d10
  real(dp), parameter :: year = 365.25d0 * 86400d0
  real(dp), parameter :: grav_default = 6.67430d-8
  real(dp), parameter :: arad_default = 7.56573325028001d-15
  real(dp), parameter :: sigma_default = arad_default * clight / 4d0
  real(dp), parameter :: msun_default = 1.98840987069805d33
  real(dp), parameter :: rsun_default = 6.957d10
  real(dp), parameter :: lsun_default = 3.828d33
  real(dp), parameter :: tiny = 1d-99
  real(dp), parameter :: atmosphere_base_tau = 2d0/3d0

  type mesa_profile
     integer :: nhead = 0
     integer :: ncol = 0
     integer :: nrow = 0
     character(len=64), allocatable :: hname(:)
     character(len=64), allocatable :: cname(:)
     character(len=128), allocatable :: htoken(:)
     real(dp), allocatable :: hval(:)
     logical, allocatable :: hnum(:)
     real(dp), allocatable :: data(:,:)
  end type mesa_profile

  type mad_model
     integer :: npi = 0
     integer :: npa = 0
     integer :: np = 0
     integer :: nz = 0
     integer :: step = 0
     real(dp) :: mass = 0d0
     real(dp) :: radius = 0d0
     real(dp) :: teff = 0d0
     real(dp) :: lum = 0d0
     real(dp) :: age = 0d0
     real(dp) :: logg = 0d0
     real(dp) :: lgOverL = 0d0
     real(dp) :: alphaConv = 0d0
     real(dp) :: alphaOver = 0d0
     real(dp) :: X0 = 0d0
     real(dp) :: Z0 = 0d0
     real(dp) :: diff = 0d0
     real(dp) :: grav = grav_default
     real(dp) :: arad = arad_default
     real(dp) :: sigma = sigma_default
     integer, allocatable :: zones(:,:)
     real(dp), allocatable :: r(:), m(:), rho(:), T(:), L(:), P(:)
     real(dp), allocatable :: Cv(:), Gam1(:), Gam31(:), Prho(:), PT(:)
     real(dp), allocatable :: Cp(:), Cprho(:), CpT(:), Q(:), Qrho(:), QT(:)
     real(dp), allocatable :: kappa(:), krho(:), kT(:), en(:)
     real(dp), allocatable :: X(:), Z(:), gradP(:), gradT(:), gradRad(:)
     real(dp), allocatable :: gradAd(:), grad(:), GamC(:), AlphaC(:)
     real(dp), allocatable :: dm1(:), dm2(:), yH1(:), yH2(:), yHe3(:), yHe4(:)
     real(dp), allocatable :: yLi7(:), yC12(:), yC13(:), yN14(:), yN15(:)
     real(dp), allocatable :: yO16(:), yO17(:), yO18(:), yNe20(:)
     real(dp), allocatable :: yOthers(:), ZOthers(:), ZZOthers(:), MOthers(:)
     real(dp), allocatable :: eg(:), mix(:), n2(:)
     real(dp), allocatable :: tau(:), rr(:), mm(:), Pg(:), Cvg(:), Gam1g(:)
     real(dp), allocatable :: Gam31g(:), Prhog(:), PTg(:), Cpg(:)
     real(dp), allocatable :: Cprhog(:), CpTg(:), Qg(:), Qrhog(:), QTg(:)
     real(dp), allocatable :: accrad(:), dTdTe(:), dTdg(:), PR(:), gradPR(:)
  end type mad_model

  type mesa_columns
     integer :: mass_grams = 0, mass_msun = 0, radius_cm = 0
     integer :: temperature = 0, rho = 0, pressure = 0, pgas = 0, prad = 0
     integer :: energy = 0, grada = 0, gamma1 = 0, chirho = 0, chit = 0
     integer :: cp = 0, cv = 0, dedrho = 0
     integer :: luminosity = 0, lum_erg_s = 0, lum_rad = 0, lum_conv = 0
     integer :: conv_l_div_l = 0, opacity = 0, krho = 0, kt = 0
     integer :: gradt = 0, gradr = 0, eps_nuc = 0, tau = 0, brunt_n2 = 0
     integer :: x = 0, y = 0, z = 0, h1 = 0, he3 = 0, he4 = 0
     integer :: c12 = 0, c13 = 0, n14 = 0, n15 = 0, o16 = 0, o17 = 0
     integer :: o18 = 0, ne20 = 0, cprho = 0, cpt = 0, qrho = 0, qt = 0
     integer :: alpha_mlt = 0
  end type mesa_columns

  character(len=512) :: infile, outfile
  logical :: adiabatic_only
  type(mesa_profile) :: profile
  type(mad_model) :: model

  call parse_args(infile, outfile, adiabatic_only)
  call ensure_distinct_files(infile, outfile)
  call read_mesa_profile(infile, profile)
  call convert_to_mad(profile, adiabatic_only, model)
  call write_madmod(outfile, model, adiabatic_only)

  write(*,'(a)') 'mesa2SISMO: wrote '//trim(outfile)
  write(*,'(a,i0,a,i0,a,i0,a,i0)') 'mesa2SISMO: points npi=', model%npi, &
       ' npa=', model%npa, ' np=', model%np, ' zones=', model%nz
  if (model%npa > 0) then
     write(*,'(a)') 'mesa2SISMO: atmosphere source=MESA profile rows above tau=2/3'
  else
     write(*,'(a)') 'mesa2SISMO: atmosphere source=none in MESA profile; converter did not add one'
  end if
  if (adiabatic_only) then
     write(*,'(a)') 'mesa2SISMO: mode=adiabatic (default compact model)'
  else
     write(*,'(a)') 'mesa2SISMO: mode=nonadiabatic (--nonad full legacy model)'
  end if
  write(*,'(a,1pe11.4,a,1pe11.4,a,1pe11.4)') 'mesa2SISMO: M=', model%mass/msun_default, &
       ' R=', model%radius/rsun_default, ' L=', model%lum/lsun_default

contains

  subroutine parse_args(infile, outfile, adiabatic_only)
    character(len=*), intent(out) :: infile, outfile
    logical, intent(out) :: adiabatic_only

    integer :: narg, i
    character(len=512) :: arg
    logical :: have_out, mode_selected

    narg = command_argument_count()
    if (narg < 1) then
       call usage()
       stop 1
    end if

    call get_command_argument(1, infile)
    if (trim(infile) == '--help' .or. trim(infile) == '-h') then
       call usage()
       stop
    end if

    outfile = default_output_name(infile)
    adiabatic_only = .true.
    have_out = .false.
    mode_selected = .false.

    i = 2
    do while (i <= narg)
       call get_command_argument(i, arg)
       select case (trim(arg))
       case ('--help', '-h')
          call usage()
          stop
       case ('--adiabatic', '--ad')
          if (mode_selected .and. .not. adiabatic_only) then
             call fatal('conflicting physics options: choose either adiabatic mode or --nonad.')
          end if
          adiabatic_only = .true.
          mode_selected = .true.
       case ('--nonad')
          if (mode_selected .and. adiabatic_only) then
             call fatal('conflicting physics options: choose either adiabatic mode or --nonad.')
          end if
          adiabatic_only = .false.
          mode_selected = .true.
       case default
          if (len_trim(arg) >= 2 .and. arg(1:2) == '--') then
             call fatal('unknown option '//trim(arg))
          end if
          if (have_out) call fatal('only one output file may be supplied.')
          outfile = arg
          have_out = .true.
       end select
       i = i + 1
    end do
  end subroutine parse_args

  subroutine usage()
    write(*,'(a)') 'Usage: mesa2SISMO profile.data [output.madmod] [--nonad]'
    write(*,'(a)') 'Default: write a compact adiabatic model.'
    write(*,'(a)') '--nonad: require and write the full legacy non-adiabatic model.'
    write(*,'(a)') '--adiabatic and --ad remain accepted for compatibility.'
    write(*,'(a)') 'The converter never adds atmosphere points; any atmosphere must already be present in the MESA input profile.'
  end subroutine usage

  function default_output_name(infile) result(outfile)
    character(len=*), intent(in) :: infile
    character(len=512) :: outfile
    integer :: n

    outfile = trim(infile)
    n = len_trim(outfile)
    if (n >= 5 .and. outfile(n-4:n) == '.data') then
       outfile = outfile(:n-5)//'.madmod'
    else
       outfile = trim(outfile)//'.madmod'
    end if
  end function default_output_name

  subroutine ensure_distinct_files(input_name, output_name)
    character(len=*), intent(in) :: input_name, output_name

    integer :: cmdstat, exitstat
    logical :: output_exists
    character(len=512) :: cmdmsg
    character(len=:), allocatable :: command

    if (trim(input_name) == trim(output_name)) then
       call fatal('input and output must be different files: '//trim(input_name))
    end if

    inquire(file=trim(output_name), exist=output_exists)
    if (.not. output_exists) return

    ! /bin/test -ef compares the underlying file identities, so this also
    ! catches relative/absolute aliases, symbolic links, and hard links.
    ! shell_quote prevents either user-supplied path from being interpreted.
    command = '/bin/test '//shell_quote(trim(input_name))//' -ef '// &
         shell_quote(trim(output_name))
    cmdmsg = ''
    call execute_command_line(command, wait=.true., exitstat=exitstat, &
         cmdstat=cmdstat, cmdmsg=cmdmsg)
    if (cmdstat /= 0) then
       call fatal('could not verify that input and output differ: '//trim(cmdmsg))
    end if
    if (exitstat == 0) then
       call fatal('input and output refer to the same file: '//trim(input_name))
    else if (exitstat /= 1) then
       call fatal('could not compare input and output paths safely.')
    end if
  end subroutine ensure_distinct_files

  function shell_quote(text) result(quoted)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: quoted

    integer :: i

    quoted = achar(39)
    do i = 1, len_trim(text)
       if (text(i:i) == achar(39)) then
          ! End the single-quoted string, quote one apostrophe, and reopen it.
          quoted = quoted//achar(39)//achar(34)//achar(39)//achar(34)//achar(39)
       else
          quoted = quoted//text(i:i)
       end if
    end do
    quoted = quoted//achar(39)
  end function shell_quote

  subroutine read_mesa_profile(filename, p)
    character(len=*), intent(in) :: filename
    type(mesa_profile), intent(out) :: p

    integer :: io, ios, i, j, ntok, nrow_header
    character(len=4096) :: line
    character(len=128) :: tokens(max_words)

    open(newunit=io, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,'(a,a,a,i0)') 'mesa2SISMO: cannot open ', trim(filename), ' iostat=', ios
       stop 1
    end if

    call read_required_line(io, line, 'header numbers')
    call read_required_line(io, line, 'header names')
    call split_words(line, tokens, ntok)
    if (ntok < 1) call fatal('empty MESA header-name line.')
    p%nhead = ntok
    allocate(p%hname(p%nhead), p%htoken(p%nhead), p%hval(p%nhead), p%hnum(p%nhead))
    do i = 1, p%nhead
       p%hname(i) = trim(tokens(i))
    end do

    call read_required_line(io, line, 'header values')
    call split_words(line, tokens, ntok)
    if (ntok < p%nhead) call fatal('MESA header has fewer values than names.')
    do i = 1, p%nhead
       p%htoken(i) = trim(tokens(i))
       read(tokens(i),*,iostat=ios) p%hval(i)
       p%hnum(i) = (ios == 0)
       if (.not. p%hnum(i)) p%hval(i) = 0d0
    end do

    call read_required_line(io, line, 'blank separator')
    call read_required_line(io, line, 'profile column numbers')
    call read_required_line(io, line, 'profile column names')
    call split_words(line, tokens, ntok)
    if (ntok < 1) call fatal('empty MESA profile-column line.')
    p%ncol = ntok
    allocate(p%cname(p%ncol))
    do i = 1, p%ncol
       p%cname(i) = trim(tokens(i))
    end do

    nrow_header = int(header_value(p, 'num_zones', -1d0))
    if (nrow_header <= 0) call fatal('MESA header does not contain a valid num_zones.')
    p%nrow = nrow_header
    allocate(p%data(p%ncol, p%nrow))

    do i = 1, p%nrow
       read(io,*,iostat=ios) (p%data(j,i), j=1,p%ncol)
       if (ios /= 0) then
          write(*,'(a,i0,a,i0)') 'mesa2SISMO: error reading profile row ', i, ' iostat=', ios
          stop 1
       end if
    end do

    close(io)
  end subroutine read_mesa_profile

  subroutine read_required_line(io, line, label)
    integer, intent(in) :: io
    character(len=*), intent(out) :: line
    character(len=*), intent(in) :: label
    integer :: ios

    read(io,'(a)',iostat=ios) line
    if (ios /= 0) call fatal('truncated MESA profile while reading '//trim(label)//'.')
  end subroutine read_required_line

  subroutine convert_to_mad(p, adiabatic_only, model)
    type(mesa_profile), intent(in) :: p
    logical, intent(in) :: adiabatic_only
    type(mad_model), intent(out) :: model

    type(mesa_columns) :: c
    integer :: n_mesa, np, i, row, npa, npi, base_row, native_npa
    real(dp) :: msun, rsun, lsun, eps_max, eps_cut, radius_scale, base_radius
    real(dp) :: pgas, prad, beta, conv_metric_i
    real(dp), allocatable :: conv_metric(:)

    call map_columns(p, c, adiabatic_only)
    call validate_positive_column(p, c%rho)
    call validate_positive_column(p, c%pressure)
    call validate_positive_column(p, c%gamma1)
    if (c%temperature > 0) call validate_positive_column(p, c%temperature)

    n_mesa = p%nrow
    if (n_mesa < 16) call fatal('MESA profile is too small for MAD/intSISMO.')

    msun = header_value(p, 'msun', msun_default)
    rsun = header_value(p, 'rsun', rsun_default)
    lsun = header_value(p, 'lsun', lsun_default)

    if (adiabatic_only) then
       base_row = 1
       native_npa = 0
    else
       base_row = atmosphere_base_row(p, c, atmosphere_base_tau)
       native_npa = max(0, base_row - 1)
    end if

    model%step = int(header_value(p, 'model_number', 0d0))
    model%mass = header_value(p, 'star_mass', 0d0) * msun
    if (model%mass <= 0d0) model%mass = col_surface(p, c%mass_grams)
    if (model%mass <= 0d0 .and. c%mass_msun > 0) model%mass = col_surface(p, c%mass_msun) * msun

    base_radius = optional_col(p, c%radius_cm, base_row, 0d0)
    if (base_radius <= 0d0) base_radius = header_value(p, 'photosphere_r', 0d0) * rsun
    if (base_radius <= 0d0) base_radius = col_surface(p, c%radius_cm)
    model%radius = base_radius

    model%lum = get_luminosity(p, c, base_row, lsun)
    if (model%lum <= 0d0) model%lum = header_value(p, 'photosphere_L', 0d0) * lsun
    if (model%lum <= 0d0) model%lum = get_luminosity(p, c, 1, lsun)

    model%teff = header_value(p, 'Teff', 0d0)
    if (model%teff <= 0d0 .and. model%radius > 0d0 .and. model%lum > 0d0) then
       model%teff = (model%lum / (4d0*pi*model%radius**2*sigma_default))**0.25d0
    end if
    model%age = header_value(p, 'star_age', 0d0) * year
    model%Z0 = header_value(p, 'initial_z', 0d0)
    model%alphaConv = col_surface(p, c%alpha_mlt)
    if (model%alphaConv <= 0d0) model%alphaConv = 2d0
    model%grav = grav_default
    model%arad = arad_default
    model%sigma = sigma_default
    model%alphaOver = 0d0
    model%diff = 0d0

    if (model%mass <= 0d0 .or. model%radius <= 0d0 .or. &
         ((.not. adiabatic_only) .and. model%lum <= 0d0)) then
       call fatal('invalid global mass/radius/luminosity in MESA profile.')
    end if
    model%logg = log10(model%grav * model%mass / model%radius**2)
    model%lgOverL = 0d0

    if (adiabatic_only) then
       npa = 0
       npi = n_mesa
       np = n_mesa
    else if (native_npa > 0) then
       npa = native_npa
       npi = n_mesa - npa
       np = n_mesa
    else
       npa = 0
       npi = n_mesa
       np = n_mesa
    end if
    model%np = np
    model%npi = npi
    model%npa = npa

    if (base_radius <= 0d0) base_radius = model%radius
    radius_scale = model%radius / base_radius

    call allocate_model_arrays(model, np)
    allocate(conv_metric(np))
    conv_metric = 0d0

    eps_max = 0d0
    if (c%eps_nuc > 0) eps_max = maxval(max(0d0, p%data(c%eps_nuc, :)))
    eps_cut = max(1d-30, eps_max * 1d-12)

    do i = 1, n_mesa
       row = n_mesa - i + 1

       model%r(i) = max(0d0, p%data(c%radius_cm,row) * radius_scale)
       model%m(i) = get_mass(p, c, row, msun)
       if (i == 1) then
          model%r(i) = 0d0
          model%m(i) = 0d0
       end if
       model%T(i) = max(tiny, optional_col(p, c%temperature, row, tiny))
       model%rho(i) = max(tiny, p%data(c%rho,row))
       model%P(i) = max(tiny, p%data(c%pressure,row))
       model%L(i) = get_luminosity(p, c, row, lsun)
       model%Cv(i) = max(tiny, optional_col(p, c%cv, row, 0d0))
       model%Gam1(i) = max(tiny, p%data(c%gamma1,row))
       model%gradAd(i) = max(tiny, optional_col(p, c%grada, row, 0.4d0))
       model%Gam31(i) = model%Gam1(i) * model%gradAd(i)
       model%Prho(i) = max(tiny, optional_col(p, c%chirho, row, 1d0))
       model%PT(i) = max(tiny, optional_col(p, c%chit, row, 0d0))
       model%Cp(i) = max(tiny, optional_col(p, c%cp, row, 0d0))
       model%Q(i) = model%PT(i) / model%Prho(i)
       model%Cprho(i) = optional_col(p, c%cprho, row, 0d0)
       model%CpT(i) = optional_col(p, c%cpt, row, 0d0)
       model%Qrho(i) = optional_col(p, c%qrho, row, 0d0)
       model%QT(i) = optional_col(p, c%qt, row, 0d0)
       model%kappa(i) = max(0d0, optional_col(p, c%opacity, row, 0d0))
       model%krho(i) = optional_col(p, c%krho, row, 0d0)
       model%kT(i) = optional_col(p, c%kt, row, 0d0)
       model%en(i) = optional_col(p, c%eps_nuc, row, 0d0)
       if (model%en(i) <= eps_cut) model%en(i) = 0d0
       model%eg(i) = 0d0
       model%grad(i) = max(0d0, optional_col(p, c%gradt, row, 0d0))
       model%gradRad(i) = max(0d0, optional_col(p, c%gradr, row, model%grad(i)))
       model%n2(i) = optional_col(p, c%brunt_n2, row, 0d0)
       model%AlphaC(i) = model%alphaConv

       model%X(i) = get_composition_x(p, c, row)
       model%Z(i) = get_composition_z(p, c, row)
       model%yH1(i) = optional_col(p, c%h1, row, model%X(i))
       model%yHe3(i) = optional_col(p, c%he3, row, 0d0)
       model%yHe4(i) = optional_col(p, c%he4, row, max(0d0, 1d0-model%X(i)-model%Z(i)))
       model%yC12(i) = optional_col(p, c%c12, row, 0d0)
       model%yC13(i) = optional_col(p, c%c13, row, 0d0)
       model%yN14(i) = optional_col(p, c%n14, row, 0d0)
       model%yN15(i) = optional_col(p, c%n15, row, 0d0)
       model%yO16(i) = optional_col(p, c%o16, row, 0d0)
       model%yO17(i) = optional_col(p, c%o17, row, 0d0)
       model%yO18(i) = optional_col(p, c%o18, row, 0d0)
       model%yNe20(i) = optional_col(p, c%ne20, row, 0d0)
       model%yOthers(i) = max(0d0, model%Z(i) - model%yC12(i) - model%yC13(i) &
            - model%yN14(i) - model%yN15(i) - model%yO16(i) - model%yO17(i) &
            - model%yO18(i) - model%yNe20(i))
       model%ZOthers(i) = model%yOthers(i)
       model%ZZOthers(i) = model%yOthers(i)
       model%MOthers(i) = model%yOthers(i)

       model%tau(i) = max(tiny, optional_col(p, c%tau, row, 1d50))
       model%rr(i) = model%r(i)
       model%mm(i) = model%m(i)
       pgas = optional_col(p, c%pgas, row, max(tiny, model%P(i) - arad_default*model%T(i)**4/3d0))
       prad = optional_col(p, c%prad, row, max(0d0, arad_default*model%T(i)**4/3d0))
       pgas = max(tiny, pgas)
       beta = min(1d0, max(tiny, pgas/model%P(i)))
       model%Pg(i) = pgas
       model%Cvg(i) = max(tiny, model%Cv(i) - 4d0*arad_default*model%T(i)**3/model%rho(i))
       model%Gam1g(i) = model%Gam1(i)
       model%Gam31g(i) = model%Gam31(i)
       model%Prhog(i) = model%Prho(i) / beta
       model%PTg(i) = (model%PT(i) - 4d0*(1d0-beta)) / beta
       model%Cpg(i) = model%Cp(i)
       model%Cprhog(i) = model%Cprho(i)
       model%CpTg(i) = model%CpT(i)
       model%Qg(i) = model%PTg(i) / max(tiny, model%Prhog(i))
       model%Qrhog(i) = model%Qrho(i)
       model%QTg(i) = model%QT(i)
       model%PR(i) = prad
       model%gradPR(i) = 0d0
       model%dTdTe(i) = 1d0
       model%dTdg(i) = 0d0

       if (i > 1 .and. model%r(i) > tiny) then
          model%gradP(i) = -model%grav * model%m(i) * model%rho(i) / model%r(i)**2
          model%gradT(i) = model%grad(i) * model%T(i) * model%gradP(i) / model%P(i)
          model%accrad(i) = (4d0/3d0) * model%arad * model%T(i)**4 &
               * model%grav * model%m(i) * max(0d0, model%grad(i)) &
               / (model%r(i)**2 * model%P(i))
       else
          model%gradP(i) = 0d0
          model%gradT(i) = 0d0
          model%accrad(i) = 0d0
       end if

       conv_metric_i = max(model%gradRad(i), model%grad(i)) - model%gradAd(i)
       if (c%conv_l_div_l > 0) conv_metric_i = max(conv_metric_i, p%data(c%conv_l_div_l,row))
       conv_metric(i) = conv_metric_i
    end do

    model%X0 = model%X(model%npi)
    if (model%Z0 <= 0d0) model%Z0 = model%Z(model%npi)

    call validate_monotonic_model(model)

    call compute_gamc(model)
    call build_zones(model, conv_metric)

    if (adiabatic_only) call neutralize_nonadiabatic_fields(model)

    deallocate(conv_metric)
  end subroutine convert_to_mad

  subroutine map_columns(p, c, adiabatic_only)
    type(mesa_profile), intent(in) :: p
    type(mesa_columns), intent(out) :: c
    logical, intent(in) :: adiabatic_only

    c%mass_grams = find_col(p, 'mass_grams')
    c%mass_msun = find_col(p, 'mass')
    if (c%mass_grams <= 0 .and. c%mass_msun <= 0) call fatal('missing mass_grams/mass column.')
    c%radius_cm = require_col(p, 'radius_cm')
    if (adiabatic_only) then
       c%temperature = find_col(p, 'temperature')
    else
       c%temperature = require_col(p, 'temperature')
    end if
    c%rho = require_col(p, 'rho')
    c%pressure = require_col(p, 'pressure')
    c%pgas = find_col(p, 'pgas')
    c%prad = find_col(p, 'prad')
    c%energy = find_col(p, 'energy')
    if (adiabatic_only) then
       c%grada = find_col(p, 'grada')
    else
       c%grada = require_col(p, 'grada')
    end if
    c%gamma1 = require_col(p, 'gamma1')
    if (adiabatic_only) then
       c%chirho = find_col(p, 'chiRho')
       c%chit = find_col(p, 'chiT')
       c%cp = find_col(p, 'cp')
       c%cv = find_col(p, 'cv')
    else
       c%chirho = require_col(p, 'chiRho')
       c%chit = require_col(p, 'chiT')
       c%cp = require_col(p, 'cp')
       c%cv = require_col(p, 'cv')
    end if
    c%dedrho = find_col(p, 'dE_dRho')
    c%luminosity = find_col(p, 'luminosity')
    c%lum_erg_s = find_col(p, 'lum_erg_s')
    if ((.not. adiabatic_only) .and. c%luminosity <= 0 .and. c%lum_erg_s <= 0) then
       call fatal('missing luminosity/lum_erg_s column.')
    end if
    c%lum_rad = find_col(p, 'lum_rad')
    c%lum_conv = find_col(p, 'lum_conv')
    c%conv_l_div_l = find_col(p, 'conv_L_div_L')
    if (adiabatic_only) then
       c%opacity = find_col(p, 'opacity')
    else
       c%opacity = require_col(p, 'opacity')
    end if
    c%krho = find_col(p, 'dkap_dlnrho_face')
    c%kt = find_col(p, 'dkap_dlnT_face')
    if (adiabatic_only) then
       c%gradt = find_col(p, 'gradT')
       c%gradr = find_col(p, 'gradr')
    else
       c%gradt = require_col(p, 'gradT')
       c%gradr = require_col(p, 'gradr')
    end if
    c%eps_nuc = find_col(p, 'eps_nuc')
    c%tau = find_col(p, 'tau')
    c%brunt_n2 = find_col(p, 'brunt_N2')
    c%x = find_first_col(p, 'x', 'x_mass_fraction_H')
    c%y = find_first_col(p, 'y', 'y_mass_fraction_He')
    c%z = find_first_col(p, 'z', 'z_mass_fraction_metals')
    if ((.not. adiabatic_only) .and. c%x <= 0) call fatal('missing x/x_mass_fraction_H column.')
    if ((.not. adiabatic_only) .and. c%z <= 0) call fatal('missing z/z_mass_fraction_metals column.')
    c%h1 = find_col(p, 'h1')
    c%he3 = find_col(p, 'he3')
    c%he4 = find_col(p, 'he4')
    c%c12 = find_col(p, 'c12')
    c%c13 = find_col(p, 'c13')
    c%n14 = find_col(p, 'n14')
    c%n15 = find_col(p, 'n15')
    c%o16 = find_col(p, 'o16')
    c%o17 = find_col(p, 'o17')
    c%o18 = find_col(p, 'o18')
    c%ne20 = find_col(p, 'ne20')
    if (adiabatic_only) then
       c%cprho = find_col(p, 'Cprho')
       c%cpt = find_col(p, 'CpT')
       c%qrho = find_col(p, 'Qrho')
       c%qt = find_col(p, 'QT')
    else
       c%cprho = require_col(p, 'Cprho')
       c%cpt = require_col(p, 'CpT')
       c%qrho = require_col(p, 'Qrho')
       c%qt = require_col(p, 'QT')
    end if
    c%alpha_mlt = find_col(p, 'alpha_mlt')
  end subroutine map_columns

  subroutine validate_positive_column(p, column)
    type(mesa_profile), intent(in) :: p
    integer, intent(in) :: column

    integer :: row
    real(dp) :: value

    if (column < 1 .or. column > p%ncol) then
       call fatal('internal error while validating a mechanical profile column.')
    end if

    do row = 1, p%nrow
       value = p%data(column,row)
       if (.not. ieee_is_finite(value)) then
          write(*,'(a,i0,a,i0,a,a,a,1pe16.8)') &
               'mesa2SISMO: invalid non-finite mechanical value at profile row ', &
               row, ', column ', column, ' (', trim(p%cname(column)), '): ', value
          stop 1
       end if
       if (value <= 0d0) then
          write(*,'(a,i0,a,i0,a,a,a,1pe16.8)') &
               'mesa2SISMO: non-positive mechanical value at profile row ', &
               row, ', column ', column, ' (', trim(p%cname(column)), '): ', value
          stop 1
       end if
    end do
  end subroutine validate_positive_column

  function get_mass(p, c, row, msun) result(val)
    type(mesa_profile), intent(in) :: p
    type(mesa_columns), intent(in) :: c
    integer, intent(in) :: row
    real(dp), intent(in) :: msun
    real(dp) :: val

    if (c%mass_grams > 0) then
       val = p%data(c%mass_grams,row)
    else
       val = p%data(c%mass_msun,row) * msun
    end if
  end function get_mass

  function get_luminosity(p, c, row, lsun) result(val)
    type(mesa_profile), intent(in) :: p
    type(mesa_columns), intent(in) :: c
    integer, intent(in) :: row
    real(dp), intent(in) :: lsun
    real(dp) :: val

    if (c%lum_erg_s > 0) then
       val = p%data(c%lum_erg_s,row)
    else if (c%luminosity > 0) then
       val = p%data(c%luminosity,row) * lsun
    else
       val = 0d0
    end if
  end function get_luminosity

  function get_composition_x(p, c, row) result(val)
    type(mesa_profile), intent(in) :: p
    type(mesa_columns), intent(in) :: c
    integer, intent(in) :: row
    real(dp) :: val

    if (c%x > 0) then
       val = max(0d0, p%data(c%x,row))
    else
       val = 0d0
    end if
  end function get_composition_x

  function get_composition_z(p, c, row) result(val)
    type(mesa_profile), intent(in) :: p
    type(mesa_columns), intent(in) :: c
    integer, intent(in) :: row
    real(dp) :: val

    if (c%z > 0) then
       val = max(0d0, p%data(c%z,row))
    else
       val = 0d0
    end if
  end function get_composition_z

  function atmosphere_base_row(p, c, tau_base) result(row_base)
    type(mesa_profile), intent(in) :: p
    type(mesa_columns), intent(in) :: c
    real(dp), intent(in) :: tau_base
    integer :: row_base

    integer :: row
    real(dp) :: tau_i

    row_base = 1
    if (c%tau <= 0) return

    do row = 1, p%nrow
       tau_i = p%data(c%tau,row)
       if (.not. is_good(tau_i)) cycle
       if (tau_i >= tau_base) then
          row_base = row
          return
       end if
    end do
  end function atmosphere_base_row

  subroutine validate_monotonic_model(model)
    type(mad_model), intent(in) :: model
    integer :: i

    do i = 2, model%np
       if (model%r(i) <= model%r(i-1)) then
          write(*,'(a,i0,2(1pe16.8))') 'mesa2SISMO: non-increasing radius at point ', &
               i, model%r(i-1), model%r(i)
          stop 1
       end if
       if (model%m(i) < model%m(i-1)) then
          write(*,'(a,i0,2(1pe16.8))') 'mesa2SISMO: decreasing mass at point ', &
               i, model%m(i-1), model%m(i)
          stop 1
       end if
    end do
  end subroutine validate_monotonic_model

  subroutine neutralize_nonadiabatic_fields(model)
    type(mad_model), intent(inout) :: model

    model%npi = model%np
    model%npa = 0

    ! The fixed MADMOD record still contains every historical field.  In
    ! adiabatic mode only geometry, P, rho, Gamma1, T, L, composition, globals,
    ! and zone boundaries are meaningful.
    model%Cv = 0d0
    model%Gam31 = 0d0
    model%Prho = 0d0
    model%PT = 0d0
    model%Cp = 0d0
    model%Cprho = 0d0
    model%CpT = 0d0
    model%Q = 0d0
    model%Qrho = 0d0
    model%QT = 0d0
    model%kappa = 0d0
    model%krho = 0d0
    model%kT = 0d0
    model%en = 0d0
    model%gradP = 0d0
    model%gradT = 0d0
    model%gradRad = 0d0
    model%gradAd = 0d0
    model%grad = 0d0
    model%GamC = 0d0
    model%AlphaC = 0d0
    model%dm1 = 0d0
    model%dm2 = 0d0
    model%eg = 0d0
    model%mix = 0d0
    ! Brunt-Vaisala frequency is an adiabatic grid-distribution aid and is
    ! retained when MESA supplies it. It defaults safely to zero otherwise.

    model%tau = 0d0
    model%rr = 0d0
    model%mm = 0d0
    model%Pg = 0d0
    model%Cvg = 0d0
    model%Gam1g = 0d0
    model%Gam31g = 0d0
    model%Prhog = 0d0
    model%PTg = 0d0
    model%Cpg = 0d0
    model%Cprhog = 0d0
    model%CpTg = 0d0
    model%Qg = 0d0
    model%Qrhog = 0d0
    model%QTg = 0d0
    model%accrad = 0d0
    model%dTdTe = 0d0
    model%dTdg = 0d0
    model%PR = 0d0
    model%gradPR = 0d0
  end subroutine neutralize_nonadiabatic_fields

  subroutine allocate_model_arrays(model, np)
    type(mad_model), intent(inout) :: model
    integer, intent(in) :: np

    allocate(model%r(np), model%m(np), model%rho(np), model%T(np), model%L(np), model%P(np))
    allocate(model%Cv(np), model%Gam1(np), model%Gam31(np), model%Prho(np), model%PT(np))
    allocate(model%Cp(np), model%Cprho(np), model%CpT(np), model%Q(np), model%Qrho(np), model%QT(np))
    allocate(model%kappa(np), model%krho(np), model%kT(np), model%en(np))
    allocate(model%X(np), model%Z(np), model%gradP(np), model%gradT(np), model%gradRad(np))
    allocate(model%gradAd(np), model%grad(np), model%GamC(np), model%AlphaC(np))
    allocate(model%dm1(np), model%dm2(np), model%yH1(np), model%yH2(np), model%yHe3(np), model%yHe4(np))
    allocate(model%yLi7(np), model%yC12(np), model%yC13(np), model%yN14(np), model%yN15(np))
    allocate(model%yO16(np), model%yO17(np), model%yO18(np), model%yNe20(np))
    allocate(model%yOthers(np), model%ZOthers(np), model%ZZOthers(np), model%MOthers(np))
    allocate(model%eg(np), model%mix(np), model%n2(np))
    allocate(model%tau(np), model%rr(np), model%mm(np), model%Pg(np), model%Cvg(np), model%Gam1g(np))
    allocate(model%Gam31g(np), model%Prhog(np), model%PTg(np), model%Cpg(np))
    allocate(model%Cprhog(np), model%CpTg(np), model%Qg(np), model%Qrhog(np), model%QTg(np))
    allocate(model%accrad(np), model%dTdTe(np), model%dTdg(np), model%PR(np), model%gradPR(np))

    model%r = 0d0; model%m = 0d0; model%rho = 0d0; model%T = 0d0; model%L = 0d0; model%P = 0d0
    model%Cv = 0d0; model%Gam1 = 0d0; model%Gam31 = 0d0; model%Prho = 0d0; model%PT = 0d0
    model%Cp = 0d0; model%Cprho = 0d0; model%CpT = 0d0; model%Q = 0d0; model%Qrho = 0d0; model%QT = 0d0
    model%kappa = 0d0; model%krho = 0d0; model%kT = 0d0; model%en = 0d0
    model%X = 0d0; model%Z = 0d0; model%gradP = 0d0; model%gradT = 0d0; model%gradRad = 0d0
    model%gradAd = 0d0; model%grad = 0d0; model%GamC = 0d0; model%AlphaC = 0d0
    model%dm1 = 0d0; model%dm2 = 0d0; model%yH1 = 0d0; model%yH2 = 0d0; model%yHe3 = 0d0
    model%yHe4 = 0d0; model%yLi7 = 0d0; model%yC12 = 0d0; model%yC13 = 0d0
    model%yN14 = 0d0; model%yN15 = 0d0; model%yO16 = 0d0; model%yO17 = 0d0
    model%yO18 = 0d0; model%yNe20 = 0d0; model%yOthers = 0d0
    model%ZOthers = 0d0; model%ZZOthers = 0d0; model%MOthers = 0d0
    model%eg = 0d0; model%mix = 0d0; model%n2 = 0d0
    model%tau = 1d50; model%rr = 0d0; model%mm = 0d0; model%Pg = 0d0; model%Cvg = 0d0
    model%Gam1g = 0d0; model%Gam31g = 0d0; model%Prhog = 0d0; model%PTg = 0d0
    model%Cpg = 0d0; model%Cprhog = 0d0; model%CpTg = 0d0; model%Qg = 0d0
    model%Qrhog = 0d0; model%QTg = 0d0; model%accrad = 0d0
    model%dTdTe = 0d0; model%dTdg = 0d0; model%PR = 0d0; model%gradPR = 0d0
  end subroutine allocate_model_arrays

  subroutine compute_gamc(model)
    type(mad_model), intent(inout) :: model

    integer :: i
    real(dp) :: aa, ti, prev, denom

    prev = 0d0
    model%GamC = 0d0
    do i = model%np, 1, -1
       if (i <= 1 .or. model%m(i) <= tiny .or. model%r(i) <= tiny) cycle
       if (model%gradRad(i) <= model%gradAd(i)) cycle

       denom = model%grav * model%m(i) * 48d0 * model%sigma * model%T(i)**3
       if (denom <= tiny) cycle
       aa = model%Q(i) * model%P(i) / (2d0*model%rho(i)) * &
            (model%Cp(i)*model%kappa(i)*model%rho(i)*model%alphaConv**2 &
            * model%P(i)*model%r(i)**2 / denom)**2
       ti = aa * (model%gradRad(i) - model%gradAd(i))
       if (ti > 0d0) then
          model%GamC(i) = gammac_root(ti, prev)
          prev = model%GamC(i)
       end if
    end do
  end subroutine compute_gamc

  function gammac_root(ti, g0) result(g)
    real(dp), intent(in) :: ti, g0
    real(dp) :: g, f, df, dg
    integer :: iter

    if (ti <= 0d0) then
       g = 0d0
       return
    end if
    g = max(0d0, g0)
    if (g <= 0d0) g = (4d0*ti/9d0)**(1d0/3d0)
    do iter = 1, 30
       f = g*((9d0/4d0)*g**2 + g + 1d0) - ti
       df = (27d0/4d0)*g**2 + 2d0*g + 1d0
       if (df <= tiny) exit
       dg = -f/df
       g = max(0d0, g + dg)
       if (abs(f) <= max(1d-12, abs(ti)*1d-10)) exit
    end do
  end function gammac_root

  subroutine build_zones(model, conv_metric)
    type(mad_model), intent(inout) :: model
    real(dp), intent(in) :: conv_metric(:)

    integer, allocatable :: class(:), ztmp(:,:)
    integer :: i, nz, start, current, maxz

    allocate(class(model%np))
    class = 1
    do i = 1, model%npi
       if (conv_metric(i) > 1d-8) class(i) = 2
    end do
    do i = model%npi + 1, model%np
       class(i) = 3
    end do

    maxz = model%np
    allocate(ztmp(3,maxz))
    nz = 0
    start = 1
    current = class(1)
    do i = 2, model%np
       if (class(i) /= current) then
          nz = nz + 1
          ztmp(1,nz) = start
          ztmp(2,nz) = i - 1
          ztmp(3,nz) = current
          start = i
          current = class(i)
       end if
    end do
    nz = nz + 1
    ztmp(1,nz) = start
    ztmp(2,nz) = model%np
    ztmp(3,nz) = current

    call merge_tiny_zones(ztmp, nz, model%npi)

    model%nz = nz
    allocate(model%zones(3,nz))
    model%zones = ztmp(:,1:nz)
    deallocate(ztmp, class)
  end subroutine build_zones

  subroutine merge_tiny_zones(zones, nz, npi)
    integer, intent(inout) :: zones(:,:), nz
    integer, intent(in) :: npi

    integer :: k, target

    k = 1
    do while (k <= nz)
       if (zones(2,k) - zones(1,k) + 1 < 3 .and. nz > 1) then
          if (zones(1,k) > npi) then
             k = k + 1
             cycle
          end if
          if (k == 1) then
             target = 2
             zones(1,target) = zones(1,k)
          else
             target = k - 1
             zones(2,target) = zones(2,k)
          end if
          if (k < nz) zones(:,k:nz-1) = zones(:,k+1:nz)
          nz = nz - 1
       else
          k = k + 1
       end if
    end do
  end subroutine merge_tiny_zones

  subroutine write_madmod(filename, model, adiabatic_only)
    character(len=*), intent(in) :: filename
    type(mad_model), intent(in) :: model
    logical, intent(in) :: adiabatic_only

    integer :: io, ios
    integer :: np1, np
    character(len=8), parameter :: compact_magic = 'SISMOAD2'

    np = model%np
    np1 = model%npi + 1
    open(newunit=io, file=trim(filename), status='replace', action='write', &
         form='unformatted', access='stream', convert='little_endian', iostat=ios)
    if (ios /= 0) then
       write(*,'(a,a,a,i0)') 'mesa2SISMO: cannot write ', trim(filename), ' iostat=', ios
       stop 1
    end if

    if (adiabatic_only) then
       if (model%npa /= 0 .or. model%npi /= model%np) then
          close(io, status='delete')
          call fatal('internal error: compact adiabatic model contains atmosphere points.')
       end if
       write(io) compact_magic
       write(io) model%npi, model%npa, model%np, model%nz, model%step
       write(io) model%mass, model%radius, model%grav
       write(io) model%zones, model%r, model%m, model%rho, model%P, &
            model%Gam1, model%n2
    else
       write(io) model%npi, model%npa, model%np, model%nz, model%step
       write(io) model%mass, model%radius, model%teff, model%lum, model%age, &
            model%logg, model%lgOverL, model%alphaConv, model%alphaOver, &
            model%X0, model%Z0, model%diff, model%grav, model%arad, model%sigma
       write(io) model%zones, model%r, model%m, model%rho, model%T, model%L, model%P, &
            model%Cv, model%Gam1, model%Gam31, model%Prho, model%PT, model%Cp, &
            model%Cprho, model%CpT, model%Q, model%Qrho, model%QT, model%kappa, &
            model%krho, model%kT, model%en, model%X, model%Z, model%gradP, &
            model%gradT, model%gradRad, model%gradAd, model%grad, model%GamC, &
            model%AlphaC, model%dm1, model%dm2, model%yH1, model%yH2, model%yHe3, &
            model%yHe4, model%yLi7, model%yC12, model%yC13, model%yN14, model%yN15, &
            model%yO16, model%yO17, model%yO18, model%yNe20, model%yOthers, &
            model%ZOthers, model%ZZOthers, model%MOthers, model%eg, model%mix, model%n2
       if (model%npa > 0) then
          write(io) model%tau(np1:np), model%rr(np1:np), model%mm(np1:np), &
               model%Pg(np1:np), model%Cvg(np1:np), model%Gam1g(np1:np), &
               model%Gam31g(np1:np), model%Prhog(np1:np), model%PTg(np1:np), &
               model%Cpg(np1:np), model%Cprhog(np1:np), model%CpTg(np1:np), &
               model%Qg(np1:np), model%Qrhog(np1:np), model%QTg(np1:np), &
               model%accrad(np1:np), model%dTdTe(np1:np), model%dTdg(np1:np), &
               model%PR(np1:np), model%gradPR(np1:np)
       end if
    end if
    close(io)
  end subroutine write_madmod

  function header_value(p, name, fallback) result(val)
    type(mesa_profile), intent(in) :: p
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: fallback
    real(dp) :: val
    integer :: i

    val = fallback
    do i = 1, p%nhead
       if (same_name(p%hname(i), name)) then
          if (p%hnum(i) .and. is_good(p%hval(i))) val = p%hval(i)
          return
       end if
    end do
  end function header_value

  function col_surface(p, idx) result(val)
    type(mesa_profile), intent(in) :: p
    integer, intent(in) :: idx
    real(dp) :: val

    if (idx > 0) then
       val = p%data(idx,1)
    else
       val = 0d0
    end if
  end function col_surface

  function optional_col(p, idx, row, fallback) result(val)
    type(mesa_profile), intent(in) :: p
    integer, intent(in) :: idx, row
    real(dp), intent(in) :: fallback
    real(dp) :: val

    if (idx > 0) then
       val = p%data(idx,row)
       if (.not. is_good(val)) val = fallback
    else
       val = fallback
    end if
  end function optional_col

  function require_col(p, name) result(idx)
    type(mesa_profile), intent(in) :: p
    character(len=*), intent(in) :: name
    integer :: idx

    idx = find_col(p, name)
    if (idx <= 0) call fatal('missing required MESA profile column: '//trim(name))
  end function require_col

  function find_first_col(p, name1, name2) result(idx)
    type(mesa_profile), intent(in) :: p
    character(len=*), intent(in) :: name1, name2
    integer :: idx

    idx = find_col(p, name1)
    if (idx <= 0) idx = find_col(p, name2)
  end function find_first_col

  function find_col(p, name) result(idx)
    type(mesa_profile), intent(in) :: p
    character(len=*), intent(in) :: name
    integer :: idx, i

    idx = 0
    do i = 1, p%ncol
       if (same_name(p%cname(i), name)) then
          idx = i
          return
       end if
    end do
  end function find_col

  logical function same_name(a, b)
    character(len=*), intent(in) :: a, b
    character(len=128) :: aa, bb

    aa = adjustl(a)
    bb = adjustl(b)
    call lower(aa)
    call lower(bb)
    same_name = trim(aa) == trim(bb)
  end function same_name

  logical function is_good(x)
    real(dp), intent(in) :: x
    is_good = ieee_is_finite(x)
  end function is_good

  subroutine split_words(line, tokens, ntok)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: tokens(:)
    integer, intent(out) :: ntok

    integer :: n, i, i1, i2

    tokens = ''
    ntok = 0
    n = len_trim(line)
    i = 1
    do while (i <= n)
       do while (i <= n .and. line(i:i) == ' ')
          i = i + 1
       end do
       if (i > n) exit
       i1 = i
       do while (i <= n .and. line(i:i) /= ' ')
          i = i + 1
       end do
       i2 = i - 1
       if (ntok < size(tokens)) then
          ntok = ntok + 1
          tokens(ntok) = adjustl(line(i1:i2))
       end if
    end do
  end subroutine split_words

  subroutine lower(str)
    character(len=*), intent(inout) :: str
    integer :: i, c

    do i = 1, len_trim(str)
       c = iachar(str(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) str(i:i) = achar(c + 32)
    end do
  end subroutine lower

  subroutine fatal(message)
    character(len=*), intent(in) :: message

    write(*,'(a,a)') 'mesa2SISMO: ', trim(message)
    stop 1
  end subroutine fatal

end program mesa2SISMO
