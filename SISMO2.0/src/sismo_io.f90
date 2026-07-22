module sismo_io
  use sismo_precision, only : dp, pi, tiny, czero
  use sismo_types, only : stellar_model, mode_frequency, reference_mode, mode_result, allocate_model
  implicit none
  private

  public :: read_intSISMO_model, read_frequency_file, read_osc_reference
  public :: open_result_file, write_result_row, write_comparison_file
  public :: write_eigenfunction_file

  integer, parameter :: line_len = 1024

contains

  subroutine read_intSISMO_model(filename, model)
    character(len=*), intent(in) :: filename
    type(stellar_model), intent(out) :: model

    integer :: unit, ios, i, layer, zeta_cols, structure_cols, pturb_cols
    real(dp) :: v3(3), v4(4), zeta_values(3)
    character(len=line_len) :: line
    logical :: old_extra_thermal, old_composition

    if (has_suffix(trim(filename), '.mod')) then
       call read_osc_binary_model(filename, model)
       return
    end if

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal('cannot open intSISMO structure '//trim(filename))

    model%grid_step = read_grid_step_sidecar(filename)
    call read_header_until_log_teff(unit, line, model)
    if (index(line, ':') > 0) then
       call read_reals_after_colon(line, v4, 4, 'stellar header')
    else
       call read_line(unit, line, 'stellar header values')
       call read_real_values(line, v4, 4, 'stellar header values')
    end if
    model%lte = v4(1)
    model%logg = v4(2)
    model%feh = v4(3)
    model%luminosity = v4(4)

    call read_labeled_real(unit, 'rayon', model%radius, 'radius')
    call read_labeled_real(unit, 'masse', model%mass, 'mass')
    call read_labeled_real(unit, 'profondeur optique raccordement', model%taurac, 'raccord optical depth')
    call read_labeled_integer(unit, 'Nbre de couches total', model%n, 'total mesh size')
    call read_labeled_integer(unit, 'Nbre de couches avec react.nuc.', model%nnuc, 'nuclear mesh size')
    call read_labeled_integer(unit, 'couche de raccordement interieur-atmosphere', model%nrac, 'raccord layer')

    call read_line(unit, line, 'omega label')
    if (contains_ci(line, 'omega, tdy')) then
       call read_line(unit, line, 'omega values')
       call read_real_values(line, v4, 4, 'omega values')
       model%omega0 = v4(1)
       model%tdy = max(tiny, v4(2))
       model%thk = v4(3)
       model%qrh0 = v4(4)
    else if (contains_ci(line, 'ordre du mode')) then
       call read_line(unit, line, 'mode degree value')
       call read_line(unit, line, 'old omega label')
       call read_line(unit, line, 'old omega value')
       read(line, *, iostat=ios) model%omega0
       if (ios /= 0) model%omega0 = 0.0_dp
       call read_line(unit, line, 'old dynamical label')
       call require_contains(line, '1/t_dyn', 'old dynamical label')
       call read_line(unit, line, 'old dynamical values')
       call read_real_values(line, v3, 2, 'old dynamical values')
       model%fdy = v3(1)
       model%tdy = 1.0d6/(max(tiny, model%fdy)*2.0_dp*pi)
       model%thk = v3(2)*model%tdy
       call read_labeled_real(unit, 'qrh0', model%qrh0, 'central qrho limit')
    else
       call fatal('unknown intSISMO header layout near '//trim(line))
    end if
    model%fdy = 1.0d6/(model%tdy*2.0_dp*pi)

    call read_line(unit, line, 'convective-zone label')
    if (.not. contains_ci(line, 'ZC:') .and. .not. contains_ci(line, 'Zones convectives')) then
       call fatal('could not find convective-zone label')
    end if
    call read_line(unit, line, 'convective-zone values')
    call parse_convective_zones(line, model)

    call skip_column_labels(unit, zeta_cols, structure_cols, old_extra_thermal, old_composition, pturb_cols)
    call allocate_model(model, model%n)

    do i = 1, model%n
       call read_integer_record(unit, layer, 'layer index')
       if (layer /= i-1) then
          write(*,*) 'SISMO reader: layer index mismatch: got ', layer, ' expected ', i-1
          stop 1
       end if

       call read_real_record(unit, v3, 3, 'x d T')
       model%x(i) = v3(1); model%d(i) = v3(2); model%temp(i) = v3(3)

       v4 = 0.0_dp
       call read_real_record(unit, v4, structure_cols, 'qx3 rho drhdxx Aosc')
       model%qx3(i) = v4(1); model%rho(i) = v4(2); model%drhdxx(i) = v4(3)
       if (structure_cols >= 4) model%aosc(i) = v4(4)

       call read_real_record(unit, v3, 3, 'P Pg Gamma1')
       model%pressure(i) = v3(1); model%pgas(i) = v3(2); model%gamma1(i) = v3(3)

       call read_real_record(unit, v3, 3, 'PT Prg PTg')
       model%pt(i) = v3(1); model%prg(i) = v3(2); model%ptg(i) = v3(3)

       call read_real_record(unit, v3, 3, 'Gamma3m1 xdlT d2lTx')
       model%gamma3m1(i) = v3(1); model%xdlt(i) = v3(2); model%d2ltx(i) = v3(3)

       call read_real_record(unit, v3, 3, 'dlTx eeps eps3rho')
       model%dltx(i) = v3(1); model%eeps(i) = v3(2); model%eps3rho(i) = v3(3)

       call read_real_record(unit, v4, 4, 'LrL1 LrL2 L krho')
       model%lrl1(i) = v4(1); model%lrl2(i) = v4(2); model%lum1(i) = v4(3); model%krho(i) = v4(4)

       call read_real_record(unit, v4, 4, 'kT ar art ki')
       model%kt(i) = v4(1); model%ar(i) = v4(2); model%art(i) = v4(3); model%ki(i) = v4(4)

       call read_real_record(unit, v3, 3, 'dTt d2Tt dTTe')
       model%dtt(i) = v3(1); model%d2tt(i) = v3(2); model%dtte(i) = v3(3)

       if (old_extra_thermal) then
          call read_real_record(unit, v3, 3, 'old d2TTe dTg d2Tg')
          model%d2tte(i) = v3(1); model%dtg(i) = v3(2); model%d2tg(i) = v3(3)
          call read_real_record(unit, v3, 3, 'old dPTe dPg dPt')
       end if
       if (old_composition) call skip_old_composition_records(unit)

       zeta_values = 0.0_dp
       call read_real_record(unit, zeta_values, zeta_cols, 'zeta/phi/ga')
       model%zeta0(i) = zeta_values(1)
       if (zeta_cols >= 3) then
          model%phi0(i) = zeta_values(2)
          model%ga0(i) = zeta_values(3)
       else
          model%phi0(i) = 0.0_dp
          model%ga0(i) = 1.0_dp
       end if

       call read_real_record(unit, v3, 3, 'lta dlt kappa')
       model%lta(i) = v3(1); model%dlt(i) = v3(2); model%kappa(i) = v3(3)

       call read_real_record(unit, v3, 3, 'cpv cprs cpT')
       model%cpv(i) = v3(1); model%cprs(i) = v3(2); model%cpt(i) = v3(3)

       call read_real_record(unit, v3, 3, 'Q Qrs QT')
       model%q(i) = v3(1); model%qrs(i) = v3(2); model%qt(i) = v3(3)

       call read_real_record(unit, v4, 4, 'Fcdsdx GamC')
       model%fcdsdx1(i) = v4(1); model%fcdsdx2(i) = v4(2); model%gamc1(i) = v4(3); model%gamc2(i) = v4(4)

       call read_real_record(unit, v3, 3, 'fconv dlcvdx')
       model%fconv1(i) = v3(1); model%fconv2(i) = v3(2); model%dlcvdx(i) = v3(3)

       v3 = 0.0_dp
       call read_real_record(unit, v3, pturb_cols, 'Pdsdx pturb')
       model%pdsdx(i) = v3(1)
       model%pturb1(i) = v3(2)
       if (pturb_cols >= 3) then
          model%pturb2(i) = v3(3)
       else
          model%pturb2(i) = v3(2)
       end if
    end do

    call compute_pressure_gradient(model)
    close(unit)
  end subroutine read_intSISMO_model

  subroutine read_osc_binary_model(filename, model)
    character(len=*), intent(in) :: filename
    type(stellar_model), intent(out) :: model

    integer :: unit, ios, n, nz, i
    integer, allocatable :: z(:,:)
    real(dp) :: globals(3), dyn_freq
    real(dp) :: r, enclosed_mass, pressure, rho, gamma1, abrunt

    open(newunit=unit, file=trim(filename), status='old', action='read', &
         form='unformatted', access='stream', convert='little_endian', iostat=ios)
    if (ios /= 0) call fatal('cannot open OSC binary model '//trim(filename))

    read(unit, iostat=ios) n, nz, globals
    if (ios /= 0) call fatal('cannot read OSC binary model header '//trim(filename))
    allocate(z(2, nz))
    read(unit, iostat=ios) z
    if (ios /= 0) call fatal('cannot read OSC binary model zones '//trim(filename))

    call allocate_model(model, n)
    ! Preserve intSISMO's GRID_STEP conversion metadata when its sidecar is
    ! available. SISMO itself always solves on the complete imported grid.
    model%grid_step = read_grid_step_sidecar(filename)
    model%radius = globals(1)
    model%mass = globals(2)
    dyn_freq = sqrt(max(tiny, globals(3)*globals(2)/globals(1)**3))
    model%fdy = dyn_freq/(2.0_dp*pi)*1.0d6
    model%tdy = 1.0d6/(max(tiny, model%fdy)*2.0_dp*pi)
    ! The OSC z payload is a generic segmentation table, not a classified
    ! convective-zone table. Do not mislabel it as model%zc metadata.
    model%nzc = 0
    model%zc = 0

    do i = 1, n
       read(unit, iostat=ios) r, enclosed_mass, pressure, rho, gamma1, abrunt
       if (ios /= 0) call fatal('cannot read OSC binary model point '//trim(filename))
       model%x(i) = r/max(tiny, model%radius)
       model%qx3(i) = model%radius**3/max(tiny, model%mass) * enclosed_mass/max(tiny, r**3)
       model%rho(i) = 4.0_dp*pi*model%radius**3/max(tiny, model%mass) * rho
       model%pressure(i) = model%radius/(globals(3)*model%mass) * pressure/max(tiny, rho) * model%rho(i)
       model%pgas(i) = model%pressure(i)
       model%gamma1(i) = gamma1
       model%aosc(i) = model%radius**2 * abrunt/max(tiny, r)
    end do
    close(unit)
    call compute_pressure_gradient(model)
  end subroutine read_osc_binary_model

  subroutine read_frequency_file(filename, modes, nmode, max_modes)
    character(len=*), intent(in) :: filename
    type(mode_frequency), allocatable, intent(out) :: modes(:)
    integer, intent(out) :: nmode
    integer, intent(in), optional :: max_modes

    integer :: unit, ios, cap, limit
    character(len=line_len) :: line
    type(mode_frequency), allocatable :: tmp(:)

    limit = huge(0)
    if (present(max_modes)) limit = max_modes

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal('cannot open OSC frequency file '//trim(filename))

    cap = max(16, min(limit, 1024))
    allocate(tmp(cap))
    nmode = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       if (nmode >= limit) exit
       if (nmode == cap) call grow_modes(tmp, cap)
       if (parse_frequency_line(line, tmp(nmode+1))) nmode = nmode + 1
    end do
    close(unit)

    allocate(modes(nmode))
    if (nmode > 0) modes = tmp(1:nmode)
  end subroutine read_frequency_file

  subroutine read_osc_reference(filename, refs, nref)
    character(len=*), intent(in) :: filename
    type(reference_mode), allocatable, intent(out) :: refs(:)
    integer, intent(out) :: nref

    type(mode_frequency), allocatable :: modes(:)
    integer :: i

    call read_frequency_file(filename, modes, nref)
    allocate(refs(nref))
    do i = 1, nref
       refs(i)%l = modes(i)%l
       refs(i)%m = modes(i)%m
       refs(i)%k = modes(i)%k
       refs(i)%n = modes(i)%n
       refs(i)%freq_com_cd = modes(i)%freq_com_cd
       refs(i)%freq_in_cd = modes(i)%freq_in_cd
       refs(i)%ltilde = modes(i)%ltilde
       refs(i)%omega = cmplx(modes(i)%omega_ad, 0.0_dp, kind=dp)
       refs(i)%omega_ad = modes(i)%omega_ad
    end do
  end subroutine read_osc_reference

  subroutine open_result_file(base, unit)
    character(len=*), intent(in) :: base
    integer, intent(out) :: unit
    integer :: ios

    open(newunit=unit, file=trim(base)//'.sismo', status='replace', action='write', iostat=ios)
    if (ios /= 0) call fatal('cannot open output '//trim(base)//'.sismo')
    write(unit,'(A)') '# l m k n iter conv stride points omega_seed omega_ad Im_omega freq_cd growth inertia work num errTakata'
  end subroutine open_result_file

  subroutine write_result_row(unit, result)
    integer, intent(in) :: unit
    type(mode_result), intent(in) :: result

    real(dp) :: freq_cd
    freq_cd = 0.0_dp
    if (result%mode%omega_ad > tiny) then
       freq_cd = result%mode%freq_com_cd * real(result%solution%omega, dp) / result%mode%omega_ad
    end if
    write(unit,'(4I5,1X,I5,1X,L1,1X,I5,1X,I8,1P,7E18.9,1X,I6,1X,E14.6)') &
         result%mode%l, result%mode%m, result%mode%k, result%mode%n, &
         result%iterations, result%converged, result%final_stride, result%solution%n, &
         result%mode%omega_ad, real(result%solution%omega, dp), aimag(result%solution%omega), &
         freq_cd, aimag(result%solution%omega), result%inertia, result%work, &
         result%num, result%err_takata
  end subroutine write_result_row

  subroutine write_comparison_file(base, results, nresult, refs, nref)
    character(len=*), intent(in) :: base
    type(mode_result), intent(in) :: results(:)
    integer, intent(in) :: nresult
    type(reference_mode), intent(in), optional :: refs(:)
    integer, intent(in), optional :: nref

    integer :: unit, i, j, ios, occurrence, seen
    logical :: have_ref, matched_ref
    complex(dp) :: omega_ref
    real(dp) :: im_ref, dre_ref, dim_ref

    have_ref = present(refs) .and. present(nref)
    open(newunit=unit, file=trim(base)//'.comparison', status='replace', action='write', iostat=ios)
    if (ios /= 0) call fatal('cannot open comparison output')
    write(unit,'(A)') '# l m k n ref_found omega_SISMO omega_OSC omega_seed d_seed d_osc rel_d_osc'
    do i = 1, nresult
       omega_ref = czero
       matched_ref = .false.
       if (have_ref) then
          occurrence = 0
          do j = 1, i
             if (results(j)%mode%l == results(i)%mode%l .and. results(j)%mode%m == results(i)%mode%m .and. &
                 results(j)%mode%k == results(i)%mode%k .and. results(j)%mode%n == results(i)%mode%n) then
                occurrence = occurrence + 1
             end if
          end do
          seen = 0
          do j = 1, nref
             if (refs(j)%l == results(i)%mode%l .and. refs(j)%m == results(i)%mode%m .and. &
                 refs(j)%k == results(i)%mode%k .and. refs(j)%n == results(i)%mode%n) then
                seen = seen + 1
                if (seen == occurrence) then
                   omega_ref = refs(j)%omega
                   matched_ref = .true.
                   exit
                end if
             end if
          end do
       end if

       im_ref = 0.0_dp
       dre_ref = 0.0_dp
       dim_ref = 0.0_dp
       if (matched_ref) then
          im_ref = real(omega_ref, dp)
          dre_ref = real(results(i)%solution%omega - omega_ref, dp)
          dim_ref = dre_ref/max(tiny, abs(real(omega_ref, dp)))
       end if

       write(unit,'(5I5,1P6E18.9)') results(i)%mode%l, results(i)%mode%m, &
            results(i)%mode%k, results(i)%mode%n, merge(1, 0, matched_ref), &
            real(results(i)%solution%omega, dp), &
            im_ref, results(i)%mode%omega_ad, &
            real(results(i)%solution%omega, dp) - results(i)%mode%omega_ad, &
            dre_ref, dim_ref
    end do
    close(unit)
    if (have_ref) call write_nearest_comparison_file(base, results, nresult, refs, nref)
  end subroutine write_comparison_file

  subroutine write_nearest_comparison_file(base, results, nresult, refs, nref)
    character(len=*), intent(in) :: base
    type(mode_result), intent(in) :: results(:)
    integer, intent(in) :: nresult
    type(reference_mode), intent(in) :: refs(:)
    integer, intent(in) :: nref

    integer :: unit, i, j, ios, best
    real(dp) :: dist, best_dist

    open(newunit=unit, file=trim(base)//'.nearest.comparison', status='replace', action='write', iostat=ios)
    if (ios /= 0) call fatal('cannot open nearest comparison output')
    write(unit,'(A)') '# l m k n nearest_n omega_SISMO omega_OSC d_osc rel_d_osc'
    do i = 1, nresult
       best = 0
       best_dist = huge(1.0_dp)
       do j = 1, nref
          if (refs(j)%l /= results(i)%mode%l) cycle
          if (refs(j)%m /= results(i)%mode%m) cycle
          if (refs(j)%k /= results(i)%mode%k) cycle
          if (results(i)%mode%n < 0 .and. refs(j)%n >= 0) cycle
          dist = abs(real(results(i)%solution%omega, dp) - real(refs(j)%omega, dp))
          if (dist < best_dist) then
             best_dist = dist
             best = j
          end if
       end do

       if (best > 0) then
          write(unit,'(5I5,1P4E18.9)') results(i)%mode%l, results(i)%mode%m, &
               results(i)%mode%k, results(i)%mode%n, refs(best)%n, &
               real(results(i)%solution%omega, dp), real(refs(best)%omega, dp), &
               real(results(i)%solution%omega - refs(best)%omega, dp), &
               real(results(i)%solution%omega - refs(best)%omega, dp)/max(tiny, abs(real(refs(best)%omega, dp)))
       else
          write(unit,'(5I5,1P4E18.9)') results(i)%mode%l, results(i)%mode%m, &
               results(i)%mode%k, results(i)%mode%n, 0, &
               real(results(i)%solution%omega, dp), 0.0_dp, 0.0_dp, 0.0_dp
       end if
    end do
    close(unit)
  end subroutine write_nearest_comparison_file

  subroutine write_eigenfunction_file(model, base, result)
    type(stellar_model), intent(in) :: model
    character(len=*), intent(in) :: base
    type(mode_result), intent(in) :: result

    integer :: unit, ios, i, n, offset, jmid
    logical :: midpoint_grid, use_result_grid
    real(dp) :: x_out, sd, sp, s2
    complex(dp) :: phi_out, theta_out
    character(len=512) :: filename

    ! Preferred path: the result carries the exact grid the solution was solved on
    ! (e.g. the strided/restricted midpoint grid), so we can write directly on it.
    use_result_grid = allocated(result%x)
    if (use_result_grid) use_result_grid = (size(result%x) == result%solution%n)

    offset = 0
    midpoint_grid = .false.
    if (.not. use_result_grid .and. result%solution%n /= model%n) then
       if (result%solution%n == model%n-1 .and. model%x(1) <= tiny) then
          offset = 1
       else if (result%solution%n == 2*model%n-1) then
          midpoint_grid = .true.
       else
          return
       end if
    end if

    write(filename,'(A,"_l",I0,"_m",I0,"_k",I0,"_n",I0,".eig")') trim(base), &
         result%mode%l, result%mode%m, result%mode%k, result%mode%n
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call fatal('cannot open eigenfunction output '//trim(filename))
    write(unit,'(A)') '# i x Re_xi_r Re_xi_h Re_dpp Re_phi Im_xi_r Im_xi_h Im_dpp Im_phi Re_theta Im_theta Re_phi_osc Im_phi_osc Re_theta_osc Im_theta_osc'
    n = result%solution%n
    do i = 1, n
       if (use_result_grid) then
          x_out = result%x(i)
       else if (midpoint_grid) then
          if (mod(i, 2) == 1) then
             x_out = model%x((i+1)/2)
          else
             jmid = i/2
             x_out = 0.5_dp*(model%x(jmid) + model%x(jmid+1))
          end if
       else
          x_out = model%x(i+offset)
       end if
       sd = regularization_scale(x_out, max(0, result%mode%l-1))
       sp = regularization_scale(x_out, max(0, result%mode%l))
       s2 = max(tiny, interp_real(model%x, model%qx3, x_out))
       phi_out = result%solution%phi(i)/cmplx(max(tiny, sp*s2), 0.0_dp, kind=dp)
       theta_out = result%solution%theta(i)/cmplx(max(tiny, sd*s2), 0.0_dp, kind=dp)
       write(unit,'(I8,1P15E18.9)') i, x_out, &
            real(result%solution%xi_r(i), dp), real(result%solution%xi_h(i), dp), &
            real(result%solution%dpp(i), dp), real(result%solution%phi(i), dp), &
            aimag(result%solution%xi_r(i)), aimag(result%solution%xi_h(i)), &
            aimag(result%solution%dpp(i)), aimag(result%solution%phi(i)), &
            real(result%solution%theta(i), dp), aimag(result%solution%theta(i)), &
            real(phi_out, dp), aimag(phi_out), real(theta_out, dp), aimag(theta_out)
    end do
    close(unit)
  end subroutine write_eigenfunction_file

  real(dp) function regularization_scale(x, power) result(scale)
    real(dp), intent(in) :: x
    integer, intent(in) :: power

    if (power <= 0) then
       scale = 1.0_dp
    else
       scale = max(x, 0.0_dp)**power
    end if
  end function regularization_scale

  real(dp) function interp_real(xgrid, values, x) result(value)
    real(dp), intent(in) :: xgrid(:), values(:), x

    integer :: i, n
    real(dp) :: dx, weight

    n = min(size(xgrid), size(values))
    if (n <= 0) then
       value = 0.0_dp
       return
    end if
    if (x <= xgrid(1)) then
       value = values(1)
       return
    end if
    if (x >= xgrid(n)) then
       value = values(n)
       return
    end if

    do i = 1, n-1
       if (x <= xgrid(i+1)) then
          dx = xgrid(i+1) - xgrid(i)
          if (abs(dx) <= tiny) then
             value = values(i)
          else
             weight = (x - xgrid(i))/dx
             value = (1.0_dp - weight)*values(i) + weight*values(i+1)
          end if
          return
       end if
    end do
    value = values(n)
  end function interp_real

  logical function parse_frequency_line(line, mode)
    character(len=*), intent(in) :: line
    type(mode_frequency), intent(out) :: mode

    integer :: ios
    real(dp) :: v(8)

    v = 0.0_dp
    read(line, *, iostat=ios) mode%l, mode%m, mode%k, mode%n, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8)
    if (ios /= 0) then
       read(line, *, iostat=ios) mode%l, mode%m, mode%k, mode%n, v(1), v(2), v(3), v(4), v(5), v(6)
    end if
    parse_frequency_line = (ios == 0)
    if (.not. parse_frequency_line) return
    mode%ltilde = v(1)
    mode%omega_ad = v(2)
    mode%freq_com_cd = v(3)
    mode%freq_in_cd = v(4)
    mode%omega_rot = v(6)
  end function parse_frequency_line

  subroutine compute_pressure_gradient(model)
    type(stellar_model), intent(inout) :: model
    integer :: i
    real(dp) :: dx

    if (model%n <= 1) return
    model%dlpdlr(1) = 0.0_dp
    do i = 2, model%n-1
       dx = max(tiny, model%x(i+1) - model%x(i-1))
       model%dlpdlr(i) = (log(max(tiny, model%pressure(i+1))) - &
            log(max(tiny, model%pressure(i-1)))) / dx
    end do
    model%dlpdlr(model%n) = model%dlpdlr(max(1, model%n-1))
  end subroutine compute_pressure_gradient

  subroutine parse_convective_zones(line, model)
    character(len=*), intent(in) :: line
    type(stellar_model), intent(inout) :: model

    integer :: ntok, ios, i, offset
    character(len=128) :: words(256)
    integer :: ints(256)

    call split_words(line, words, ntok)
    if (ntok < 1) call fatal('empty convective-zone values')
    do i = 1, ntok
       read(words(i), *, iostat=ios) ints(i)
       if (ios /= 0) call fatal('invalid convective-zone integer: '//trim(words(i)))
    end do
    model%nzc = ints(1)
    if (model%nzc < 0 .or. model%nzc > 50) call fatal('invalid convective-zone count')
    model%zc = 0
    if (model%nzc == 0) return

    if (ntok >= 2 + 2*model%nzc) then
       offset = 2
    else if (ntok >= 1 + 2*model%nzc) then
       offset = 1
    else
       call fatal('too few convective-zone bounds')
    end if
    do i = 1, model%nzc
       model%zc(i,1) = ints(offset + 2*i - 1)
       model%zc(i,2) = ints(offset + 2*i)
       if (model%zc(i,1) < 0 .or. model%zc(i,2) < model%zc(i,1) .or. &
            model%zc(i,2) > model%n) then
          call fatal('invalid convective-zone bounds')
       end if
    end do
  end subroutine parse_convective_zones

  subroutine skip_column_labels(unit, zeta_cols, structure_cols, old_extra_thermal, old_composition, pturb_cols)
    integer, intent(in) :: unit
    integer, intent(out) :: zeta_cols, structure_cols
    logical, intent(out) :: old_extra_thermal, old_composition
    integer, intent(out) :: pturb_cols

    character(len=line_len) :: line

    zeta_cols = 3
    structure_cols = 3
    old_extra_thermal = .false.
    old_composition = .false.
    pturb_cols = 3
    do
       call read_line(unit, line, 'column labels')
       if (len_trim(line) > 0) exit
    end do
    call require_contains(line, 'x1', 'first column label')
    do
       if (index(line, 'zeta0') > 0) then
          if (index(line, 'phi0') > 0 .or. index(line, 'ga0') > 0) then
             zeta_cols = 3
          else
             zeta_cols = 1
          end if
       end if
       if (contains_ci(line, 'qx3') .and. contains_ci(line, 'Aosc')) structure_cols = 4
       if (contains_ci(line, 'd2TTe') .or. contains_ci(line, 'dPTe')) old_extra_thermal = .true.
       if (contains_ci(line, 'H1') .or. contains_ci(line, 'He4') .or. contains_ci(line, 'elec')) old_composition = .true.
       if (contains_ci(line, 'Pdsdx1') .and. contains_ci(line, 'Pturb') .and. &
           .not. contains_ci(line, 'Pturb2')) pturb_cols = 2
       call read_line(unit, line, 'column labels')
       if (len_trim(line) == 0) exit
    end do
  end subroutine skip_column_labels

  subroutine skip_old_composition_records(unit)
    integer, intent(in) :: unit
    real(dp) :: v3(3)

    call read_real_record(unit, v3, 2, 'old H1 He4')
    call read_real_record(unit, v3, 3, 'old He3 H2 Be7')
    call read_real_record(unit, v3, 3, 'old Li7 C12 N13')
    call read_real_record(unit, v3, 3, 'old C13 N14 O15')
    call read_real_record(unit, v3, 3, 'old N15 O16 F17')
    call read_real_record(unit, v3, 3, 'old O17 elec eps')
  end subroutine skip_old_composition_records

  subroutine read_labeled_real(unit, expected, value, context)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: expected, context
    real(dp), intent(out) :: value
    character(len=line_len) :: line
    integer :: ios

    call read_line(unit, line, trim(context)//' label')
    call require_contains(line, expected, context)
    call read_line(unit, line, context)
    read(line, *, iostat=ios) value
    if (ios /= 0) call fatal('could not read '//trim(context))
  end subroutine read_labeled_real

  subroutine read_labeled_integer(unit, expected, value, context)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: expected, context
    integer, intent(out) :: value
    character(len=line_len) :: line
    integer :: ios

    call read_line(unit, line, trim(context)//' label')
    call require_contains(line, expected, context)
    call read_line(unit, line, context)
    read(line, *, iostat=ios) value
    if (ios /= 0) call fatal('could not read '//trim(context))
  end subroutine read_labeled_integer

  subroutine read_integer_record(unit, value, context)
    integer, intent(in) :: unit
    integer, intent(out) :: value
    character(len=*), intent(in) :: context
    character(len=line_len) :: line
    integer :: ios

    call read_line(unit, line, context)
    read(line, *, iostat=ios) value
    if (ios /= 0) call fatal('could not read '//trim(context))
  end subroutine read_integer_record

  subroutine read_real_record(unit, values, n, context)
    integer, intent(in) :: unit, n
    real(dp), intent(inout) :: values(:)
    character(len=*), intent(in) :: context
    character(len=line_len) :: line

    call read_line(unit, line, context)
    call read_real_values(line, values, n, context)
  end subroutine read_real_record

  subroutine read_reals_after_colon(line, values, n, context)
    character(len=*), intent(in) :: line, context
    real(dp), intent(inout) :: values(:)
    integer, intent(in) :: n
    integer :: pos

    pos = index(line, ':')
    if (pos <= 0) call fatal('missing colon in '//trim(context))
    call read_real_values(line(pos+1:), values, n, context)
  end subroutine read_reals_after_colon

  subroutine read_real_values(line, values, n, context)
    character(len=*), intent(in) :: line, context
    real(dp), intent(inout) :: values(:)
    integer, intent(in) :: n
    integer :: ios, j

    if (size(values) < n) call fatal('internal reader buffer too small')
    read(line, *, iostat=ios) (values(j), j=1,n)
    if (ios /= 0) call fatal('could not read '//trim(context))
  end subroutine read_real_values

  subroutine read_line(unit, line, context)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: line
    character(len=*), intent(in) :: context
    integer :: ios

    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) call fatal('unexpected end of file while reading '//trim(context))
  end subroutine read_line

  subroutine require_contains(line, expected, context)
    character(len=*), intent(in) :: line, expected, context

    if (.not. contains_ci(line, expected)) then
       write(*,*) 'SISMO reader expected ', trim(expected), ' while reading ', trim(context)
       write(*,*) 'got: ', trim(line)
       stop 1
    end if
  end subroutine require_contains

  subroutine read_header_until_log_teff(unit, line, model)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: line
    type(stellar_model), intent(inout) :: model

    do
       call read_line(unit, line, 'stellar header')
       call parse_grid_step_line(line, model%grid_step)
       if (contains_ci(line, 'log Teff')) return
    end do
  end subroutine read_header_until_log_teff

  subroutine parse_grid_step_line(line, grid_step)
    character(len=*), intent(in) :: line
    integer, intent(inout) :: grid_step

    integer :: ios, pos, value

    if (.not. contains_ci(line, 'MULTI_P')) return
    pos = index(line, ':')
    if (pos <= 0) return
    read(line(pos+1:), *, iostat=ios) value
    if (ios == 0 .and. value >= 1) grid_step = value
  end subroutine parse_grid_step_line

  integer function read_grid_step_sidecar(filename) result(grid_step)
    character(len=*), intent(in) :: filename

    character(len=512) :: sidecar

    call zone_filename_for_structure(filename, sidecar)
    grid_step = read_grid_step_from_sidecar(sidecar, .true.)
    if (grid_step > 1) return

    call grid_filename_for_structure(filename, sidecar)
    grid_step = read_grid_step_from_sidecar(sidecar, .false.)
  end function read_grid_step_sidecar

  integer function read_grid_step_from_sidecar(sidecar, needs_label) result(grid_step)
    character(len=*), intent(in) :: sidecar
    logical, intent(in) :: needs_label

    character(len=512) :: line
    integer :: unit, ios, value
    logical :: exists, want_values

    grid_step = 1
    inquire(file=trim(sidecar), exist=exists)
    if (.not. exists) return

    open(newunit=unit, file=trim(sidecar), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    want_values = .not. needs_label
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (want_values) then
          if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
          read(line, *, iostat=ios) value
          if (ios == 0 .and. value >= 1) grid_step = value
          exit
       end if
       if (contains_ci(line, 'MULTI_P')) want_values = .true.
    end do
    close(unit)
  end function read_grid_step_from_sidecar

  subroutine zone_filename_for_structure(filename, zonefile)
    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: zonefile

    integer :: lt

    zonefile = trim(filename)//'.zones.d'
    lt = len_trim(filename)
    if (has_suffix(trim(filename), '.sta1.d')) then
       if (lt > 7) zonefile = filename(1:lt-7)//'.zones.d'
    else if (has_suffix(trim(filename), '.d')) then
       if (lt > 2) zonefile = filename(1:lt-2)//'.zones.d'
    end if
  end subroutine zone_filename_for_structure

  subroutine grid_filename_for_structure(filename, gridfile)
    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: gridfile

    integer :: lt

    gridfile = trim(filename)//'.grid.d'
    lt = len_trim(filename)
    if (has_suffix(trim(filename), '.osc.mod')) then
       if (lt > 8) gridfile = filename(1:lt-8)//'.grid.d'
    else if (has_suffix(trim(filename), '.sta1.d')) then
       if (lt > 7) gridfile = filename(1:lt-7)//'.grid.d'
    else if (has_suffix(trim(filename), '.d')) then
       if (lt > 2) gridfile = filename(1:lt-2)//'.grid.d'
    end if
  end subroutine grid_filename_for_structure

  pure logical function contains_ci(line, expected) result(found)
    character(len=*), intent(in) :: line, expected
    character(len=len(line)) :: a
    character(len=len(expected)) :: b

    a = line
    b = expected
    call lower(a)
    call lower(b)
    found = index(a, trim(b)) > 0
  end function contains_ci

  pure logical function has_suffix(text, suffix) result(found)
    character(len=*), intent(in) :: text, suffix

    integer :: lt, ls

    lt = len_trim(text)
    ls = len_trim(suffix)
    found = lt >= ls
    if (found) found = text(lt-ls+1:lt) == suffix(1:ls)
  end function has_suffix

  pure subroutine lower(str)
    character(len=*), intent(inout) :: str
    integer :: i, c

    do i = 1, len_trim(str)
       c = iachar(str(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) str(i:i) = achar(c + 32)
    end do
  end subroutine lower

  subroutine split_words(line, words, ntok)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: words(:)
    integer, intent(out) :: ntok

    integer :: i, n, start
    logical :: in_word

    ntok = 0
    n = len_trim(line)
    in_word = .false.
    start = 1
    do i = 1, n+1
       if (i <= n) then
          if (line(i:i) /= ' ' .and. line(i:i) /= char(9)) then
             if (.not. in_word) then
                start = i
                in_word = .true.
             end if
             cycle
          end if
       end if
       if (in_word) then
          if (ntok >= size(words)) call fatal('too many values in one input record')
          if (i - start > len(words(1))) call fatal('convective-zone token is too long')
          ntok = ntok + 1
          words(ntok) = line(start:i-1)
          in_word = .false.
       end if
    end do
  end subroutine split_words

  subroutine grow_modes(modes, cap)
    type(mode_frequency), allocatable, intent(inout) :: modes(:)
    integer, intent(inout) :: cap
    type(mode_frequency), allocatable :: larger(:)

    allocate(larger(2*cap))
    larger(1:cap) = modes
    call move_alloc(larger, modes)
    cap = size(modes)
  end subroutine grow_modes

  subroutine fatal(message)
    character(len=*), intent(in) :: message
    write(*,*) 'SISMO: ', trim(message)
    stop 1
  end subroutine fatal

end module sismo_io
