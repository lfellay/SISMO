module sismo_types
  use sismo_precision, only : dp, czero
  implicit none
  private

  public :: stellar_model, mode_frequency, reference_mode, source_terms, mode_solution
  public :: iteration_config, mode_result
  public :: allocate_model, allocate_solution, allocate_sources
  public :: zero_sources, copy_model_metadata

  type :: stellar_model
     integer :: n = 0
     integer :: grid_step = 1
     integer :: nrac = 0
     integer :: nnuc = 0
     integer :: nzc = 0
     integer :: zc(50,2) = 0
     real(dp) :: lte = 0.0_dp, logg = 0.0_dp, feh = 0.0_dp
     real(dp) :: luminosity = 0.0_dp, mass = 0.0_dp, radius = 0.0_dp
     real(dp) :: taurac = 0.0_dp, omega0 = 0.0_dp
     real(dp) :: tdy = 1.0_dp, thk = 1.0_dp, qrh0 = 0.0_dp, fdy = 0.0_dp
     real(dp), allocatable :: x(:), d(:), temp(:), qx3(:), rho(:), drhdxx(:), aosc(:)
     real(dp), allocatable :: pressure(:), pgas(:), gamma1(:)
     real(dp), allocatable :: pt(:), prg(:), ptg(:), gamma3m1(:)
     real(dp), allocatable :: xdlt(:), d2ltx(:), dltx(:)
     real(dp), allocatable :: eeps(:), eps3rho(:), lrl1(:), lrl2(:), lum1(:)
     real(dp), allocatable :: krho(:), kt(:), ar(:), art(:), ki(:)
     real(dp), allocatable :: dtt(:), d2tt(:), dtte(:), d2tte(:), dtg(:), d2tg(:)
     real(dp), allocatable :: zeta0(:), phi0(:), ga0(:)
     real(dp), allocatable :: lta(:), dlt(:), kappa(:)
     real(dp), allocatable :: cpv(:), cprs(:), cpt(:), q(:), qrs(:), qt(:)
     real(dp), allocatable :: fcdsdx1(:), fcdsdx2(:), gamc1(:), gamc2(:)
     real(dp), allocatable :: fconv1(:), fconv2(:), dlcvdx(:)
     real(dp), allocatable :: pdsdx(:), pturb1(:), pturb2(:), dlpdlr(:)
  end type stellar_model

  type :: mode_frequency
     integer :: l = 0, m = 0, k = 0, n = 0
     real(dp) :: ltilde = 0.0_dp
     real(dp) :: omega_ad = 0.0_dp
     real(dp) :: freq_com_cd = 0.0_dp
     real(dp) :: freq_in_cd = 0.0_dp
     real(dp) :: omega_rot = 0.0_dp
  end type mode_frequency

  type :: reference_mode
     integer :: l = 0, m = 0, k = 0, n = 0
     real(dp) :: freq_com_cd = 0.0_dp
     real(dp) :: freq_in_cd = 0.0_dp
     real(dp) :: ltilde = 0.0_dp
     complex(dp) :: omega = czero
     real(dp) :: omega_ad = 0.0_dp
  end type reference_mode

  type :: source_terms
     integer :: n = 0
     complex(dp), allocatable :: radial(:), continuity(:), horizontal(:)
  end type source_terms

  type :: mode_solution
     integer :: n = 0
     integer :: l = 0
     complex(dp) :: omega = czero
     complex(dp), allocatable :: xi_r(:), dpp(:), xi_h(:)
     complex(dp), allocatable :: phi(:), theta(:), dlum(:)
  end type mode_solution

  type :: iteration_config
     integer :: max_iter = 150                  ! outer (frequency) iterations per mode
     real(dp) :: tol = 1.0d-8                   ! relative frequency convergence
     real(dp) :: eps_g = 2.5_dp                 ! asymptotic phase (scan WINDOW only, never seeds)
     integer :: sturm_scan_points = 2000         ! endpoint-inclusive samples per degree
     real(dp) :: poisson_relax = 0.4_dp         ! Picard damping of the frozen-phi refresh
     integer :: inner_poisson_iters = 1         ! phi refreshes per frequency step (cap)
     real(dp) :: poisson_inner_tol = 0.0_dp     ! >0: rate-corrected phi-convergence criterion
     logical :: takata_closure = .true.         ! l=1: first-integral (J=0) amplitude closure
     logical :: use_poisson = .true.            ! split self-gravity on/off (off = Cowling)
     logical :: write_eigenfunctions = .false.  ! retain and write mechanical mode shapes
  end type iteration_config

  type :: mode_result
     type(mode_frequency) :: mode
     type(mode_solution) :: solution
     integer :: iterations = 0
     integer :: final_stride = 1
     logical :: converged = .false.
     real(dp) :: inertia = 0.0_dp
     real(dp) :: work = 0.0_dp
     integer :: num = 0              ! classified radial order (Takata/Lee/Scuflaire)
     real(dp) :: err_takata = 0.0_dp ! Takata dipole-identity residual (l=1 quality check)
     real(dp), allocatable :: x(:)   ! grid the solution is on (for eigenfunction output)
  end type mode_result

contains

  subroutine allocate_model(model, n)
    type(stellar_model), intent(inout) :: model
    integer, intent(in) :: n

    model%n = n
    if (allocated(model%x)) then
       deallocate(model%x, model%d, model%temp, model%qx3, model%rho, model%drhdxx, model%aosc)
       deallocate(model%pressure, model%pgas, model%gamma1)
       deallocate(model%pt, model%prg, model%ptg, model%gamma3m1)
       deallocate(model%xdlt, model%d2ltx, model%dltx)
       deallocate(model%eeps, model%eps3rho, model%lrl1, model%lrl2, model%lum1)
       deallocate(model%krho, model%kt, model%ar, model%art, model%ki)
       deallocate(model%dtt, model%d2tt, model%dtte, model%d2tte, model%dtg, model%d2tg)
       deallocate(model%zeta0, model%phi0, model%ga0)
       deallocate(model%lta, model%dlt, model%kappa)
       deallocate(model%cpv, model%cprs, model%cpt, model%q, model%qrs, model%qt)
       deallocate(model%fcdsdx1, model%fcdsdx2, model%gamc1, model%gamc2)
       deallocate(model%fconv1, model%fconv2, model%dlcvdx)
       deallocate(model%pdsdx, model%pturb1, model%pturb2, model%dlpdlr)
    end if
    allocate(model%x(n), model%d(n), model%temp(n), model%qx3(n), model%rho(n), model%drhdxx(n), model%aosc(n))
    allocate(model%pressure(n), model%pgas(n), model%gamma1(n))
    allocate(model%pt(n), model%prg(n), model%ptg(n), model%gamma3m1(n))
    allocate(model%xdlt(n), model%d2ltx(n), model%dltx(n))
    allocate(model%eeps(n), model%eps3rho(n), model%lrl1(n), model%lrl2(n), model%lum1(n))
    allocate(model%krho(n), model%kt(n), model%ar(n), model%art(n), model%ki(n))
    allocate(model%dtt(n), model%d2tt(n), model%dtte(n), model%d2tte(n), model%dtg(n), model%d2tg(n))
    allocate(model%zeta0(n), model%phi0(n), model%ga0(n))
    allocate(model%lta(n), model%dlt(n), model%kappa(n))
    allocate(model%cpv(n), model%cprs(n), model%cpt(n), model%q(n), model%qrs(n), model%qt(n))
    allocate(model%fcdsdx1(n), model%fcdsdx2(n), model%gamc1(n), model%gamc2(n))
    allocate(model%fconv1(n), model%fconv2(n), model%dlcvdx(n))
    allocate(model%pdsdx(n), model%pturb1(n), model%pturb2(n), model%dlpdlr(n))

    model%x = 0.0_dp; model%d = 0.0_dp; model%temp = 0.0_dp
    model%qx3 = 0.0_dp; model%rho = 0.0_dp; model%drhdxx = 0.0_dp; model%aosc = 0.0_dp
    model%pressure = 0.0_dp; model%pgas = 0.0_dp; model%gamma1 = 0.0_dp
    model%pt = 0.0_dp; model%prg = 0.0_dp; model%ptg = 0.0_dp; model%gamma3m1 = 0.0_dp
    model%xdlt = 0.0_dp; model%d2ltx = 0.0_dp; model%dltx = 0.0_dp
    model%eeps = 0.0_dp; model%eps3rho = 0.0_dp; model%lrl1 = 0.0_dp
    model%lrl2 = 0.0_dp; model%lum1 = 0.0_dp; model%krho = 0.0_dp
    model%kt = 0.0_dp; model%ar = 0.0_dp; model%art = 0.0_dp; model%ki = 0.0_dp
    model%dtt = 0.0_dp; model%d2tt = 0.0_dp; model%dtte = 0.0_dp
    model%d2tte = 0.0_dp; model%dtg = 0.0_dp; model%d2tg = 0.0_dp
    model%zeta0 = 0.0_dp; model%phi0 = 0.0_dp; model%ga0 = 0.0_dp
    model%lta = 0.0_dp; model%dlt = 0.0_dp; model%kappa = 0.0_dp
    model%cpv = 0.0_dp; model%cprs = 0.0_dp; model%cpt = 0.0_dp
    model%q = 0.0_dp; model%qrs = 0.0_dp; model%qt = 0.0_dp
    model%fcdsdx1 = 0.0_dp; model%fcdsdx2 = 0.0_dp
    model%gamc1 = 0.0_dp; model%gamc2 = 0.0_dp
    model%fconv1 = 0.0_dp; model%fconv2 = 0.0_dp; model%dlcvdx = 0.0_dp
    model%pdsdx = 0.0_dp; model%pturb1 = 0.0_dp; model%pturb2 = 0.0_dp
    model%dlpdlr = 0.0_dp
  end subroutine allocate_model

  subroutine allocate_solution(solution, n)
    type(mode_solution), intent(inout) :: solution
    integer, intent(in) :: n

    solution%n = n
    if (allocated(solution%xi_r)) then
       deallocate(solution%xi_r, solution%dpp, solution%xi_h)
       deallocate(solution%phi, solution%theta, solution%dlum)
    end if
    allocate(solution%xi_r(n), solution%dpp(n), solution%xi_h(n))
    allocate(solution%phi(n), solution%theta(n), solution%dlum(n))
    solution%xi_r = czero; solution%dpp = czero; solution%xi_h = czero
    solution%phi = czero; solution%theta = czero; solution%dlum = czero
  end subroutine allocate_solution

  subroutine allocate_sources(sources, n)
    type(source_terms), intent(inout) :: sources
    integer, intent(in) :: n

    sources%n = n
    if (allocated(sources%radial)) deallocate(sources%radial, sources%continuity, sources%horizontal)
    allocate(sources%radial(n), sources%continuity(n), sources%horizontal(n))
    call zero_sources(sources)
  end subroutine allocate_sources

  subroutine zero_sources(sources)
    type(source_terms), intent(inout) :: sources

    if (allocated(sources%radial)) sources%radial = czero
    if (allocated(sources%continuity)) sources%continuity = czero
    if (allocated(sources%horizontal)) sources%horizontal = czero
  end subroutine zero_sources

  subroutine copy_model_metadata(src, dst)
    type(stellar_model), intent(in) :: src
    type(stellar_model), intent(inout) :: dst

    dst%nrac = src%nrac
    dst%grid_step = src%grid_step
    dst%nnuc = src%nnuc
    dst%nzc = src%nzc
    dst%zc = src%zc
    dst%lte = src%lte
    dst%logg = src%logg
    dst%feh = src%feh
    dst%luminosity = src%luminosity
    dst%mass = src%mass
    dst%radius = src%radius
    dst%taurac = src%taurac
    dst%omega0 = src%omega0
    dst%tdy = src%tdy
    dst%thk = src%thk
    dst%qrh0 = src%qrh0
    dst%fdy = src%fdy
  end subroutine copy_model_metadata

end module sismo_types
