! =============================================================================
! SISMO 2.0 -- Split Inhomogeneous Solver for Modelling Oscillations.
! Adiabatic stellar oscillations with a split (inhomogeneous) Poisson
! treatment on a scalar second-order Cowling core.
!
! Pipeline per degree l:
!   1. STURM SCAN of the Cowling operator (parallel sweep + bisection of the
!      parity of negative pivots + negative face denominators) -> the complete,
!      labeled spectrum (g, f, p) on the full model grid.  No seeding, no
!      asymptotics (the asymptotic relation only sets the frequency window).
!   2. SPLIT REFINEMENT of each mode: the frozen potential phi is recomputed
!      from the mechanics by the inhomogeneous shooting Poisson solve at every
!      step (for l=1 closed by Takata's first integral when takata_closure is
!      enabled in sismo.conf) and drives the Cowling operator through the
!      right-hand side; a clamped Newton update converges the frequency.
! =============================================================================
program sismo
  use sismo_types, only : stellar_model, iteration_config
  use sismo_config_io, only : find_default_sismo_config, load_sismo_config
  use sismo_io, only : read_intSISMO_model
  use sismo_poisson, only : set_takata_closure
  use sismo_core2, only : run_core2
  implicit none

  character(len=2048) :: structure_file, output_base, config_file, arg
  type(stellar_model) :: model
  type(iteration_config) :: config
  integer :: narg
  integer :: l_min, l_max, g_min, g_max

  narg = command_argument_count()
  if (narg == 1) then
     call get_command_argument(1, arg)
     if (trim(arg) == '--help' .or. trim(arg) == '-h') call usage(0)
  end if

  if (narg < 2 .or. narg > 3) then
     if (narg > 3) then
        write(*,*) 'SISMO: runtime options have moved to the configuration file.'
     end if
     call usage(1)
  end if

  call get_command_argument(1, structure_file)
  call get_command_argument(2, output_base)

  if (narg == 3) then
     call get_command_argument(3, config_file)
     if (index(trim(config_file), '--') == 1) then
        write(*,*) 'SISMO: runtime options have moved to the configuration file.'
        call usage(1)
     end if
  else
     call find_default_sismo_config(config_file)
  end if

  call load_sismo_config(config_file, config, l_min, l_max, g_min, g_max)
  call read_intSISMO_model(structure_file, model)
  call set_takata_closure(config%takata_closure)
  call run_core2(model, config, l_min, l_max, g_min, g_max, output_base)

contains

  subroutine usage(status)
    integer, intent(in) :: status

    write(*,*) 'Usage: sismo <model.osc.mod|sta1.d> <output_base> [config_file]'
    write(*,*) 'All degree, order, scan, Poisson, and convergence settings are read'
    write(*,*) 'from the configuration file. With no explicit file, SISMO searches:'
    write(*,*) '  1. SISMO_CONFIG'
    write(*,*) '  2. ./sismo.conf'
    write(*,*) '  3. sismo.conf beside the executable'
    write(*,*) 'Only -h and --help remain as command-line options.'
    if (status == 0) stop
    stop 1, quiet=.true.
  end subroutine usage

end program sismo
