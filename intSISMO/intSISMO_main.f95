program intSISMO
  use madmodlib
  use intSISMO_lib, only : write_osc_model

  implicit none

  integer, parameter :: UNIT_STDOUT = 6

  integer :: narg, ios, ldat
  integer :: N6, N6_requested, NN1, NN1a, nz_osc, GRID_STEP
  character(len=256) :: dat, basefile, madmod, oscfile
  character(len=32) :: CN6, CSTEP, CMODE, GRID_MODE
  logical :: exst
  type(tMadmod) :: star
  integer, allocatable :: zones_osc(:,:)

  write(UNIT_STDOUT,*) '---------------------------------------------------'
  write(UNIT_STDOUT,*) '**** intSISMO ****'

  narg = command_argument_count()
  if (narg == 1) then
     call get_command_argument(1, dat, status=ios)
     if (ios == 0 .and. (trim(dat) == '--help' .or. trim(dat) == '-h')) then
        call print_usage()
        stop
     end if
  end if
  if (narg < 2 .or. narg > 4) then
     call print_usage()
     stop 1
  end if

  call get_command_argument(1, dat, status=ios)
  if (ios /= 0 .or. len_trim(dat) == 0) then
     write(UNIT_STDOUT,*) 'intSISMO: invalid first command-line argument.'
     stop 1
  end if

  basefile = trim(dat)
  ldat = len_trim(basefile)
  if (ends_with(basefile, '.osc.mod')) then
     if (ldat <= len('.osc.mod')) call invalid_basefile(trim(dat))
     basefile = basefile(1:ldat-len('.osc.mod'))
  else if (ends_with(basefile, '.madmod')) then
     if (ldat <= len('.madmod')) call invalid_basefile(trim(dat))
     basefile = basefile(1:ldat-len('.madmod'))
  else if (ends_with(basefile, '.mod')) then
     if (ldat <= len('.mod')) call invalid_basefile(trim(dat))
     basefile = basefile(1:ldat-len('.mod'))
  end if

  madmod = trim(basefile)//'.madmod'
  oscfile = trim(basefile)//'.osc.mod'

  inquire(file=madmod, exist=exst)
  if (.not. exst) then
     write(UNIT_STDOUT,*) '** intSISMO: file ', trim(madmod), ' does not exist'
     stop 1
  end if

  call get_command_argument(2, CN6, status=ios)
  if (ios /= 0 .or. len_trim(CN6) == 0) then
     write(UNIT_STDOUT,*) 'intSISMO: missing/invalid OSC grid-size argument.'
     stop 1
  end if
  read(CN6, *, iostat=ios) N6
  if (ios /= 0 .or. N6 < 8) then
     write(UNIT_STDOUT,*) 'intSISMO: invalid OSC grid size: ', trim(CN6)
     stop 1
  end if
  N6_requested = N6

  GRID_STEP = 1
  GRID_MODE = 'radial'
  if (narg >= 3) then
     call get_command_argument(3, CSTEP, status=ios)
     if (ios /= 0 .or. len_trim(CSTEP) == 0) then
        write(UNIT_STDOUT,*) 'intSISMO: missing/invalid GRID_STEP argument.'
        stop 1
     end if
     read(CSTEP, *, iostat=ios) GRID_STEP
     if (ios /= 0 .or. GRID_STEP < 1) then
        write(UNIT_STDOUT,*) 'intSISMO: GRID_STEP must be an integer >= 1: ', trim(CSTEP)
       stop 1
     end if
  end if
  if (narg >= 4) then
     call get_command_argument(4, CMODE, status=ios)
     if (ios /= 0 .or. len_trim(CMODE) == 0) then
        write(UNIT_STDOUT,*) 'intSISMO: missing/invalid GRID_MODE argument.'
        stop 1
     end if
     GRID_MODE = normalize_grid_mode(CMODE)
     if (len_trim(GRID_MODE) == 0) then
        write(UNIT_STDOUT,*) 'intSISMO: GRID_MODE must be radial or bv: ', trim(CMODE)
        stop 1
     end if
  end if
  N6 = nearest_multiple(N6_requested, GRID_STEP)

  write(UNIT_STDOUT,*) '** intSISMO: reading .madmod file ..'
  call readmadmod(madmod, star)
  write(UNIT_STDOUT,*) '** intSISMO: reading .madmod file: done **'

  NN1 = star%npi
  NN1a = star%np
  nz_osc = star%nz
  if (nz_osc < 1) then
     write(UNIT_STDOUT,*) 'intSISMO: input model has no zone table.'
     call DeleteMadmod(star)
     stop 1
  end if

  allocate(zones_osc(3, nz_osc))
  zones_osc = star%zones(:,1:nz_osc)

  write(UNIT_STDOUT,*) 'OSC requested grid points = ', N6_requested
  write(UNIT_STDOUT,*) 'OSC target grid points    = ', N6
  write(UNIT_STDOUT,*) 'Grid point multiple       = ', GRID_STEP
  write(UNIT_STDOUT,*) 'Grid distribution mode    = ', trim(GRID_MODE)
  write(UNIT_STDOUT,*) 'OSC output file        = ', trim(oscfile)

  call write_osc_model(oscfile, star, NN1, NN1a, N6, zones_osc, nz_osc, GRID_MODE)
  call write_grid_step_metadata(basefile, GRID_STEP, N6_requested, N6, GRID_MODE)

  deallocate(zones_osc)
  call DeleteMadmod(star)

contains

  subroutine print_usage()
    implicit none

    write(UNIT_STDOUT,*) 'Usage: intSISMO <input base/.madmod/.osc.mod/.mod> <OSC grid size> [GRID_STEP] [GRID_MODE]'
    write(UNIT_STDOUT,*) '       GRID_STEP: any integer >= 1 used as the output point-count multiple'
    write(UNIT_STDOUT,*) '       GRID_MODE: radial (default) or bv'
  end subroutine print_usage

  logical function ends_with(text, suffix)
    implicit none
    character(len=*), intent(in) :: text, suffix
    integer :: lt, ls

    lt = len_trim(text)
    ls = len(suffix)
    ends_with = .false.
    if (lt < ls) return
    ends_with = text(lt-ls+1:lt) == suffix
  end function ends_with


  subroutine invalid_basefile(argument)
    implicit none
    character(len=*), intent(in) :: argument

    write(UNIT_STDOUT,*) 'intSISMO: cannot derive a model basename from ', trim(argument)
    stop 1
  end subroutine invalid_basefile

  integer function nearest_multiple(value, step)
    implicit none
    integer, intent(in) :: value, step
    integer :: lower, upper

    if (step <= 1) then
       nearest_multiple = max(8, value)
       return
    end if
    lower = max(step, (value/step)*step)
    upper = lower
    if (upper < value .and. upper <= huge(upper)-step) upper = upper + step
    if (value - lower <= upper - value) then
       nearest_multiple = lower
    else
       nearest_multiple = upper
    end if
    if (nearest_multiple < 8) nearest_multiple = ((8 + step - 1)/step)*step
  end function nearest_multiple

  function normalize_grid_mode(argument) result(mode)
    implicit none
    character(len=*), intent(in) :: argument
    character(len=32) :: mode
    character(len=32) :: lower

    lower = lowercase(adjustl(argument))
    mode = ''
    select case (trim(lower))
    case ('radial', 'legacy', 'default')
       mode = 'radial'
    case ('bv', 'brunt', 'brunt-vaisala', 'brunt_vaisala', 'bruntvaisala')
       mode = 'bv'
    end select
  end function normalize_grid_mode

  function lowercase(text) result(out)
    implicit none
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code

    out = text
    do i = 1, len(text)
       code = iachar(out(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) then
          out(i:i) = achar(code + iachar('a') - iachar('A'))
       end if
    end do
  end function lowercase

  subroutine write_grid_step_metadata(base, grid_step, requested_grid, target_grid, grid_mode)
    implicit none
    character(len=*), intent(in) :: base
    integer, intent(in) :: grid_step, requested_grid, target_grid
    character(len=*), intent(in) :: grid_mode

    integer :: unit, ios

    open(newunit=unit, file=trim(base)//'.grid.d', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       write(UNIT_STDOUT,*) 'intSISMO: WARNING - could not write grid metadata for ', trim(base)
       return
    end if
    write(unit,'(A)') '# intSISMO grid metadata'
    write(unit,'(A)') '# GRID_STEP, requested_target, adjusted_target, GRID_MODE'
    write(unit,*) grid_step, requested_grid, target_grid, trim(grid_mode)
    close(unit)
  end subroutine write_grid_step_metadata

end program intSISMO
