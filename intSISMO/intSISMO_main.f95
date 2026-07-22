program intSISMO
  use iso_c_binding, only : c_int
  use madmodlib
  use intSISMO_lib, only : write_osc_model

  implicit none

  interface
     function c_getpid() bind(C, name='getpid') result(pid)
       import :: c_int
       integer(c_int) :: pid
     end function c_getpid
  end interface

  integer, parameter :: UNIT_STDOUT = 6
  integer, parameter :: PATH_LEN = 4096
  integer, parameter :: MAX_OSC_POINTS = 2000000

  integer :: narg, ios, ldat, argument_length, points_written, clock_count
  integer :: N6, N6_requested, NN1, NN1a, nz_osc, GRID_STEP
  character(len=PATH_LEN) :: dat, basefile, madmod, oscfile
  character(len=32) :: CN6, CSTEP, CMODE, GRID_MODE
  character(len=96) :: transaction_tag
  character(len=:), allocatable :: gridfile, osc_temp, grid_temp, lock_path
  logical :: exst, osc_ok, grid_ok, grid_owned, commit_ok, path_ok
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

  call get_command_argument(1, length=argument_length, status=ios)
  if (ios /= 0 .or. argument_length < 1) then
     write(UNIT_STDOUT,*) 'intSISMO: invalid first command-line argument.'
     stop 1
  end if
  if (argument_length > len(dat)) then
     write(UNIT_STDOUT,*) 'intSISMO: input path is too long: ', argument_length, &
          ' characters; maximum is ', len(dat)
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

  if (len_trim(basefile) > len(madmod)-len('.madmod') .or. &
       len_trim(basefile) > len(oscfile)-len('.osc.mod')) then
     write(UNIT_STDOUT,*) 'intSISMO: model basename is too long to append output suffixes.'
     stop 1
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
  if (ios /= 0 .or. N6 < 8 .or. N6 > MAX_OSC_POINTS) then
     write(UNIT_STDOUT,*) 'intSISMO: invalid OSC grid size: ', trim(CN6)
     write(UNIT_STDOUT,*) 'intSISMO: supported range is 8..', MAX_OSC_POINTS
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
     if (ios /= 0 .or. GRID_STEP < 1 .or. GRID_STEP > MAX_OSC_POINTS) then
        write(UNIT_STDOUT,*) 'intSISMO: GRID_STEP must be in 1..', &
             MAX_OSC_POINTS, ': ', trim(CSTEP)
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
  if (N6 > MAX_OSC_POINTS) then
     write(UNIT_STDOUT,*) 'intSISMO: adjusted OSC grid exceeds the safety limit: ', N6
     stop 1
  end if

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

  gridfile = trim(basefile)//'.grid.d'
  call system_clock(count=clock_count)
  write(transaction_tag,'(A,I0,A,I0)') '.sismo-', int(c_getpid()), '-', clock_count
  osc_temp = trim(oscfile)//trim(transaction_tag)//'.tmp'
  grid_temp = gridfile//trim(transaction_tag)//'.tmp'
  lock_path = trim(basefile)//'.sismo-output.lock'

  call require_absent_transaction_path(osc_temp, path_ok)
  if (.not. path_ok) then
     deallocate(zones_osc)
     call DeleteMadmod(star)
     stop 1
  end if
  call require_absent_transaction_path(grid_temp, path_ok)
  if (.not. path_ok) then
     deallocate(zones_osc)
     call DeleteMadmod(star)
     stop 1
  end if

  call write_osc_model(osc_temp, star, NN1, NN1a, N6, zones_osc, nz_osc, &
       GRID_MODE, osc_ok, points_written)
  if (.not. osc_ok) then
     deallocate(zones_osc)
     call DeleteMadmod(star)
     stop 1
  end if

  call write_grid_step_metadata(grid_temp, GRID_STEP, N6_requested, N6, &
       GRID_MODE, grid_ok, grid_owned)
  if (.not. grid_ok) then
    call remove_regular_file(osc_temp, path_ok)
     if (.not. path_ok) write(UNIT_STDOUT,*) &
         'intSISMO: WARNING - could not remove temporary file ', trim(osc_temp)
     if (grid_owned) then
        call remove_regular_file(grid_temp, path_ok)
        if (.not. path_ok) write(UNIT_STDOUT,*) &
             'intSISMO: WARNING - could not remove temporary file ', trim(grid_temp)
     end if
     deallocate(zones_osc)
     call DeleteMadmod(star)
     stop 1
  end if

  call commit_output_pair(osc_temp, trim(oscfile), grid_temp, gridfile, &
       lock_path, trim(transaction_tag), commit_ok)
  if (.not. commit_ok) then
     call remove_regular_file(osc_temp, path_ok)
     if (.not. path_ok) write(UNIT_STDOUT,*) &
          'intSISMO: WARNING - could not remove temporary file ', trim(osc_temp)
     call remove_regular_file(grid_temp, path_ok)
     if (.not. path_ok) write(UNIT_STDOUT,*) &
          'intSISMO: WARNING - could not remove temporary file ', trim(grid_temp)
     deallocate(zones_osc)
     call DeleteMadmod(star)
     stop 1
  end if

  write(UNIT_STDOUT,*) 'OSC model written: ', trim(oscfile), &
       ' points=', points_written, ' zones=1'
  write(UNIT_STDOUT,*) 'Grid metadata written: ', trim(gridfile)

  deallocate(zones_osc)
  call DeleteMadmod(star)

contains

  subroutine print_usage()
    implicit none

    write(UNIT_STDOUT,*) 'Usage: intSISMO <input base/.madmod/.osc.mod/.mod> <OSC grid size> [GRID_STEP] [GRID_MODE]'
    write(UNIT_STDOUT,*) '       GRID_STEP: integer in 1..', MAX_OSC_POINTS, &
         ' used as the output point-count multiple'
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

  subroutine write_grid_step_metadata(filename, grid_step, requested_grid, &
       target_grid, grid_mode, output_ok, output_owned)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: grid_step, requested_grid, target_grid
    character(len=*), intent(in) :: grid_mode
    logical, intent(out) :: output_ok
    logical, intent(out) :: output_owned

    integer :: unit, ios
    character(len=512) :: iomsg
    logical :: cleanup_ok

    output_ok = .false.
    output_owned = .false.
    iomsg = ''
    open(newunit=unit, file=trim(filename), status='new', action='write', &
         iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       write(UNIT_STDOUT,'(A,A,A,I0,A,A)') 'intSISMO: ERROR - cannot open grid metadata ', &
            trim(filename), ' (iostat=', ios, '): ', trim(iomsg)
       return
    end if
    output_owned = .true.

    write(unit,'(A)',iostat=ios,iomsg=iomsg) '# intSISMO grid metadata'
    if (ios /= 0) then
       call grid_metadata_failure(unit, filename, 'header', ios, iomsg, cleanup_ok)
       output_owned = .not. cleanup_ok
       return
    end if
    write(unit,'(A)',iostat=ios,iomsg=iomsg) &
         '# GRID_STEP, requested_target, adjusted_target, GRID_MODE'
    if (ios /= 0) then
       call grid_metadata_failure(unit, filename, 'column description', ios, iomsg, cleanup_ok)
       output_owned = .not. cleanup_ok
       return
    end if
    write(unit,*,iostat=ios,iomsg=iomsg) grid_step, requested_grid, target_grid, trim(grid_mode)
    if (ios /= 0) then
       call grid_metadata_failure(unit, filename, 'data record', ios, iomsg, cleanup_ok)
       output_owned = .not. cleanup_ok
       return
    end if

    close(unit, iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       call grid_metadata_failure(unit, filename, 'close', ios, iomsg, cleanup_ok)
       output_owned = .not. cleanup_ok
       return
    end if
    output_ok = .true.
  end subroutine write_grid_step_metadata


  subroutine grid_metadata_failure(unit, filename, context, ios, iomsg, removed)
    implicit none

    integer, intent(in) :: unit, ios
    character(len=*), intent(in) :: filename, context, iomsg
    logical, intent(out) :: removed

    call discard_grid_metadata(unit, filename, removed)
    write(UNIT_STDOUT,'(A,A,A,A,A,I0,A,A)') 'intSISMO: ERROR - grid metadata ', &
         trim(context), ' failed for ', trim(filename), ' (iostat=', ios, '): ', trim(iomsg)
  end subroutine grid_metadata_failure


  subroutine discard_grid_metadata(unit, filename, removed)
    implicit none

    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename
    logical, intent(out) :: removed

    integer :: cleanup_unit, cleanup_ios
    logical :: exists, opened

    removed = .false.
    inquire(unit=unit, opened=opened)
    if (opened) close(unit, status='delete', iostat=cleanup_ios)

    inquire(file=trim(filename), exist=exists)
    if (.not. exists) then
       removed = .true.
       return
    end if
    open(newunit=cleanup_unit, file=trim(filename), status='old', &
         action='readwrite', iostat=cleanup_ios)
    if (cleanup_ios == 0) close(cleanup_unit, status='delete', iostat=cleanup_ios)
    inquire(file=trim(filename), exist=exists)
    removed = .not. exists
    if (.not. removed) write(UNIT_STDOUT,*) &
         'intSISMO: WARNING - retained incomplete grid metadata ', trim(filename)
  end subroutine discard_grid_metadata


  subroutine require_absent_transaction_path(path, available)
    implicit none

    character(len=*), intent(in) :: path
    logical, intent(out) :: available

    logical :: exists, regular, query_ok

    call path_information(path, exists, regular, query_ok)
    available = query_ok .and. .not. exists
    if (.not. query_ok) then
       write(UNIT_STDOUT,*) 'intSISMO: ERROR - cannot inspect transaction path ', trim(path)
    else if (exists) then
       write(UNIT_STDOUT,*) 'intSISMO: ERROR - transaction path already exists: ', trim(path)
    end if
  end subroutine require_absent_transaction_path


  subroutine commit_output_pair(osc_tmp, osc_final, grid_tmp, grid_final, &
       lock_filename, tag, committed)
    implicit none

    character(len=*), intent(in) :: osc_tmp, osc_final, grid_tmp, grid_final
    character(len=*), intent(in) :: lock_filename, tag
    logical, intent(out) :: committed

    character(len=:), allocatable :: osc_backup, grid_backup
    character(len=512) :: iomsg
    integer :: lock_unit, ios
    logical :: osc_exists, grid_exists, regular, query_ok
    logical :: backup_exists, backed_osc, backed_grid
    logical :: installed_osc, installed_grid, cleanup_ok, rollback_ok

    committed = .false.
    backed_osc = .false.
    backed_grid = .false.
    installed_osc = .false.
    installed_grid = .false.
    osc_backup = trim(osc_final)//trim(tag)//'.bak'
    grid_backup = trim(grid_final)//trim(tag)//'.bak'

    iomsg = ''
    open(newunit=lock_unit, file=trim(lock_filename), status='new', &
         action='write', iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
       write(UNIT_STDOUT,'(A,A,A,I0,A,A)') &
            'intSISMO: ERROR - cannot acquire output lock ', trim(lock_filename), &
            ' (iostat=', ios, '): ', trim(iomsg)
       write(UNIT_STDOUT,*) &
            'intSISMO: remove a stale lock only after confirming no conversion is running.'
       return
    end if
    write(lock_unit,'(A)',iostat=ios,iomsg=iomsg) trim(tag)
    if (ios /= 0) then
       write(UNIT_STDOUT,*) 'intSISMO: ERROR - cannot initialize output lock: ', trim(iomsg)
       call release_output_lock(lock_unit, lock_filename, cleanup_ok)
       return
    end if

    call path_information(osc_final, osc_exists, regular, query_ok)
    if (.not. query_ok .or. (osc_exists .and. .not. regular)) then
       write(UNIT_STDOUT,*) &
            'intSISMO: ERROR - OSC destination is not an ordinary file: ', trim(osc_final)
       call release_output_lock(lock_unit, lock_filename, cleanup_ok)
       return
    end if
    call path_information(grid_final, grid_exists, regular, query_ok)
    if (.not. query_ok .or. (grid_exists .and. .not. regular)) then
       write(UNIT_STDOUT,*) &
            'intSISMO: ERROR - grid destination is not an ordinary file: ', trim(grid_final)
       call release_output_lock(lock_unit, lock_filename, cleanup_ok)
       return
    end if

    call path_information(osc_backup, backup_exists, regular, query_ok)
    if (.not. query_ok .or. backup_exists) then
       write(UNIT_STDOUT,*) &
            'intSISMO: ERROR - OSC backup path is unavailable: ', trim(osc_backup)
       call release_output_lock(lock_unit, lock_filename, cleanup_ok)
       return
    end if
    call path_information(grid_backup, backup_exists, regular, query_ok)
    if (.not. query_ok .or. backup_exists) then
       write(UNIT_STDOUT,*) &
            'intSISMO: ERROR - grid backup path is unavailable: ', trim(grid_backup)
       call release_output_lock(lock_unit, lock_filename, cleanup_ok)
       return
    end if

    if (osc_exists) then
       call rename_regular_file(osc_final, osc_backup, ios)
       if (ios /= 0) then
          write(UNIT_STDOUT,*) 'intSISMO: ERROR - cannot back up ', trim(osc_final)
          call release_output_lock(lock_unit, lock_filename, cleanup_ok)
          return
       end if
       backed_osc = .true.
    end if
    if (grid_exists) then
       call rename_regular_file(grid_final, grid_backup, ios)
       if (ios /= 0) then
          write(UNIT_STDOUT,*) 'intSISMO: ERROR - cannot back up ', trim(grid_final)
          call rollback_output_pair(osc_final, grid_final, osc_backup, grid_backup, &
               backed_osc, backed_grid, installed_osc, installed_grid, rollback_ok)
          call release_output_lock(lock_unit, lock_filename, cleanup_ok)
          return
       end if
       backed_grid = .true.
    end if

    call rename_regular_file(osc_tmp, osc_final, ios)
    if (ios /= 0) then
       write(UNIT_STDOUT,*) 'intSISMO: ERROR - cannot install ', trim(osc_final)
       call rollback_output_pair(osc_final, grid_final, osc_backup, grid_backup, &
            backed_osc, backed_grid, installed_osc, installed_grid, rollback_ok)
       call release_output_lock(lock_unit, lock_filename, cleanup_ok)
       return
    end if
    installed_osc = .true.

    call rename_regular_file(grid_tmp, grid_final, ios)
    if (ios /= 0) then
       write(UNIT_STDOUT,*) 'intSISMO: ERROR - cannot install ', trim(grid_final)
       call rollback_output_pair(osc_final, grid_final, osc_backup, grid_backup, &
            backed_osc, backed_grid, installed_osc, installed_grid, rollback_ok)
       call release_output_lock(lock_unit, lock_filename, cleanup_ok)
       return
    end if
    installed_grid = .true.

    if (backed_osc) then
       call remove_regular_file(osc_backup, cleanup_ok)
       if (.not. cleanup_ok) write(UNIT_STDOUT,*) &
            'intSISMO: WARNING - retained old OSC backup at ', trim(osc_backup)
    end if
    if (backed_grid) then
       call remove_regular_file(grid_backup, cleanup_ok)
       if (.not. cleanup_ok) write(UNIT_STDOUT,*) &
            'intSISMO: WARNING - retained old grid backup at ', trim(grid_backup)
    end if
    call release_output_lock(lock_unit, lock_filename, cleanup_ok)
    if (.not. cleanup_ok) write(UNIT_STDOUT,*) &
         'intSISMO: WARNING - output lock could not be removed: ', trim(lock_filename)

    committed = .true.
  end subroutine commit_output_pair


  subroutine rollback_output_pair(osc_final, grid_final, osc_backup, grid_backup, &
       backed_osc, backed_grid, installed_osc, installed_grid, rollback_ok)
    implicit none

    character(len=*), intent(in) :: osc_final, grid_final, osc_backup, grid_backup
    logical, intent(inout) :: backed_osc, backed_grid, installed_osc, installed_grid
    logical, intent(out) :: rollback_ok

    integer :: ios
    logical :: removed

    rollback_ok = .true.
    if (installed_grid) then
       call remove_regular_file(grid_final, removed)
       if (removed) then
          installed_grid = .false.
       else
          rollback_ok = .false.
          write(UNIT_STDOUT,*) &
               'intSISMO: ERROR - could not remove incomplete output ', trim(grid_final)
       end if
    end if
    if (installed_osc) then
       call remove_regular_file(osc_final, removed)
       if (removed) then
          installed_osc = .false.
       else
          rollback_ok = .false.
          write(UNIT_STDOUT,*) &
               'intSISMO: ERROR - could not remove incomplete output ', trim(osc_final)
       end if
    end if

    if (backed_grid) then
       call rename_regular_file(grid_backup, grid_final, ios)
       if (ios == 0) then
          backed_grid = .false.
       else
          rollback_ok = .false.
          write(UNIT_STDOUT,*) &
               'intSISMO: ERROR - old grid output remains at ', trim(grid_backup)
       end if
    end if
    if (backed_osc) then
       call rename_regular_file(osc_backup, osc_final, ios)
       if (ios == 0) then
          backed_osc = .false.
       else
          rollback_ok = .false.
          write(UNIT_STDOUT,*) &
               'intSISMO: ERROR - old OSC output remains at ', trim(osc_backup)
       end if
    end if

    if (.not. rollback_ok) write(UNIT_STDOUT,*) &
         'intSISMO: ERROR - rollback was incomplete; retained backup paths are reported above.'
  end subroutine rollback_output_pair


  subroutine rename_regular_file(old_path, new_path, status)
    implicit none

    character(len=*), intent(in) :: old_path, new_path
    integer, intent(out) :: status

    ! GNU Fortran RENAME is sufficient here because intSISMO is built with
    ! gfortran. Keeping it in one wrapper makes that dependency explicit.
    call rename(trim(old_path), trim(new_path), status)
  end subroutine rename_regular_file


  subroutine remove_regular_file(path, removed)
    implicit none

    character(len=*), intent(in) :: path
    logical, intent(out) :: removed

    integer :: ios
    logical :: exists, regular, query_ok

    call path_information(path, exists, regular, query_ok)
    if (.not. query_ok) then
       removed = .false.
       return
    end if
    if (.not. exists) then
       removed = .true.
       return
    end if
    if (.not. regular) then
       removed = .false.
       return
    end if

    ! GNU Fortran UNLINK removes the directory entry directly, so cleanup does
    ! not depend on reopening a readable file. The project is gfortran-only.
    call unlink(trim(path), ios)
    removed = ios == 0
  end subroutine remove_regular_file


  subroutine release_output_lock(unit, lock_filename, released)
    implicit none

    integer, intent(in) :: unit
    character(len=*), intent(in) :: lock_filename
    logical, intent(out) :: released

    integer :: ios
    logical :: opened

    inquire(unit=unit, opened=opened)
    ios = 0
    if (opened) close(unit, status='delete', iostat=ios)
    released = ios == 0
    if (.not. released) write(UNIT_STDOUT,*) &
         'intSISMO: WARNING - lock remains at ', trim(lock_filename)
  end subroutine release_output_lock


  subroutine path_information(path, exists, regular, query_ok)
    implicit none

    character(len=*), intent(in) :: path
    logical, intent(out) :: exists, regular, query_ok

    logical :: is_link, is_present

    call path_test('-L', path, is_link, query_ok)
    if (.not. query_ok) then
       exists = .false.
       regular = .false.
       return
    end if
    if (is_link) then
       exists = .true.
       regular = .false.
       return
    end if

    call path_test('-e', path, is_present, query_ok)
    if (.not. query_ok) then
       exists = .false.
       regular = .false.
       return
    end if
    exists = is_present
    if (.not. exists) then
       regular = .false.
       return
    end if

    call path_test('-f', path, regular, query_ok)
  end subroutine path_information


  subroutine path_test(predicate, path, matches, query_ok)
    implicit none

    character(len=*), intent(in) :: predicate, path
    logical, intent(out) :: matches, query_ok

    integer :: cmdstat, exitstat
    character(len=512) :: cmdmsg
    character(len=:), allocatable :: command

    command = '/bin/test '//trim(predicate)//' '//shell_quote(trim(path))
    cmdmsg = ''
    call execute_command_line(command, wait=.true., exitstat=exitstat, &
         cmdstat=cmdstat, cmdmsg=cmdmsg)
    query_ok = cmdstat == 0 .and. (exitstat == 0 .or. exitstat == 1)
    matches = query_ok .and. exitstat == 0
    if (.not. query_ok) write(UNIT_STDOUT,*) &
         'intSISMO: ERROR - path inspection failed: ', trim(cmdmsg)
  end subroutine path_test


  function shell_quote(text) result(quoted)
    implicit none

    character(len=*), intent(in) :: text
    character(len=:), allocatable :: quoted

    integer :: i

    quoted = achar(39)
    do i = 1, len_trim(text)
       if (text(i:i) == achar(39)) then
          quoted = quoted//achar(39)//achar(34)//achar(39)//achar(34)//achar(39)
       else
          quoted = quoted//text(i:i)
       end if
    end do
    quoted = quoted//achar(39)
  end function shell_quote

end program intSISMO
