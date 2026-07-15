module sismo_config_io
  use, intrinsic :: iso_fortran_env, only : error_unit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use sismo_precision, only : dp
  use sismo_types, only : iteration_config
  implicit none
  private

  integer, parameter :: path_len = 2048
  integer, parameter :: line_len = 2048
  integer, parameter :: key_count = 14

  public :: find_default_sismo_config, load_sismo_config

contains

  subroutine find_default_sismo_config(filename)
    character(len=*), intent(out) :: filename

    character(len=path_len) :: env_value, candidate
    integer :: status
    logical :: exists, found

    filename = ''
    env_value = ''
    call get_environment_variable('SISMO_CONFIG', env_value, status=status)
    if (status == 0 .and. len_trim(env_value) > 0) then
       filename = trim(env_value)
       return
    else if (status == -1) then
       call config_error('', 0, 'SISMO_CONFIG is longer than the supported path length')
    end if

    candidate = 'sismo.conf'
    inquire(file=trim(candidate), exist=exists)
    if (exists) then
       filename = trim(candidate)
       return
    end if

    call executable_config_path(candidate, found)
    if (found) then
       filename = trim(candidate)
       return
    end if

    call config_error('', 0, 'no configuration file found; pass CONFIG as the third argument, '// &
         'set SISMO_CONFIG, create ./sismo.conf, or install sismo.conf beside the executable')
  end subroutine find_default_sismo_config

  subroutine executable_config_path(config_path, found)
    character(len=*), intent(out) :: config_path
    logical, intent(out) :: found

    character(len=path_len) :: executable, candidate, directory, token
    character(len=:), allocatable :: path_value
    integer :: slash, status, path_length, start, relative_colon
    logical :: exists, executable_exists, finished

    config_path = ''
    found = .false.
    executable = ''
    call get_command_argument(0, executable, status=status)
    if (status /= 0 .or. len_trim(executable) == 0) return

    slash = scan(trim(executable), '/', back=.true.)
    if (slash > 0) then
       candidate = executable(:slash)//'sismo.conf'
       inquire(file=trim(candidate), exist=exists)
       if (exists) then
          config_path = trim(candidate)
          found = .true.
       end if
       return
    end if

    call get_environment_variable('PATH', length=path_length, status=status)
    if (status /= 0 .or. path_length <= 0) return
    allocate(character(len=path_length) :: path_value)
    call get_environment_variable('PATH', path_value, status=status)
    if (status /= 0) return

    start = 1
    do
       if (start > len(path_value)) exit
       relative_colon = index(path_value(start:), ':')
       if (relative_colon == 0) then
          token = path_value(start:)
          finished = .true.
       else
          if (relative_colon == 1) then
             token = '.'
          else
             token = path_value(start:start+relative_colon-2)
          end if
          start = start + relative_colon
          finished = .false.
       end if
       if (len_trim(token) == 0) token = '.'

       directory = trim(token)
       candidate = trim(directory)//'/'//trim(executable)
       inquire(file=trim(candidate), exist=executable_exists)
       if (executable_exists) then
          candidate = trim(directory)//'/sismo.conf'
          inquire(file=trim(candidate), exist=exists)
          if (exists) then
             config_path = trim(candidate)
             found = .true.
             return
          end if
          ! PATH resolves an executable to its first match. Do not continue and
          ! accidentally load a configuration beside a different SISMO binary.
          return
       end if
       if (finished) exit
    end do
  end subroutine executable_config_path

  subroutine load_sismo_config(filename, config, l_min, l_max, g_min, g_max)
    character(len=*), intent(in) :: filename
    type(iteration_config), intent(out) :: config
    integer, intent(out) :: l_min, l_max, g_min, g_max

    character(len=line_len) :: line, value
    character(len=128) :: key
    character(len=512) :: iomsg
    logical :: seen(key_count), logical_value, ok
    integer :: unit, ios, line_number, equals, integer_value, key_index
    real(dp) :: real_value

    call set_release_defaults(config, l_min, l_max, g_min, g_max)
    seen = .false.
    open(newunit=unit, file=trim(filename), status='old', action='read', form='formatted', &
         iostat=ios, iomsg=iomsg)
    if (ios /= 0) call config_error(filename, 0, 'cannot open file: '//trim(iomsg))

    line_number = 0
    do
       read(unit, '(A)', iostat=ios, iomsg=iomsg) line
       if (ios < 0) exit
       line_number = line_number + 1
       if (ios > 0) call config_error(filename, line_number, 'cannot read line: '//trim(iomsg))
       call normalize_config_whitespace(line)
       call strip_config_comment(line)
       line = adjustl(line)
       if (len_trim(line) == 0) cycle

       equals = index(line, '=')
       if (equals <= 1 .or. equals >= len_trim(line)) then
          call config_error(filename, line_number, 'expected key = value')
       end if
       if (index(line(equals+1:), '=') > 0) then
          call config_error(filename, line_number, 'more than one equals sign')
       end if
       key = lowercase(trim(adjustl(line(:equals-1))))
       value = trim(adjustl(line(equals+1:)))
       if (len_trim(key) == 0 .or. len_trim(value) == 0) then
          call config_error(filename, line_number, 'empty key or value')
       end if

       select case (trim(key))
       case ('l_min')
          key_index = 1
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_integer(value, integer_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'integer')
          l_min = integer_value
       case ('l_max')
          key_index = 2
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_integer(value, integer_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'integer')
          l_max = integer_value
       case ('g_min')
          key_index = 3
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_integer(value, integer_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'integer')
          g_min = integer_value
       case ('g_max')
          key_index = 4
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_integer(value, integer_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'integer')
          g_max = integer_value
       case ('use_poisson')
          key_index = 5
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_logical(value, logical_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'boolean')
          config%use_poisson = logical_value
       case ('takata_closure')
          key_index = 6
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_logical(value, logical_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'boolean')
          config%takata_closure = logical_value
       case ('write_eigenfunctions')
          key_index = 7
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_logical(value, logical_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'boolean')
          config%write_eigenfunctions = logical_value
       case ('scan_points')
          key_index = 8
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_integer(value, integer_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'integer')
          config%sturm_scan_points = integer_value
       case ('eps_g')
          key_index = 9
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_real(value, real_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'real number')
          config%eps_g = real_value
       case ('max_iterations')
          key_index = 10
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_integer(value, integer_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'integer')
          config%max_iter = integer_value
       case ('tolerance')
          key_index = 11
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_real(value, real_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'real number')
          config%tol = real_value
       case ('inner_poisson_iterations')
          key_index = 12
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_integer(value, integer_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'integer')
          config%inner_poisson_iters = integer_value
       case ('poisson_relaxation')
          key_index = 13
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_real(value, real_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'real number')
          config%poisson_relax = real_value
       case ('poisson_inner_tolerance')
          key_index = 14
          call mark_key(seen, key_index, filename, line_number, key)
          call parse_real(value, real_value, ok)
          if (.not. ok) call invalid_value(filename, line_number, key, value, 'real number')
          config%poisson_inner_tol = real_value
       case default
          call config_error(filename, line_number, 'unknown configuration key: '//trim(key))
       end select
    end do
    close(unit)

    call validate_configuration(filename, config, l_min, l_max, g_min, g_max)
    write(*,'(A)') ' SISMO: loaded configuration '//trim(filename)
    write(*,'(A,I0,A,I0,A,I0,A,I0)') ' SISMO: l=', l_min, '..', l_max, '  orders=', g_min, '..', g_max
  end subroutine load_sismo_config

  subroutine set_release_defaults(config, l_min, l_max, g_min, g_max)
    type(iteration_config), intent(out) :: config
    integer, intent(out) :: l_min, l_max, g_min, g_max

    l_min = 1
    l_max = 2
    g_min = 1
    g_max = 108
    config%use_poisson = .true.
    config%takata_closure = .true.
    config%write_eigenfunctions = .false.
    config%sturm_scan_points = 2000
    config%eps_g = 2.5_dp
    config%max_iter = 150
    config%tol = 1.0d-8
    config%inner_poisson_iters = 1
    config%poisson_relax = 0.4_dp
    config%poisson_inner_tol = 0.0_dp
  end subroutine set_release_defaults

  subroutine validate_configuration(filename, config, l_min, l_max, g_min, g_max)
    character(len=*), intent(in) :: filename
    type(iteration_config), intent(in) :: config
    integer, intent(in) :: l_min, l_max, g_min, g_max

    if (l_min < 1) call config_error(filename, 0, 'l_min must be at least 1')
    if (l_max < l_min) call config_error(filename, 0, 'l_max must be greater than or equal to l_min')
    if (g_min < 1) call config_error(filename, 0, 'g_min must be at least 1')
    if (g_max < g_min) call config_error(filename, 0, 'g_max must be greater than or equal to g_min')
    if (config%sturm_scan_points < 2) call config_error(filename, 0, 'scan_points must be at least 2')
    if (config%eps_g < 0.0_dp) call config_error(filename, 0, 'eps_g must be non-negative')
    if (config%max_iter < 1) call config_error(filename, 0, 'max_iterations must be at least 1')
    if (config%tol <= 0.0_dp) call config_error(filename, 0, 'tolerance must be positive')
    if (config%inner_poisson_iters < 1) then
       call config_error(filename, 0, 'inner_poisson_iterations must be at least 1')
    end if
    if (config%poisson_relax <= 0.0_dp .or. config%poisson_relax > 1.0_dp) then
       call config_error(filename, 0, 'poisson_relaxation must be greater than 0 and at most 1')
    end if
    if (config%poisson_inner_tol < 0.0_dp) then
       call config_error(filename, 0, 'poisson_inner_tolerance must be non-negative')
    end if
  end subroutine validate_configuration

  subroutine mark_key(seen, index_value, filename, line_number, key)
    logical, intent(inout) :: seen(:)
    integer, intent(in) :: index_value, line_number
    character(len=*), intent(in) :: filename, key

    if (seen(index_value)) call config_error(filename, line_number, 'duplicate key: '//trim(key))
    seen(index_value) = .true.
  end subroutine mark_key

  subroutine parse_integer(text, value, ok)
    character(len=*), intent(in) :: text
    integer, intent(out) :: value
    logical, intent(out) :: ok

    character(len=1) :: extra
    integer :: ios, trailing_ios

    value = 0
    read(text, *, iostat=ios) value
    if (ios /= 0) then
       ok = .false.
       return
    end if
    read(text, *, iostat=trailing_ios) value, extra
    ok = trailing_ios < 0
  end subroutine parse_integer

  subroutine parse_real(text, value, ok)
    character(len=*), intent(in) :: text
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    character(len=1) :: extra
    integer :: ios, trailing_ios

    value = 0.0_dp
    read(text, *, iostat=ios) value
    if (ios /= 0) then
       ok = .false.
       return
    end if
    read(text, *, iostat=trailing_ios) value, extra
    ok = trailing_ios < 0 .and. ieee_is_finite(value)
  end subroutine parse_real

  subroutine parse_logical(text, value, ok)
    character(len=*), intent(in) :: text
    logical, intent(out) :: value, ok

    character(len=line_len) :: normalized

    normalized = lowercase(trim(adjustl(text)))
    select case (trim(normalized))
    case ('true', '.true.', 'yes', 'on', '1')
       value = .true.
       ok = .true.
    case ('false', '.false.', 'no', 'off', '0')
       value = .false.
       ok = .true.
    case default
       value = .false.
       ok = .false.
    end select
  end subroutine parse_logical

  subroutine strip_config_comment(line)
    character(len=*), intent(inout) :: line

    integer :: hash, bang, cut

    cut = len_trim(line) + 1
    hash = index(line, '#')
    bang = index(line, '!')
    if (hash > 0) cut = min(cut, hash)
    if (bang > 0) cut = min(cut, bang)
    if (cut <= 1) then
       line = ''
    else
       line = line(:cut-1)
    end if
  end subroutine strip_config_comment

  subroutine normalize_config_whitespace(line)
    character(len=*), intent(inout) :: line

    integer :: i

    do i = 1, len(line)
       if (iachar(line(i:i)) == 9 .or. iachar(line(i:i)) == 13) line(i:i) = ' '
    end do
  end subroutine normalize_config_whitespace

  pure function lowercase(text) result(lowered)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lowered

    integer :: i, code

    lowered = text
    do i = 1, len(text)
       code = iachar(lowered(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) lowered(i:i) = achar(code + 32)
    end do
  end function lowercase

  subroutine invalid_value(filename, line_number, key, value, expected)
    character(len=*), intent(in) :: filename, key, value, expected
    integer, intent(in) :: line_number

    call config_error(filename, line_number, trim(key)//' expects a valid '//trim(expected)// &
         ', got: '//trim(value))
  end subroutine invalid_value

  subroutine config_error(filename, line_number, message)
    character(len=*), intent(in) :: filename, message
    integer, intent(in) :: line_number

    if (len_trim(filename) > 0 .and. line_number > 0) then
       write(error_unit,'(A,A,A,I0,A,A)') 'SISMO configuration error in ', trim(filename), &
            ':', line_number, ': ', trim(message)
    else if (len_trim(filename) > 0) then
       write(error_unit,'(A,A,A,A)') 'SISMO configuration error in ', trim(filename), ': ', trim(message)
    else
       write(error_unit,'(A,A)') 'SISMO configuration error: ', trim(message)
    end if
    stop 2, quiet=.true.
  end subroutine config_error

end module sismo_config_io
