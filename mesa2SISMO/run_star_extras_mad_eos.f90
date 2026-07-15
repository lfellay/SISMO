! Example MESA src/run_star_extras.f90 additions for MAD/nonadiabatic output.
!
! This file is not meant to replace an entire MESA work directory blindly.
! Copy the pieces below into your work directory's src/run_star_extras.f90,
! preserving any existing extras hooks you already use.
!
! The important part is:
!   s% how_many_extra_profile_columns => how_many_extra_profile_columns
!   s% data_for_extra_profile_columns => data_for_extra_profile_columns
!
! If this is the only extras customization in the work directory, also set
! warn_run_star_extras = .false. in &star_job.  Otherwise MESA may stop because
! the unused standard extras hooks are still assigned to null warning routines.
!
! The columns added are:
!   Cprho = (d ln Cp / d ln rho)_T
!   CpT   = (d ln Cp / d ln T)_rho
!   Qrho  = (d ln Q / d ln rho)_T, Q = chiT/chiRho
!   QT    = (d ln Q / d ln T)_rho, Q = chiT/chiRho

module run_star_extras

   use star_lib
   use star_def
   use const_def, only: dp
   use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, &
      i_Cp, i_chiRho, i_chiT
   use eos_lib, only: eosDT_get

   implicit none

contains

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! Keep any other hooks you already use in this routine, then add:
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns
   end subroutine extras_controls

   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: s

      how_many_extra_profile_columns = 0
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      how_many_extra_profile_columns = 4
   end function how_many_extra_profile_columns

   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character(len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz,n)
      integer, intent(out) :: ierr

      type(star_info), pointer :: s
      integer :: k
      real(dp) :: logRho, logT
      real(dp) :: res(num_eos_basic_results)
      real(dp) :: d_dlnd(num_eos_basic_results)
      real(dp) :: d_dlnT(num_eos_basic_results)
      real(dp), allocatable :: d_dxa(:,:)

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (n /= 4) stop 'data_for_extra_profile_columns: expected 4 MAD EOS columns'

      ! Do not add these names to profile_columns.list.  MESA appends them here.
      names(1) = 'Cprho'
      names(2) = 'CpT'
      names(3) = 'Qrho'
      names(4) = 'QT'

      allocate(d_dxa(num_eos_d_dxa_results, s% species))

      do k = 1, nz
         logRho = log10(s% rho(k))
         logT = log10(s% T(k))

         call eosDT_get( &
            s% eos_handle, s% species, s% chem_id, s% net_iso, &
            s% xa(1:s% species,k), s% rho(k), logRho, s% T(k), logT, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
         if (ierr /= 0) return

         vals(k,1) = safe_log_deriv(d_dlnd(i_Cp), res(i_Cp))
         vals(k,2) = safe_log_deriv(d_dlnT(i_Cp), res(i_Cp))
         vals(k,3) = safe_log_deriv(d_dlnd(i_chiT), res(i_chiT)) &
                   - safe_log_deriv(d_dlnd(i_chiRho), res(i_chiRho))
         vals(k,4) = safe_log_deriv(d_dlnT(i_chiT), res(i_chiT)) &
                   - safe_log_deriv(d_dlnT(i_chiRho), res(i_chiRho))
      end do

      deallocate(d_dxa)
   end subroutine data_for_extra_profile_columns

   real(dp) function safe_log_deriv(dval_dlnx, val)
      real(dp), intent(in) :: dval_dlnx, val

      if (abs(val) > tiny(1d0)) then
         safe_log_deriv = dval_dlnx/val
      else
         safe_log_deriv = 0d0
      end if
   end function safe_log_deriv

end module run_star_extras
