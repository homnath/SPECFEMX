!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


!-------------------------------------------------------------------------------
!  locate MPI slice which contains the point and bcast to all
!-------------------------------------------------------------------------------
module locate_mpi_domain
implicit none
contains
!-------------------------------------------------------------------------------
subroutine locate_MPI_slice(npoint_subset,npoint_already_done, &
                            ielmt_selected_subset, &
                            x_found_subset, y_found_subset, z_found_subset, &
                            xi_subset,eta_subset,gamma_subset, &
                            final_distance_subset, &
                            npoints_total, ielmt_selected, islice_selected, &
                            x_found,y_found,z_found, &
                            xi_point, eta_point, gamma_point, &
                            final_distance)

! locates subset of points in all slices

use math_constants,only: DHUGEVAL
use global,only: NDIM,NPROC,myrank
use mpi_library
use, intrinsic :: ieee_arithmetic
integer, intent(in) :: npoint_subset,npoint_already_done
integer, dimension(npoint_subset), intent(in)  :: ielmt_selected_subset
double precision, dimension(npoint_subset), intent(in) :: x_found_subset, y_found_subset, z_found_subset
double precision, dimension(npoint_subset), intent(in) :: xi_subset, eta_subset, gamma_subset
double precision, dimension(npoint_subset), intent(in) :: final_distance_subset

integer, intent(in) :: npoints_total
integer, dimension(npoints_total), intent(inout)  :: ielmt_selected, islice_selected
double precision, dimension(npoints_total), intent(inout)  :: x_found, y_found, z_found
double precision, dimension(npoints_total), intent(inout)  :: xi_point, eta_point, gamma_point
double precision, dimension(npoints_total), intent(inout)  :: final_distance

! local parameters
integer :: ipoint,ipoint_in_this_subset,iproc
double precision :: distmin

! gather arrays
integer, dimension(npoint_subset,0:NPROC-1) :: ielmt_selected_all
double precision, dimension(npoint_subset,0:NPROC-1) :: xi_all,eta_all,gamma_all
double precision, dimension(npoint_subset,0:NPROC-1) :: x_found_all,y_found_all,z_found_all
double precision, dimension(npoint_subset,0:NPROC-1) :: final_distance_all
! initializes with dummy values
ielmt_selected_all(:,:) = -1
xi_all(:,:) = 0.d0
eta_all(:,:) = 0.d0
gamma_all(:,:) = 0.d0
x_found_all(:,:) = 0.d0
y_found_all(:,:) = 0.d0
z_found_all(:,:) = 0.d0
final_distance_all(:,:) = DHUGEVAL

! gather all (on main process)
call gather_all_i(ielmt_selected_subset,npoint_subset,ielmt_selected_all,npoint_subset,NPROC)
call gather_all_dp(x_found_subset,npoint_subset,x_found_all,npoint_subset,NPROC)
call gather_all_dp(y_found_subset,npoint_subset,y_found_all,npoint_subset,NPROC)
call gather_all_dp(z_found_subset,npoint_subset,z_found_all,npoint_subset,NPROC)
call gather_all_dp(xi_subset,npoint_subset,xi_all,npoint_subset,NPROC)
call gather_all_dp(eta_subset,npoint_subset,eta_all,npoint_subset,NPROC)
call gather_all_dp(gamma_subset,npoint_subset,gamma_all,npoint_subset,NPROC)
call gather_all_dp(final_distance_subset,npoint_subset,final_distance_all,npoint_subset,NPROC)
! find the slice and element to put the source
if (myrank == 0) then

  ! loops over subset
  do ipoint_in_this_subset = 1,npoint_subset

    ! mapping from station/source number in current subset to real station/source number in all the subsets
    ipoint = ipoint_in_this_subset + npoint_already_done
    distmin = DHUGEVAL
    do iproc = 0,NPROC-1
      if (final_distance_all(ipoint_in_this_subset,iproc) < distmin) then
        distmin =  final_distance_all(ipoint_in_this_subset,iproc)

        islice_selected(ipoint) = iproc
        ielmt_selected(ipoint) = ielmt_selected_all(ipoint_in_this_subset,iproc)

        xi_point(ipoint)    = xi_all(ipoint_in_this_subset,iproc)
        eta_point(ipoint)   = eta_all(ipoint_in_this_subset,iproc)
        gamma_point(ipoint) = gamma_all(ipoint_in_this_subset,iproc)

        x_found(ipoint) = x_found_all(ipoint_in_this_subset,iproc)
        y_found(ipoint) = y_found_all(ipoint_in_this_subset,iproc)
        z_found(ipoint) = z_found_all(ipoint_in_this_subset,iproc)
      endif
    enddo
    final_distance(ipoint) = distmin
  enddo

endif ! end of section executed by main process only

end subroutine locate_MPI_slice
!-------------------------------------------------------------------------------

end module locate_mpi_domain
!===============================================================================
