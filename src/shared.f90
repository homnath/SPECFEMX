! REVISION:
!   HNG, June 02,2023
module shared
contains
!_______________________________________________________________________________
subroutine check_allocate(ierr,errsrc)
implicit none
integer,intent(in) :: ierr
character(len=500),intent(in) :: errsrc
if(ierr.ne.0)then
    write(*,*)'ERROR: cannot allocate array/s!'
    write(*,*)'Code: ',ierr
    write(*,*)'Source: '//trim(errsrc)
    stop
endif
end subroutine check_allocate
!-------------------------------------------------------------------------------

! get processor tag
function proc_tag() result(ptag)
use global,only:ismpi,myrank,nproc
implicit none
character(len=20) :: format_str,ptag

if (ismpi) then
  write(format_str,*)ceiling(log10(real(nproc)+1.))
  format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

  write(ptag,fmt=format_str)'_proc',myrank
else
  ptag=''
endif
return
end function
!-------------------------------------------------------------------------------

! write error and stop
subroutine control_error(errcode,errtag,stdout)
implicit none
integer,intent(in) :: errcode
character(len=*),intent(in) :: errtag
integer,intent(in) :: stdout
integer :: ierr

if(errcode.eq.0)return
! print error message and stop execution
write(stdout,'(a)')trim(errtag)
flush(stdout)
stop
end subroutine control_error
!=======================================================
end module shared
!===============================================================================
