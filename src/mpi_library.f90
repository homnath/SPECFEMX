! routines requied for MPI proecsses
! REVISION:
!   HNG, Jul 11,2011; Apr 09,2010
module mpi_library
use mpi
contains

! start MPI processes
subroutine start_process()
use global,only:ismpi,myrank,nproc,stdout
implicit none
integer :: errcode
ismpi=.true. ! parallel

call MPI_INIT(errcode)
if(errcode /= 0) call mpierror('ERROR: cannot initialize MPI!',errcode,stdout)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,errcode)
if(errcode /= 0) call mpierror('ERROR: cannot find processor ID (rank)!',errcode,stdout)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,errcode)
if(errcode /= 0) call mpierror('ERROR: cannot find number of processors!',errcode,stdout)

if(myrank==0)then
    write(*,*)'* Running in parallel... '
    write(*,*)'  --> number of processors: ', nproc
  endif

return
end subroutine start_process
!=======================================================

! close all MPI processes
subroutine close_process()
implicit none
integer :: errcode
call MPI_FINALIZE(errcode)
stop
return
end subroutine close_process
!=======================================================

! MPI error
subroutine mpierror(errtag,errcode,stdout)
implicit none
character(len=*) :: errtag
integer :: errcode,stdout
write(stdout,'(a,a,i4,a)')errtag,' (MPI ERROR code:',errcode,')!'
stop
end subroutine mpierror
!=======================================================

! syncronize all MPI processes
subroutine sync_process()
implicit none
integer :: errcode

call MPI_BARRIER(MPI_COMM_WORLD,errcode)

end subroutine sync_process
!=======================================================

end module mpi_library
!===============================================================================
