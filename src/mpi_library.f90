! routines requied for MPI proecsses
! REVISION:
!   HNG, Jul 11,2011; Apr 09,2010
module my_mpi
! main parameter module for specfem simulations
use mpi
implicit none

! my MPI groups
integer :: my_local_mpi_comm_world
integer :: my_local_mpi_comm_for_bcast
integer :: my_local_mpi_comm_inter        ! MPI subgroup for hdf5 i/o server

end module my_mpi
!===============================================================================

module mpi_library
use mpi
contains
!-------------------------------------------------------------------------------
! start MPI processes
subroutine start_process()
use global,only:ismpi,myrank,nproc,stdout
use my_mpi,only: my_local_mpi_comm_world
implicit none
integer :: errcode
ismpi=.true. ! parallel

call MPI_INIT(errcode)
if(errcode /= 0) call mpierror('ERROR: cannot initialize MPI!',errcode,stdout)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,errcode)
if(errcode /= 0) call mpierror('ERROR: cannot find processor ID (rank)!',errcode,stdout)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,errcode)
if(errcode /= 0) call mpierror('ERROR: cannot find number of processors!',errcode,stdout)

! broadcast parameters read from main to all processes
my_local_mpi_comm_world = MPI_COMM_WORLD

if(myrank==0)then
    write(*,*)'* Running in parallel... '
    write(*,*)'  --> number of processors: ', nproc
  endif

return
end subroutine start_process
!-------------------------------------------------------------------------------

! close all MPI processes
subroutine close_process()
implicit none
integer :: errcode
call MPI_FINALIZE(errcode)
stop
return
end subroutine close_process
!-------------------------------------------------------------------------------

! MPI error
subroutine mpierror(errtag,errcode,stdout)
implicit none
character(len=*) :: errtag
integer :: errcode,stdout
write(stdout,'(a,a,i4,a)')errtag,' (MPI ERROR code:',errcode,')!'
stop
end subroutine mpierror
!-------------------------------------------------------------------------------

! syncronize all MPI processes
subroutine sync_process()
implicit none
integer :: errcode

call MPI_BARRIER(MPI_COMM_WORLD,errcode)

end subroutine sync_process
!-------------------------------------------------------------------------------

subroutine bcast_all_i(buffer, countval)

use my_mpi

implicit none

integer :: countval
integer, dimension(countval) :: buffer

integer :: ier

! checks if anything to do
if (countval == 0) return
call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_i
!-------------------------------------------------------------------------------

subroutine bcast_all_singlei(buffer)

use my_mpi

implicit none

integer :: buffer

integer :: ier

call MPI_BCAST(buffer,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_singlei
!-------------------------------------------------------------------------------

subroutine bcast_all_singlel(buffer)

use my_mpi

implicit none

logical :: buffer

integer :: ier

call MPI_BCAST(buffer,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_singlel
!-------------------------------------------------------------------------------

subroutine bcast_all_cr(buffer, countval)

use my_mpi
use set_precision, only: kreal
use set_precision_mpi, only: MPI_KREAL

implicit none

integer :: countval
real(kind=kreal), dimension(countval) :: buffer

integer :: ier

! checks if anything to do
if (countval == 0) return

call MPI_BCAST(buffer,countval,MPI_KREAL,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_cr
!-------------------------------------------------------------------------------

subroutine bcast_all_singlecr(buffer)

use my_mpi
use set_precision, only: kreal
use set_precision_mpi, only: MPI_KREAL

implicit none

real(kind=kreal) :: buffer

integer :: ier

call MPI_BCAST(buffer,1,MPI_KREAL,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_singlecr
!-------------------------------------------------------------------------------

subroutine bcast_all_r(buffer, countval)

use my_mpi

implicit none

integer :: countval
real, dimension(countval) :: buffer

integer :: ier

! checks if anything to do
if (countval == 0) return

call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_r
!-------------------------------------------------------------------------------

subroutine bcast_all_dp(buffer, countval)

use my_mpi

implicit none

integer :: countval
double precision, dimension(countval) :: buffer

integer :: ier

! checks if anything to do
if (countval == 0) return

call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_dp
!-------------------------------------------------------------------------------

subroutine bcast_all_singledp(buffer)

use my_mpi

implicit none

double precision :: buffer

integer :: ier

call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_singledp
!-------------------------------------------------------------------------------

subroutine bcast_all_ch_array(buffer,countval,STRING_LEN)

use my_mpi

implicit none

integer :: countval, STRING_LEN

character(len=STRING_LEN), dimension(countval) :: buffer

integer :: ier

! checks if anything to do
if (countval == 0) return

call MPI_BCAST(buffer,STRING_LEN*countval,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_ch_array
!-------------------------------------------------------------------------------

subroutine bcast_all_l_array(buffer, countval)

use my_mpi

implicit none

integer :: countval
logical, dimension(countval) :: buffer
integer :: ier

! checks if anything to do
if (countval == 0) return

call MPI_BCAST(buffer,countval,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

end subroutine bcast_all_l_array
!-------------------------------------------------------------------------------
!
! MPI gather helper
!
!-------------------------------------------------------------------------------

subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

use my_mpi

implicit none

integer :: sendcnt, recvcount, NPROC
integer, dimension(sendcnt) :: sendbuf
integer, dimension(recvcount,0:NPROC-1) :: recvbuf

integer :: ier

call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                recvbuf,recvcount,MPI_INTEGER, &
                0,my_local_mpi_comm_world,ier)

end subroutine gather_all_i
!-------------------------------------------------------------------------------

subroutine gather_all_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

use my_mpi

implicit none

integer :: sendcnt, recvcount, NPROC
integer, dimension(sendcnt) :: sendbuf
integer, dimension(recvcount,0:NPROC-1) :: recvbuf

integer :: ier

call MPI_ALLGATHER(sendbuf,sendcnt,MPI_INTEGER, &
                   recvbuf,recvcount,MPI_INTEGER, &
                   my_local_mpi_comm_world,ier)

end subroutine gather_all_all_i
!-------------------------------------------------------------------------------

subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

use my_mpi

implicit none

integer :: NPROC
integer :: sendbuf
integer, dimension(0:NPROC-1) :: recvbuf

integer :: ier

call MPI_GATHER(sendbuf,1,MPI_INTEGER, &
                recvbuf,1,MPI_INTEGER, &
                0,my_local_mpi_comm_world,ier)

end subroutine gather_all_singlei
!-------------------------------------------------------------------------------

subroutine gather_all_all_singlei(sendbuf, recvbuf, NPROC)

use my_mpi

implicit none

integer :: NPROC
integer :: sendbuf
integer, dimension(0:NPROC-1) :: recvbuf

integer :: ier

call MPI_ALLGATHER(sendbuf,1,MPI_INTEGER, &
                   recvbuf,1,MPI_INTEGER, &
                   my_local_mpi_comm_world,ier)

end subroutine gather_all_all_singlei
!-------------------------------------------------------------------------------

subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

use my_mpi
use set_precision, only: kreal
use set_precision_mpi, only: MPI_KREAL

implicit none

integer :: sendcnt, recvcount, NPROC
real(kind=kreal), dimension(sendcnt) :: sendbuf
real(kind=kreal), dimension(recvcount,0:NPROC-1) :: recvbuf

integer :: ier

call MPI_GATHER(sendbuf,sendcnt,MPI_KREAL, &
                recvbuf,recvcount,MPI_KREAL, &
                0,my_local_mpi_comm_world,ier)

end subroutine gather_all_cr
!-------------------------------------------------------------------------------

subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

use my_mpi

implicit none

integer :: sendcnt, recvcount, NPROC
double precision, dimension(sendcnt) :: sendbuf
double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

integer :: ier

call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
                recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
                0,my_local_mpi_comm_world,ier)

end subroutine gather_all_dp
!-------------------------------------------------------------------------------

end module mpi_library
!===============================================================================
