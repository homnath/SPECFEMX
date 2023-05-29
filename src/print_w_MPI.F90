module print_w_MPI

contains
    subroutine print_completion_details(cpu_tstart,cpu_tend,telap,&
    max_telap,mean_telap, format_str)
! USES
use global
#if(USE_MPI)
use mpi_library
use math_library_mpi
#else 
use serial_library
use math_library_serial
#endif 


implicit none 
! IO variables
real(kind=kreal) :: cpu_tstart,cpu_tend,telap,max_telap,mean_telap
character(len=20) ::format_str

! Local variables
logical :: isopen

! compute elapsed time
call cpu_time(cpu_tend)
telap=cpu_tend-cpu_tstart
max_telap=maxscal(telap)
mean_telap=sumscal(telap)/real(nproc,kreal)

if(myrank==0)then
write(format_str,*)ceiling(log10(real(max_telap)+1.))+5 ! 1 . and 4 decimals
format_str='(3(f'//trim(adjustl(format_str))//'.4,1X))'
write(logunit,'(a)')'ELAPSED TIME, MAX ELAPSED TIME, MEAN ELAPSED TIME'
write(logunit,fmt=format_str)telap,max_telap,mean_telap
write(logunit,'(a)')'--------------------------------------------'
flush(logunit)
close(logunit)
endif
!-----------------------------------

if(myrank==0)then
inquire(stdout,opened=isopen)
if(isopen)close(stdout)
endif


end subroutine print_completion_details
end module print_w_MPI