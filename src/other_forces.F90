! Holds all the bits for forces that arent traction: 


module other_forces 

implicit none

contains




subroutine compute_cmt_load(extload, freq)

use global
use output_to_user
use earthquake
use source_function
#if(USE_MPI)
use mpi_library
use math_library_mpi
use ghost_library_mpi
#else
use serial_library
use math_library_serial
#endif
use cmtsolution,only:source_tshift,source_hdur
    
implicit none 

real(kind=kreal),allocatable :: extload(:)
real(kind=kreal)             :: freq
integer                      :: errcode
character(len=250)           :: errtag

log_msg = trim(' Earthquake source type: moment-density tensor')
call write_ifproc0(logunit)

call earthquake_load(neq,extload,errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)

if(steptype==FREQSTEP)then
  !WARNING: make it general for nsrc
  extload=extload*source_frequency_function_complex(freq,source_hdur(1))
endif
end subroutine compute_cmt_load
    
  
  


subroutine compute_magnetic_traction(errcode, errtag, extload)

use global
use output_to_user
use mpi_library ! but what about serial version 
use mtraction

integer :: errcode
character(len=250) :: errtag
real(kind=kreal),allocatable ::extload(:)



! apply magnetic traction
log_msg = trim('applying magnetic traction...') ;   call write_ifproc0(logunit)
call apply_mtraction(extload,errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)
if(myrank==0)then
  write(logunit,*)'complete!',maxval(abs(extload))
  flush(logunit)
endif


end subroutine compute_magnetic_traction


subroutine compute_split_node_load(t, i_step, sfac, slipload, extload, storekmat, errcode, errtag)
    ! uses 
    use global ! itaper_slip, divide_slip, iseqsource, eqsource_type, srate
    use fault 
    use math_constants
 
    ! SHOULD BE SERIAL LIBRARY.F90 if running serial compilation 
    use mpi_library

    implicit none 
    ! IO Variables
    real(kind=kreal)              :: t
    integer                       :: i_step 
    real(kind=kreal)              :: sfac
    real(kind=kreal), allocatable :: slipload(:),extload(:)
    real(kind=kreal), allocatable :: storekmat(:,:,:)
    character(len=250) :: errtag ! error message
    integer :: errcode

    ! Local variables: NONE

  ! compute load contributed by the earthquake slip
  ! split-node apparoch: prescribe the slip on the fault explicitly
      if(i_step==1)then
        if(myrank==0)then
          write(logunit,'(a)')'  Earthquake source type: slip with split node'
          write(logunit,'(a,1x,i2)')'  Slip taper option: ',itaper_slip
          flush(logunit)
        endif
        if(divide_slip)then
          sfac=HALF
          ! Plus side
          call compute_fault_slip_load(-1,sfac,storekmat,slipload,errcode,errtag)
          ! Minus side
          call compute_fault_slip_load(1,sfac,storekmat,slipload,errcode,errtag)
        else
          sfac=ONE
          ! Plus side
          call compute_fault_slip_load(1,sfac,storekmat,slipload,errcode,errtag)
        endif
        ! The "sync" here is very important because some processors arrive this 
        ! stage faster than other. This may hang going to control_error routine!
        call sync_process
        call control_error(errcode,errtag,stdout,myrank)
      endif

      if(srate)then
        extload=t*slipload
      else
        extload=slipload
      endif

    return 

end subroutine compute_split_node_load




  
end module other_forces
