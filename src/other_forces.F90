! Holds all the bits for forces that aren't traction: 
module other_forces 
use shared
implicit none

contains
!_______________________________________________________________________________

subroutine compute_cmt_load(istep,freq)

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
integer,intent(in)           :: istep
real(kind=kreal)             :: freq
integer                      :: errcode
character(len=250)           :: errtag

! Note: initial istep=0 for the frequency domain and 1 for the time domain.
if( (steptype.eq.TIMESTEP .and. istep.eq.1) .or. &
    (steptype.eq.FREQSTEP .and. istep.eq.0) )then
  log_msg = trim(' Earthquake source type: moment-density tensor')
  call write_ifproc0
  ! First compute only the time/frequency independent factor, eqload0.
  call earthquake_load(neq,errcode,errtag)
  call sync_process
  call control_error(errcode,errtag,stdout)
endif
if(steptype==FREQSTEP)then
  ! Frequency domain
  !WARNING: make it general for nsrc
  extload=eqload0*source_frequency_function(freq,source_hdur(1))
else
  ! Time domain
  ! Coseismic or Postseismic. Load at the beginning only.
  if(istep.eq.1)then
    extload=extload+eqload0
  endif
endif
end subroutine compute_cmt_load
!===============================================================================    

subroutine compute_magnetic_traction(errcode, errtag)

use global
use output_to_user
!use mpi_library ! but what about serial version 
use mtraction

integer :: errcode
character(len=250) :: errtag

! apply magnetic traction
log_msg = trim('applying magnetic traction...') ;   call write_ifproc0
call apply_mtraction(errcode,errtag)
!call sync_process
call control_error(errcode,errtag,stdout)
if(myrank==0)then
  write(logunit,*)'complete!',maxval(abs(extload))
  flush(logunit)
endif

end subroutine compute_magnetic_traction
!===============================================================================    

subroutine compute_electrical_load(errcode, errtag)

use global
use output_to_user
!use mpi_library ! but what about serial version 
use electrical

integer :: errcode
character(len=250) :: errtag

log_msg = trim('computing electrical load...') ;   call write_ifproc0
call electrical_load(errcode,errtag)
!call sync_process
call control_error(errcode,errtag,stdout)
!print*,'in electrical:',maxval(abs(extload))
end subroutine compute_electrical_load
!===============================================================================    

subroutine compute_split_node_load(t, i_step, sfac, errcode, errtag)
    ! uses 
    use global ! itaper_slip, divide_slip, iseqsource, eqsource_type, srate
    use fault 
    use math_constants
 
    implicit none 
    ! IO Variables
    real(kind=kreal)              :: t
    integer                       :: i_step 
    real(kind=kreal)              :: sfac
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
          call compute_fault_slip_load(-1,sfac,errcode,errtag)
          ! Minus side
          call compute_fault_slip_load(1,sfac,errcode,errtag)
        else
          sfac=ONE
          ! Plus side
          call compute_fault_slip_load(1,sfac,errcode,errtag)
        endif
        ! The "sync" here is very important because some processors arrive this 
        ! stage faster than other. This may hang going to control_error routine!
        !call sync_process
        call control_error(errcode,errtag,stdout)
      endif

      if(srate)then
        extload=extload+t*slipload
      else
        extload=extload+slipload
      endif

    return 

end subroutine compute_split_node_load
!===============================================================================    
  
end module other_forces
!===============================================================================    
