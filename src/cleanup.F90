module cleanup 

contains 

subroutine run_cleanup(errtag, errcode)

use global 
use gll_library, only: cleanup_gll1d
use element
use model
use integration
use free_surface
#if (USE_MPI)
use mpi_library
#else
use serial_library
#endif
implicit none 

character(len=250) :: errtag ! error message
integer :: errcode



! Clean up and deallocate 
if(ISDISP_DOF)then
    deallocate(edofu)
  endif
  if(ISPOT_DOF)then
    deallocate(edofphi)
  endif
  ! clean up                                                                       
  call cleanup_model(errcode,errtag)
  call control_error(errcode,errtag,stdout,myrank)
  
  call cleanup_gll1d()
  
  call cleanup_hexface(errcode,errtag)                                             
  call control_error(errcode,errtag,stdout,myrank)
  
  call cleanup_integration(errcode,errtag)
  call control_error(errcode,errtag,stdout,myrank)
  
  call cleanup_integration2d(errcode,errtag)
  call control_error(errcode,errtag,stdout,myrank)
  
  call cleanup_free_surface()
  
  return 
end subroutine run_cleanup 
end module cleanup