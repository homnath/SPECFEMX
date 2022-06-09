module cleanup 

contains 

subroutine run_cleanup()

use global 
use gll_library, only: cleanup_gll1d
implicit none 

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
  

end subroutine run_cleanup 