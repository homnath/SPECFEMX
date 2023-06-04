module cleanup 

contains 
!_______________________________________________________________________________

subroutine run_cleanup_specfemx(errtag, errcode)

use global
use local 
use shared
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
  call control_error(errcode,errtag,stdout)
  
  call cleanup_gll1d()
  
  call cleanup_hexface(errcode,errtag)                                             
  call control_error(errcode,errtag,stdout)
  
  call cleanup_integration(errcode,errtag)
  call control_error(errcode,errtag,stdout)
  
  call cleanup_integration2d(errcode,errtag)
  call control_error(errcode,errtag,stdout)
  
  call cleanup_free_surface()
  
  return 
end subroutine run_cleanup_specfemx 
!-------------------------------------------------------------------------------

subroutine run_cleanup_specfem3d()
! USES
use global
use local
use fault
#if (USE_MPI)
use mpi_library
use ghost_library_mpi
use math_library_mpi
use sparse
use parsolver
#if (USE_COMPLEX)
use parsolver_petsc_complex
#else
use parsolver_petsc
#endif
#else
use serial_library
use math_library_serial
use sparse_serial
use solver
use solver_petsc
#endif

implicit none 
! IO variables

  !CODE:

  ! cleanup solver
  if(savedata%strain)then
    close(77)
    deallocate(strain_elmt,strain_nodal)
  endif


  if(solver_type.eq.petsc_solver)then
    call petsc_destroy_vector()                                                      
    call petsc_destroy_matrix()                                                      
    call petsc_destroy_solver()                                                      
    call petsc_finalize()
  endif
  
  call cleanup_fault()
  !deallocate(egdof,egdofu)
  !if(allocated(gdofu))deallocate(gdofu)
  deallocate(extload,load,resload,rhoload,ubcload)
  deallocate(du,u)
  deallocate(nodalu,bcnodalv)
  if(ISPOT_DOF)then
    deallocate(nodalphi)
  endif
  
  if(isplastic)then
    deallocate(olddu,evpt)
  endif
  if(solver_type.eq.builtin_solver .and.solver_diagscale)then
    deallocate(dprecon,ndscale)
  endif
  deallocate(mat_id,mat_domain,gam_blk,ym_blk,coh_blk,nu_blk,phi_blk,psi_blk,srf)
  if(allocated(imat_to_imatve))deallocate(imat_to_imatve)
  if(allocated(imatve_to_imat))deallocate(imatve_to_imat)
  deallocate(g_coord,g_num)
  deallocate(node_valency)
  deallocate(bmat,deriv,eld,num)
  deallocate(kmat)
  if(allocated(infinite_iface))deallocate(infinite_iface)
  if(allocated(infinite_face_idir))deallocate(infinite_face_idir)
  call cleanup_ghost()

end subroutine run_cleanup_specfem3d
!-------------------------------------------------------------------------------

end module cleanup
!===============================================================================
