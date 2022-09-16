module stiffness_matrix 

    implicit none 
    contains 


subroutine get_stiffness_matrix_freq()

use global 
implicit none

! compute elastic stiffness matrix for time = 0
    if(i_step==0)then
        call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
      endif
      
      ! Set Petsc stiffness matrix
      if(solver_type.eq.petsc_solver)then
        if (ISSL_DOF)then 
          call set_petsc_stiffness_SL(isscale_ang_freq, storekmat, QSL,& 
          storeRu, storeRphi, storemmat, ang_freq, scale_ang_freq2,    &
          reuse_pc_bool=.false.,freq_bool=.true.)   
        else 
          call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
          ang_freq, scale_ang_freq2, reuse_pc_bool=.false.,freq_bool=.true.)   
        endif 
      endif

end subroutine get_stiffness_matrix_freq



subroutine get_stiffness_matrix_elastic


if(i_step==1)then 
    ! compute elastic stiffness matrix for time = 0
    call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
    
    if(solver_type.eq.petsc_solver)then
        if (ISSL_DOF)then 
        call set_petsc_stiffness_SL(isscale_ang_freq, storekmat, QSL,& 
        storeRu, storeRphi, storemmat, ang_freq, scale_ang_freq2,    &
        reuse_pc_bool=.false.,freq_bool=.false.)   
        else 
        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
        ang_freq, scale_ang_freq2, reuse_pc_bool=.false.,freq_bool=.false.)   
        endif 
    endif


end subroutine get_stiffness_matrix_elastic




end module stiffness_matrix