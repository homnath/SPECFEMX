module prepare_solver_mod
contains 

subroutine prepare_petsc_solver()
! USES: 
use global 
use output_to_user

#if(USE_MPI)
use prepare_sparse

#if(USE_COMPLEX)
use parsolver_petsc_complex
#else
use parsolver_petsc
#endif

#else
use prepare_sparse_serial
#endif



    ! IO Variables:

    ! Local variables: 


    ! Code: 
    ! Prepare sparsity of the stiffness matrix (working out size etc)
        call prepare_sparse() 
    ! petsc solver
        call petsc_initialize() 
        log_msg = 'petsc_initialize: SUCCESS!' ; call write_ifproc0()
    
    
    ! Create sparse vector, matrix, and preallocate                                                         
    ! TODO: following call is not necessary for RECYCLE                            
        call petsc_create_vector()                                                     
        call petsc_matrix_preallocate_size()                                           
        call petsc_create_matrix()                                                     
        call petsc_create_solver()                                                     
        log_msg = 'petsc_preallocate_matrix_size: SUCCESS!' ; call write_ifproc0()
        
        

end subroutine prepare_petsc_solver

end module prepare_solver_mod 
