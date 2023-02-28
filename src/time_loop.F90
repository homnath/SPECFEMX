module time_loop 
    use global 
    use sea_level
    use set_precision
    implicit none 

    contains 


    subroutine calc_time_step(i_step, t, dt, freq, ang_freq, scale_ang_freq2)
        ! Calculates the time or frequency step for the time marching.
        ! USES
        use global 
        use set_precision
        use math_constants 
        use dimensionless
        ! IO variables
        real(kind=kreal) :: freq, ang_freq, scale_ang_freq2
        real(kind=kreal) :: t, dt, step
        integer :: i_step

        step=step0+dstep*real(i_step,kreal)

        if(steptype.eq.TIMESTEP)then
            ! Time step.
            t=step
            dt=dstep
            if(myrank==0)then
            write(logunit,'(a,i0,a,g0.6)')'step: ',i_step,' t: ',t
            flush(logunit)
            endif

        
        elseif(steptype.eq.FREQSTEP)then
            ! Frequency step.
            freq=step
            if(myrank==0)then
            write(logunit,'(a,i0,a,g0.6)')'step: ',i_step,' f: ',freq
            flush(logunit)
            endif

            if(devel_nondim)then
            ang_freq=TWO*freq*DIM_T
            else
            ang_freq=TWO*PI*freq
            endif
            scale_ang_freq2=ONE/(ang_freq*ang_freq)
        endif
        
    end subroutine calc_time_step

 

subroutine reset_nodal_arrays_loads(nodalu,ubcload,rhoload,nodalphi,nodalslrate)
    use global
    use math_constants
    real(kind=kreal),allocatable :: nodalu(:,:)
    real(kind=kreal),allocatable :: ubcload(:)
    real(kind=kreal),allocatable :: rhoload(:)
    real(kind=kreal),allocatable :: nodalphi(:)
    real(kind=kreal),allocatable :: nodalslrate(:)

    nodalu  = ZERO
    ubcload = ZERO
    rhoload = ZERO
  
    if(ISPOT_DOF)then
      nodalphi=ZERO
    endif
  
    if(ISSL_DOF)then
      nodalslrate = ZERO 
    endif
end subroutine reset_nodal_arrays_loads




subroutine set_elasto_visco_stiffness_matrix(i_step, dt, storekmat, storemmat, rhoload, isscale_ang_freq, & 
                                                ang_freq, scale_ang_freq2, nelmt_viscoelas, & 
                                                eid_viscoelas, relaxtime)
use global 
use matrix_vector
#if (USE_MPI)
use parsolver
use mpi_library
#if (USE_COMPLEX)
use parsolver_petsc_complex
#else
use parsolver_petsc
#endif
#else
use sparse_serial
use serial_library
use solver_petsc
#endif
        
! IO variables: 
integer :: i_step ,i 
real(kind=kreal), allocatable :: storekmat(:,:,:), storemmat(:,:),  rhoload(:)
logical            :: isscale_ang_freq
real(kind=kreal)   :: ang_freq, scale_ang_freq2, dt
real(kind=kreal), allocatable :: relaxtime(:,:) 
integer :: nelmt_viscoelas
integer,allocatable :: eid_viscoelas(:)

! Local: 
integer            :: errcode
character(len=250) :: errtag 


! Code: 
if(steptype.eq.FREQSTEP)then
    ! FREQUENCY STIFFNESS MATRIX 

    if(i_step==0)then
        call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
    endif
        
    ! Set Petsc stiffness matrix
    if(solver_type.eq.petsc_solver)then
        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
        ang_freq, scale_ang_freq2, reuse_pc_bool=.false.,freq_bool=.true.)  
    endif

    else ! TIMESTEPPING

        if(i_step.eq.1.or.i_step.eq.nstep)then
            if(myrank.eq.0)then
                write(*,*)
                write(*,*)'CALCULATING KMAT AT EVERY TIMESTEP!'
                write(*,*)
            endif
        endif 



        if(myrank.eq.0)then
            write(*,*)'Calculating elastic stiffness matrix...'
        endif
                
        call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
        if(myrank.eq.0)then
            write(*,*)' ✓ Done'
            write(*,*)
        endif
        
        if(myrank.eq.0)then
            write(*,*)'Setting PETSC stiffness matrix...'
        endif

        if(solver_type.eq.petsc_solver)then
            call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
            ang_freq, scale_ang_freq2, reuse_pc_bool=.false.,freq_bool=.false.)   
        endif

        if(myrank.eq.0)then
            write(*,*)' ✓ Done'
            write(*,*)
        endif
        





        ! ORIGINAL VERSION - BUT FOR SL WE NEED ADAPTIVE KMAT 
    !if(i_step==1)then 
    !    if(myrank.eq.0)then
    !        write(*,*)'Calculating elastic stiffness matrix...'
    !    endif
    !    call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
    !    if(myrank.eq.0)then
    !        write(*,*)' ✓ Done'
    !        write(*,*)
    !    endif

    !    if(myrank.eq.0)then
    !        write(*,*)'Setting PETSC stiffness matrix...'
    !    endif

    !    if(solver_type.eq.petsc_solver)then
    !        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
    !        ang_freq, scale_ang_freq2, reuse_pc_bool=.false.,freq_bool=.false.)   
    !    endif

    !    if(myrank.eq.0)then
    !        write(*,*)' ✓ Done'
    !        write(*,*)
    !    endif

    !elseif(i_step==2)then
        ! Since we use a uniform dt, following routine has to be called only once 
        ! for a linear viscoelastic model. For nonlinear or nonuniform time steps
        ! it has to be called for every time steps or every changing time step.
        ! This will simply overwrite the storekmat for viscoelastic elements.
    !    call compute_stiffness_viscoelastic(nelmt_viscoelas,             &   
    !                                        eid_viscoelas, dt, relaxtime,&
    !                                        storekmat, errcode, errtag)
    
        ! If using PETSC solver set stiffness matric                                      
    !    if(solver_type.eq.petsc_solver)then

    !        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
    !        ang_freq, scale_ang_freq2, reuse_pc_bool=.true.,freq_bool=.false.)   

    !    endif ! Petsc solver
    !endif ! istep
        
        



endif ! if(steptype.eq.FREQSTEP)

end subroutine set_elasto_visco_stiffness_matrix



end module