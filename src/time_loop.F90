module time_loop 
use global 
use sea_level
use set_precision
implicit none 

contains 
!_______________________________________________________________________________

    subroutine calc_time_step(i_step, t, dt, freq, ang_freq, scale_ang_freq2)
        ! Calculates the time or frequency step for the time marching.
        ! USES
        use global 
        use set_precision
        use math_constants 
        use nondimensionpar
        ! IO variables
        real(kind=kreal) :: freq, ang_freq, scale_ang_freq2
        real(kind=kreal) :: t, dt
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
!-------------------------------------------------------------------------------
 

subroutine reset_nodal_arrays_loads(nodalslrate)
    use global
    use math_constants
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
!-------------------------------------------------------------------------------

subroutine set_elasto_visco_stiffness_matrix(i_step, dt, isscale_ang_freq, & 
                                             ang_freq, scale_ang_freq2,istep0)
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
integer :: i_step ,i , istep0
logical            :: isscale_ang_freq
real(kind=kreal)   :: ang_freq, scale_ang_freq2, dt

! Local: 
integer            :: errcode
character(len=250) :: errtag 
logical            :: reuse_pc_bool,freq_bool  

! Code: 
if(steptype.eq.FREQSTEP)then
    ! FREQUENCY STIFFNESS MATRIX 

    if(i_step==0)then
        call compute_stiffness_elastic(errcode,errtag)
    endif
        
    ! Set Petsc stiffness matrix
    if(solver_type.eq.petsc_solver)then
      reuse_pc_bool=.false.
      freq_bool=.true.
      call set_petsc_stiffness(isscale_ang_freq, &  
      ang_freq, scale_ang_freq2, reuse_pc_bool,freq_bool)  
    endif

else ! TIMESTEPPING

        ! If it is the first timestep we need the elastic stiffness matrix (storekmat)
        ! It should be constant so we dont need to change it
        if(i_step.eq.istep0)then
            if(myrank.eq.0.and.verbose_bool)then
                write(*,*)'Calculating elastic stiffness matrix...'
            endif
            call compute_stiffness_elastic(errcode,errtag)
            if(myrank.eq.0.and.verbose_bool)then
                write(*,*)' ✓ Done'
                write(*,*)
            endif
        endif 



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

        
        ! At all timesteps we need the SL contribution to the kmat: 
        ! Combine with the main storekmat
        ! At the end of the timestep we will then remove this so that
        ! we recover the original storekmat
        ! I dont think we will have the memory to store two kmats
        ! we may even need to get rid of the storekmatSL at some point
        ! And directly add/remove
        if (ISSL_DOF) then
            if(myrank.eq.0.and.verbose_bool)then 
                write(*,*)' --> Calculating sea-level stiffness matrix'
              endif 
            call compute_storekmatSL(storekmatSL, kSL)
            storekmat = storekmat + storekmatSL
            if(myrank.eq.0.and.verbose_bool)then
                write(*,*)' --> Combined SL and normal Kmats'
              endif
        endif
    
        ! Add the matrix to petsc: 
        if(myrank.eq.0.and.verbose_bool)then
            write(*,*)'Setting PETSC stiffness matrix...'
        endif

        if(solver_type.eq.petsc_solver)then
          reuse_pc_bool=.false.
          freq_bool=.false.
            call set_petsc_stiffness(isscale_ang_freq, &  
            ang_freq, scale_ang_freq2, reuse_pc_bool,freq_bool)   
        endif

        if(myrank.eq.0.and.verbose_bool)then
            write(*,*)' ✓ Done'
            write(*,*)
        endif

        ! DO NOT DELETE - HAS VISCOELASTIC CONTENT
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
!-------------------------------------------------------------------------------

subroutine run_convergence_loop(nodalsl, nodalslrate, nodalice, nodalicerate)
    use global 
    use set_precision
#if (USE_MPI)
use math_library_mpi
use mpi
#else
use math_library_serial
#endif

    ! IO variables: 
    real(kind=kreal), allocatable :: nodalsl(:)
    real(kind=kreal), allocatable :: nodalslrate(:), nodalice(:), nodalicerate(:)

    ! Local variables: 
    real(kind=kreal)   :: mass_imbalance_0, dt_cloop, mass_imbalance
    integer            :: i_cloop
    character(len=250) :: errtag 
    integer            :: errcode
    errtag=""; errcode=-1

    ! Initialise
    cloop_converged = .false.
    dt_cloop = 1.0_kreal

    ! Run a loop to find the best update (dt) that converges towards mass balance
    cloop: do i_cloop=1,ncloop_MAX
      ! WE: Update our vectors with a timestep: 
      ! WE: For sea level we solve for the time derivatives of the system so then we use
      ! WE: u(t + dt) = u(t) + dt * f(t) where f is the time derivative
      ! WE: These time derivatives are stored in nodalu, nodalphi, nodalsl 
      nodalustore   = nodalustore   +  dt_cloop*nodalu
      nodalphistore = nodalphistore +  dt_cloop*nodalphi
      nodalsl       = nodalsl       + (dt_cloop*nodalslrate)
      nodalice      = nodalice      +  nodalicerate

      ! Update ocean function and ocean area/volume 
      call update_ocean_function(nodalice, nodalsl, errcode, errtag, verbose=.false.)  

      ! Note that for first timestep we need to update the area of the ocean
      if(i_cloop.eq.1)then
        call update_SL_area(nodalsl, overwrite_old=.true., verbose=.false.)
      else 
        call update_SL_area(nodalsl, overwrite_old=.false., verbose=.false.)
      endif 

      ! Update the mass imbalance
      mass_imbalance = SLsummasschange + icemasschange_per_ts


      ! Check convergence 
      if (ABS(mass_imbalance).lt.CLOOP_CONV_THRESH)then
        cloop_converged = .true.
      endif 


      ! store initial value for reference, avoiding division by zero
      if(i_cloop.eq.1)then 
        if(icemasschange_per_ts.ne.zero)then 
          mass_imbalance_0 = (mass_imbalance/ABS(icemasschange_per_ts))*100
        else
          mass_imbalance_0=zero
        endif 
      endif 


      ! Evaluate the ocean nodes as proportion of overall FS nodes
      totaloceannodes = sumscal(oceannodes);

      ! Output to user if final timestep 
      if(cloop_converged.or.i_cloop.eq.ncloop_MAX)then
        if(myrank.eq.0)then
            write(*,*)
            write(*,*)'                    COMPLETED CLOOPING:'
            write(*,*)'------------------------------------------------------------'
            if(i_cloop.eq.ncloop_MAX)then 
                write(*,*)' WARNING - DIDNT CONVERGE WITHIN THRESHOLD: ', CLOOP_CONV_THRESH
            endif 
            write(*,'(a,i0)')' Number of loops              : ', i_cloop
            write(*,'(a,g0.4)')' Timestep                     : ',  dt_cloop
            write(*,'(a,i0, a,i0)')'  Total ocean nodes           : ', totaloceannodes,'/', allnodesfs
            write(*,'(a,g0.4)')' Mass of ice (kg)             : ',  icemasschange_per_ts
            write(*,'(a,g0.4)')' Mass of water (kg)           : ',  SLsummasschange
            write(*,'(a,g0.4)')' Mass imbalance (kg)          : ',  mass_imbalance
            write(*,'(a,g0.4)')'  Original mass imbalance (%) : ', mass_imbalance_0

            if (icemasschange_per_ts.ne.zero)then
            write(*,'(a,g0.4)')'  Mass imbalance          (%) : ', (mass_imbalance/ABS(icemasschange_per_ts))*100
            else
            write(*,'(a,g0.4)')'  Mass imbalance          (%) : 0 since no ice change'
            endif 
            write(*,*)
        endif

        exit cloop

      ! Recover original values for next cloop unless final loopstep
      else
        nodalustore   = nodalustore   -  dt_cloop*nodalu
        nodalphistore = nodalphistore -  dt_cloop*nodalphi
        nodalsl       = nodalsl       - (dt_cloop*nodalslrate)
        nodalice      = nodalice      -  nodalicerate

        ! Update the dt for next timestep
        if (icemasschange_per_ts.ne.zero)then
          dt_cloop = dt_cloop + mass_imbalance/icemasschange_per_ts
        endif
      endif

    enddo cloop

end subroutine run_convergence_loop
!-------------------------------------------------------------------------------

end module
