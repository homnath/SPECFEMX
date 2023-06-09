module output_to_user
    
contains
!_______________________________________________________________________________

    subroutine write_ifproc0(unit)
        ! Writes to log file only if processor rank = 0
        use global, only: logunit, myrank, log_msg

        implicit none 
        integer :: unit

        if(myrank==0)then
            write(unit,'(a)') trim(log_msg)
        endif
        flush(unit)

    end subroutine
!-------------------------------------------------------------------------------

    subroutine print_model_details
    use global

    if(myrank==0)then
        write(logunit,'(a)')'------------- MODEL DETAILS -------------'

        if(ISDISP_DOF)then
            write(logunit,'(a)')' * Displacement DOF: '//'ON'
            if(isplastic)then
                write(logunit,'(a)')' * Plasticity      : '//'ON'
            else
                write(logunit,'(a)')' * Plasticity      : '//'OFF'
            endif
        else
            write(logunit,'(a)')' * Displacement DOF: '//'OFF'
        endif
        write(logunit,*)


        if(ISPOT_DOF)then
            write(logunit,'(a)')' * Potential DOF   : '//'ON'
            write(logunit,'(a)')'   --> DOF type    : '//trim(POT_STRING)
        else
            write(logunit,'(a)')' * Potential DOF: '//'OFF'
        endif

        write(logunit,*)
        write(logunit,*)
        flush(logunit)

    endif
    
    ngllxy=ngllx*nglly                                                               
    ngllyz=nglly*ngllz                                                               
    ngllzx=ngllz*ngllx                                                               
    maxngll2d=max(ngllxy,ngllyz,ngllzx)
    
    end subroutine print_model_details 
!-------------------------------------------------------------------------------



    subroutine log_ksp_iteration(maxdu, ksp_iter, ksp_convreason)
        use global 
        use set_precision

        implicit none 
        real(kind=kreal) :: maxdu
        integer :: ksp_iter, ksp_convreason

        if(myrank==0)then
            write(logunit,'(a,i0,1x,a,g0.6)')' KSP iters: ',ksp_iter, &
            'max du: ',maxdu
            if(solver_type.eq.petsc_solver)then
              write(logunit,'(a,i0)')' convergence reason: ',ksp_convreason
            endif
            flush(logunit)
        endif

    end subroutine log_ksp_iteration
!-------------------------------------------------------------------------------


! Prints error if the memory can not be allocated
subroutine check_memory_alloc(istat, fileval)
    implicit none 
    integer :: istat, fileval

    if (istat/=0)then
        write(fileval,*)'ERROR: cannot allocate memory!'
        flush(fileval)
        stop
    endif
end subroutine check_memory_alloc
!-------------------------------------------------------------------------------


subroutine log_KSP_summary()
    
    use global
    use ksp_constants
    implicit none 
    ! open summary file
    if(myrank==0)then
        write(logunit,'(a)')'KSP_MAXITER, KSP_TOL, NL_MAXITER, NL_TOL'
        write(logunit,'(i0,1x,g0.6,1x,i0,1x,g0.6)')KSP_MAXITER,KSP_RTOL,NL_MAXITER,NL_TOL
        write(logunit,'(a,i0)')'Number of time steps:',nstep
        flush(logunit)
    endif
end subroutine log_KSP_summary
!-------------------------------------------------------------------------------

end module output_to_user
!===============================================================================
