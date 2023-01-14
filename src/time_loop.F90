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

 

end module