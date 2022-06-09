module output_to_user
    
    contains 
    subroutine write_ifproc0()
        ! Writes to log file only if processor rank = 0
        use global, only: logunit, myrank, log_msg

        implicit none 

        if(myrank==0)then
            write(logunit,'(a)') trim(log_msg)
        endif
        flush(logunit)

    end subroutine



    subroutine print_completion_details(cpu_tstart,cpu_tend,telap,max_telap,mean_telap)
        ! USES
        use global

        implicit none 
        ! IO variables
        real(kind=kreal) :: cpu_tstart,cpu_tend,telap,max_telap,mean_telap

        ! Local variables
        ! compute elapsed time
        call cpu_time(cpu_tend)
        telap=cpu_tend-cpu_tstart
        max_telap=maxscal(telap)
        mean_telap=sumscal(telap)/real(nproc,kreal)

        if(myrank==0)then
        write(format_str,*)ceiling(log10(real(max_telap)+1.))+5 ! 1 . and 4 decimals
        format_str='(3(f'//trim(adjustl(format_str))//'.4,1X))'
        write(logunit,'(a)')'ELAPSED TIME, MAX ELAPSED TIME, MEAN ELAPSED TIME'
        write(logunit,fmt=format_str)telap,max_telap,mean_telap
        write(logunit,'(a)')'--------------------------------------------'
        flush(logunit)
        close(logunit)
        endif
        !-----------------------------------

        if(myrank==0)then
        inquire(stdout,opened=isopen)
        if(isopen)close(stdout)
        endif


    subroutine print_model_details()

    use global

    if(myrank==0)then
        write(logunit,'(a)')'-----------------------------------------------'
        if(ISDISP_DOF)then
        write(logunit,'(a)')'Displacement DOF: '//'ON'
        write(logunit,'(a,i0)')' NNDOFU: ',nndofu
        if(isplastic)then
            write(logunit,'(a)')' Plasticity: ON'
        else
            write(logunit,'(a)')' Plasticity: OFF'
        endif
        else
        write(logunit,'(a)')'Displacement DOF: '//'OFF'
        write(logunit,'(a,i0)')' NNDOFU: ',nndofu
        endif
        if(ISPOT_DOF)then
        write(logunit,'(a)')'Potential DOF: '//'ON'
        write(logunit,'(a,i0)')' NNDOPHI: ',nndofphi
        write(logunit,'(a)')' Potential DOF type: '//trim(POT_STRING)
        else
        write(logunit,'(a)')'Potential DOF: '//'OFF'
        write(logunit,'(a,i0)')' NNDOPHI: ',nndofphi
        endif
        write(logunit,'(a,i0)')'Total DOFs per node: ',nndof
        write(logunit,'(a,3(i0,1x))')'Displacement DOF indices: ',idofu
        write(logunit,'(a,i0)')'Potential DOF indices: ',idofphi
        write(logunit,'(a,i0)')'NEDOFU: ',nedofu
        write(logunit,'(a,i0)')'NEDOFPHI: ',nedofphi
        write(logunit,'(a,i0)')'Total DOFs per element: ',nedof
        write(logunit,'(a)')'-----------------------------------------------'
        flush(logunit)
    endif
    
    write(logunit,*)'Sorting 2D mesh stuff...'
    ngllxy=ngllx*nglly                                                               
    ngllyz=nglly*ngllz                                                               
    ngllzx=ngllz*ngllx                                                               
    maxngll2d=max(ngllxy,ngllyz,ngllzx)
    
    write(logunit,*)'   ngllxy: ', ngllxy
    write(logunit,*)'   ngllyz: ', ngllyz
    write(logunit,*)'   ngllzx: ', ngllzx
    write(logunit,*)'   max GLL 2D: ', maxngll2d
    
  
    end subroutine print_model_details 

end module output_to_user