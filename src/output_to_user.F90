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