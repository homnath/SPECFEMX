module output_to_user
    
    contains 
    subroutine write_ifproc0()
        ! Writes to log file only if processor rank = 0
        use global, only: logunit, myrank, log_msg

        implicit none 

        if(myrank==0)then
            write(logunit,'(a)') trim(log_msg)
          endif
    end subroutine


end module output_to_user