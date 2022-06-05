! Holds all the bits for forces that arent traction: 


module other_forces 

    contains
    subroutine compute_magnetic_traction(errcode, errtag, extload)

    use global
    use output_to_user
    use mpi_library ! but what about serial version 
    use mtraction

    integer :: errcode
    character(len=250) :: errtag
    real(kind=kreal),allocatable ::extload(:)



    ! apply magnetic traction
    if(ismtraction)then
        log_msg = trim('applying magnetic traction...') ;   call write_ifproc0()
        call apply_mtraction(extload,errcode,errtag)
        call sync_process
        call control_error(errcode,errtag,stdout,myrank)
        if(myrank==0)then
        write(logunit,*)'complete!',maxval(abs(extload))
        flush(logunit)
        endif
    endif


    end subroutine compute_magnetic_traction

end module other_forces