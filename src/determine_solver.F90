module det_solver 
contains
    subroutine determine_solver(errcode, errtag)

        ! USES 
use global
use shared
#if(USE_MPI)
use mpi_library
use math_library_mpi
#else 
use serial_library
use math_library_serial
#endif 
        
        implicit none 
        character(len=250) :: errtag ! error message
        integer :: errcode
        
        ! solver type
        if(solver_type.eq.smart_solver)then
          if(ismpi)then
            solver_type=petsc_solver
          else
            solver_type=builtin_solver
          endif
        endif
        if(ismpi)then
          if(solver_type.eq.builtin_solver)then
            if(myrank==0)write(logunit,'(a)')'* Solver type: builtin parallel solver'
          elseif(solver_type.eq.petsc_solver)then
            if(myrank==0)write(logunit,'(a)')'* Solver type: PETSc parallel solver'
          else
            write(errtag,'(a,i4)')'ERROR: invalid solver type:',solver_type
            call control_error(errcode,errtag,logunit)
          endif
        else
          if(solver_type.eq.builtin_solver)then
            write(logunit,'(a)')'* Solver type: builtin serial solver'
          elseif(solver_type.eq.petsc_solver)then
            write(errtag,'(a,i4)')'ERROR: PETSc solver must be run with MPI:',solver_type
            call control_error(errcode,errtag,logunit)
          else
            write(errtag,'(a,i4)')'ERROR: invalid solver type:',solver_type
            call control_error(errcode,errtag,logunit)
          endif
        endif
        
        
    end subroutine determine_solver
        


end module det_solver 
