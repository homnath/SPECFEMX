module prepare_solver 
contains 




! WE - prepare solver; originally in nonlinear loop
subroutine prep_inbuilt_solver(dprecon, egdof, storekmat, ndscale,     &
                               nzero_dprecon, nelmt_elas, eid_elas,    &
                               nelmt_viscoelas, eid_viscoelas)
    ! USES: 
use global
use math_constants
use set_precision
#if (USE_MPI)
use ghost_library_mpi
use mpi_library
#else
use serial_library
#endif

    implicit none 
    ! IO 
    real(kind=kreal),allocatable :: dprecon(:), storekmat(:,:,:), ndscale(:)
    integer,allocatable :: egdof(:), eid_elas(:), eid_viscoelas(:)
    integer :: nzero_dprecon, nelmt_elas, nelmt_viscoelas

    ! Local 
    integer :: i_elmt, ielmt, j, i
  
    ! CODE: 
    if(solver_type.eq.builtin_solver)then
      ! Compute diagonal precoditioner
      dprecon=ZERO
      do i_elmt=1,nelmt
        ielmt=i_elmt ! all elements
        egdof=gdof_elmt(:,ielmt)
        do j=1,nedof
          dprecon(egdof(j))=dprecon(egdof(j))+storekmat(j,j,ielmt)
        enddo
      enddo ! i_elmt
      dprecon(0)=ZERO
  
      ! assemble diagonal preconditioner
      if(nproc.gt.1)then
        call assemble_ghosts(nndof,neq,dprecon,dprecon)
      endif
      call sync_process()
      dprecon(0)=ZERO
      nzero_dprecon=count(dprecon(1:).eq.ZERO)
      if(nzero_dprecon.gt.0)then
        write(logunit,*)'WARNING: nzero in dprecon:', &
        nzero_dprecon,minval(abs(dprecon(1:))),maxval(abs(dprecon(1:)))
        flush(logunit)
      endif
  
      if(solver_diagscale)then
        ! regularize linear equations,
        ndscale=ONE
        do i=1,neq
          ndscale(i)=one/sqrt(abs(dprecon(i)))
        enddo
        ! nondimensionalize stiffness matrix
        ! elastic region
        do i_elmt=1,nelmt_elas
          ielmt=eid_elas(i_elmt)
          egdof=gdof_elmt(:,ielmt)
          do i=1,nedof
            do j=1,nedof
              storekmat(i,j,ielmt)=ndscale(egdof(i))*storekmat(i,j,ielmt)*       &
              ndscale(egdof(j))
            enddo
          enddo
        enddo
        ! viscoelastic region
        do i_elmt=1,nelmt_viscoelas
          ielmt=eid_viscoelas(i_elmt)
          egdof=gdof_elmt(:,ielmt)
          do i=1,nedof
            do j=1,nedof
              storekmat(i,j,ielmt)=ndscale(egdof(i))*storekmat(i,j,ielmt)*       &
              ndscale(egdof(j))
            enddo
          enddo
        enddo
      else
        dprecon(1:)=one/dprecon(1:)
      endif
    endif !(solver_type.eq.builtin_solver)
  
  
  end subroutine prep_inbuilt_solver
  
  
  

end module prepare_solver