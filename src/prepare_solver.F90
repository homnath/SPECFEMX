module prepare_solver 
contains 
!_______________________________________________________________________________

! WE - prepare solver; originally in nonlinear loop
subroutine prep_inbuilt_solver(nzero_dprecon, nelmt_elas, eid_elas,    &
                               nelmt_viscoelas, eid_viscoelas)
    ! USES: 
use global
use local
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
    integer,allocatable :: eid_elas(:), eid_viscoelas(:)
    integer :: nzero_dprecon, nelmt_elas, nelmt_viscoelas

    ! Local 
    integer :: i_elmt, ielmt, j_dof, i_eq, i_dof
  
    ! CODE: 
      ! Compute diagonal precoditioner
      dprecon=ZERO
      do i_elmt=1,nelmt
        ielmt=i_elmt ! all elements
        egdof=gdof_elmt(:,ielmt)
        do j_dof=1,nedof
          dprecon(egdof(j_dof))=dprecon(egdof(j_dof))+storekmat(j_dof,j_dof,ielmt)
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
        do i_eq=1,neq
          ndscale(i_eq)=one/sqrt(abs(dprecon(i_eq)))
        enddo
        ! nondimensionalize stiffness matrix
        ! elastic region
        do i_elmt=1,nelmt_elas
          ielmt=eid_elas(i_elmt)
          egdof=gdof_elmt(:,ielmt)
          do i_dof=1,nedof
            do j_dof=1,nedof
              storekmat(i_dof,j_dof,ielmt)=ndscale(egdof(i_dof))*storekmat(i_dof,j_dof,ielmt)*       &
              ndscale(egdof(j_dof))
            enddo
          enddo
        enddo
        ! viscoelastic region
        do i_elmt=1,nelmt_viscoelas
          ielmt=eid_viscoelas(i_elmt)
          egdof=gdof_elmt(:,ielmt)
          do i_dof=1,nedof
            do j_dof=1,nedof
              storekmat(i_dof,j_dof,ielmt)=ndscale(egdof(i_dof))*storekmat(i_dof,j_dof,ielmt)*       &
              ndscale(egdof(j_dof))
            enddo
          enddo
        enddo
      else
        dprecon(1:)=one/dprecon(1:)
      endif
  
  
  end subroutine prep_inbuilt_solver
!-------------------------------------------------------------------------------  
  

end module prepare_solver
!===============================================================================
