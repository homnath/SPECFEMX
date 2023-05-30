module bcs_and_dof
implicit none

contains
!_______________________________________________________________________________
subroutine sort_gdofs_and_bc()

    ! This was previously a large part of the specfem3d script. Overall this
    ! section does the following: 
    ! 1) Allocate necessary boundary condition/gdof arrays 
    ! 2) Work out which DOFs should be on based on what is being solved for
    !    (in activate_dof)
    ! 3) Ensures that the global IDs for these DOFs are the same across
    !    processors if using parallel run 
    ! 4) Looks at BCs and determines which DOFs no longer need to be solved
    !    for because they are determined by the BCs 
    ! 5) Create global indexing for Petsc if required

use global
use math_constants
use set_precision
use output_to_user
#if(USE_MPI)
use mpi_library
use math_library_mpi
use ghost_library_mpi
#else
use serial_library
use math_library_serial
#endif
use global_dof
use bc
use dof

implicit none 

! Local 
integer :: istat, errcode, i_elmt
character(len=250) :: errtag

write(logunit,*)' -------------------------------------------------'
write(logunit,*)' Sorting GDOFs and BCs   (sort_gdofs_and_bc)'
write(logunit,*)

! Initialise boundary conditions
bcnodalv=ZERO

! Initialise infinite element faces
allocate(infinite_iface(6,nelmt),infinite_face_idir(6,nelmt))
infinite_iface=.false.
infinite_face_idir=-9999

! Activate the degrees of freedom 
call activate_dof(errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)

! Ensure that gdof IDs are same in the finite/infinite interface
! nodes if they lie across different processors.
! This can also be done if we explicitly define the gdof ON/OFF state on those
! inteface nodes across all the processors.
call assemble_ghosts_gdof(nndof,gdof,gdof)

where(gdof>0)gdof=1
call sync_process

! Apply Dirichlet boundary conditions
call apply_bc(errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)

! Finalise the GDOF after BCs have been applied 
call finalize_gdof(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)
log_msg = 'complete!' ; call write_ifproc0(logunit)

!call modify_ghost_gdof(num, egdof, egdofu, coord, deriv, jac, bmat, &
!        eld, eload, bload, vload, nodalu, nodalphi, nodalg, nodalB, nodalphistore, nodalustore)

!-------------------------------------
! modify ghost GDOFs
if(nproc.gt.1)then
  call prepare_ghost_gdof()
endif

! store elemental global degrees of freedoms from nodal gdof
! this removes the repeated use of reshape later but it has larger size than
! gdof!!!
gdof_elmt=0
do i_elmt=1,nelmt
  gdof_elmt(:,i_elmt)=reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
enddo

! global indexing
! this process is necessary for petsc implementation and is done in a single
! processor
call sync_process
if(ismpi .and. nproc.gt.1 .and. myrank.eq.0)then 
  call gindex()
  write(logunit,*)'CALLED GINDEX'
endif 
call sync_process


! Output details to user 
tot_neq=sumscal(neq); max_neq=maxscal(neq); min_neq=minscal(neq)
if(myrank==0)then
  write(logunit,'(a,i0,a,i0,a,i0)')'degrees of freedoms => total:',tot_neq,&
                                  ' max:',max_neq,' min:',min_neq
  flush(logunit)
endif

write(logunit,*)' ✓  Finished sorting GDOFs with BCs'
write(logunit,*)


end subroutine sort_gdofs_and_bc
! ______________________________________________________________________

end module  bcs_and_dof
