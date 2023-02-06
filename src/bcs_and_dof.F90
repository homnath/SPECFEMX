module bcs_and_dof
implicit none

contains
! ______________________________________________________________________
subroutine sort_gdofs_and_bc(bcnodalv, num, egdof, egdofu, coord, deriv,&
    eld, eload, bload, vload, rhoload, resload, jac, bmat, nodalu, & 
    nodalg, nodalphi, nodalB, tot_neq, max_neq, min_neq, currentu)

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
use ghost
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

real(kind=kreal), allocatable :: bcnodalv(:,:), nodalu(:,:), currentu(:,:), nodalphi(:),nodalg(:,:), nodalB(:,:)
integer,allocatable::num(:)
integer,allocatable :: egdof(:),egdofu(:)
integer :: tot_neq,max_neq,min_neq
real(kind=kreal),allocatable :: bmat(:,:),coord(:,:),deriv(:,:),    &
jac(:,:)
real(kind=kreal),allocatable :: eld(:),eload(:),bload(:),   &
vload(:),rhoload(:),resload(:)


! Local 
integer :: istat, errcode, i_elmt
character(len=250) :: errtag

write(logunit,*)' -------------------------------------------------'
write(logunit,*)' Sorting GDOFs and BCs   (sort_gdofs_and_bc)'
write(logunit,*)

! Initialise boundary conditions
allocate(bcnodalv(nndof,nnode))
bcnodalv=ZERO

! Initialise global degrees of freedom database
allocate(gdof(nndof,nnode),stat=istat)
call check_memory_alloc(istat, logunit)

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
call apply_bc(bcnodalv,errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)


! Finalise the GDOF after BCs have been applied 
call finalize_gdof(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)
log_msg = 'complete!' ; call write_ifproc0(logunit)

call modify_ghost_gdof(num, egdof, egdofu, coord, deriv, jac, bmat, &
        eld, eload, bload, vload, nodalu, nodalphi, nodalg, nodalB, currentu)


! store elemental global degrees of freedoms from nodal gdof
! this removes the repeated use of reshape later but it has larger size than
! gdof!!!
allocate(gdof_elmt(nedof,nelmt))
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








end module 