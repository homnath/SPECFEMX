module save_variables
implicit none
contains 
!_______________________________________________________________________________

subroutine save_potential_variables(i_step)
! USES
use global 
use postprocess
use free_surface
use math_constants
use set_precision
use nondimensionpar
#if (USE_MPI)
use ghost_library_mpi
#else
use serial_library
#endif 
implicit none
! IO variables: 
integer,intent(in) :: i_step

! Local variables: 
integer :: i_comp

! CODE: 
! plot gravity potential
if(savedata%gpot)then
  call write_scalar_to_file(nnode,DIM_GPOT*nodalphi,ext='gpot',istep=i_step) 
  ! On the free surface
  if(savedata%fsplot)then
    call write_scalar_to_file_freesurf(nnode_fs, DIM_GPOT*nodalphi(gnode_fs),&
    ext='gpot',istep=i_step) 
  endif
  if(savedata%fsplot_plane)then
    call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs),&
    ext='gpot',istep=i_step,plane=.true.) 
  endif
  ! Print confirmation
  if(myrank.eq.0.and.verbose_save_var)then
    write(*,'(a,i6)')'[OK] Saved gravity potential for step ', i_step
    write(*,*)
  endif 
endif
if(savedata%mpot)then
  call write_scalar_to_file(nnode,DIM_MPOT*nodalphi,ext='mpot',istep=i_step) 
  ! On the free surface
  if(savedata%fsplot)then
    call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs),&
    ext='mpot',istep=i_step) 
  endif
  if(savedata%fsplot_plane)then
    call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs),&
    ext='mpot',istep=i_step,plane=.true.) 
  endif
  ! Print confirmation
  if(myrank.eq.0.and.verbose_save_var)then
    write(*,'(a,i6)')'[OK] Saved magnetic potential for step ', i_step
    write(*,*)
  endif 
endif
if(savedata%epot)then
  call write_scalar_to_file(nnode,DIM_EPOT*nodalphi,ext='epot',istep=i_step)
  ! On the free surface
  if(savedata%fsplot)then
    call write_scalar_to_file_freesurf(nnode_fs,DIM_EPOT*nodalphi(gnode_fs), &
    ext='epot',istep=i_step) 
  endif
  if(savedata%fsplot_plane)then
    call write_scalar_to_file_freesurf(nnode_fs,DIM_EPOT*nodalphi(gnode_fs), &
    ext='epot',istep=i_step,plane=.true.) 
  endif
endif

! Gravitational
if(savedata%agrav)then
  ! Compute acceleration due to gravity
  call compute_gradient_of_scalar(nodalphi,nodalg)
  if(nproc.gt.1)then
    call assemble_ghosts_nodal_vector(nodalg,nodalg)
  endif
  ! compute average on the sharing nodes
  do i_comp=1,ndim
    nodalg(i_comp,:)=nodalg(i_comp,:)/real(node_valency,kreal)
  enddo
  ! Plot gravity accelration
  if(savedata%agrav)then
    call write_vector_to_file(nnode,DIM_G*nodalg,ext='grav',istep=i_step)
    ! On the free surface
    if(savedata%fsplot)then
      call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs),&
      ext='grav',istep=i_step)
    endif
    if(savedata%fsplot_plane)then
      call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs),&
      ext='grav',istep=i_step,plane=.true.)
    endif
  endif
  ! Print confirmation
  if(myrank.eq.0.and.verbose_save_var)then
    write(*,'(a,i6)')'[OK] Saved grav. acceleration for step ', i_step
    write(*,*)
  endif 
endif
! Magnetic
if(savedata%magb)then
  ! Compute magnetic field
  call compute_premagnetic_field(nodalphi,nodalB)
  if(nproc.gt.1)then
    call assemble_ghosts_nodal_vector(nodalB,nodalB)
  endif
  ! Compute average on the sharing nodes
  do i_comp=1,ndim
    nodalB(i_comp,:)=nodalB(i_comp,:)/real(node_valency,kreal)
  enddo
  ! Multiply by \mu_0
  nodalB=MAG_CONS*nodalB
  ! Plot magnetic field
  if(savedata%magb)then
    call write_vector_to_file(nnode,DIM_B*nodalB,ext='magb',istep=i_step)
    ! On the free surface
    if(savedata%fsplot)then
      call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs),&
      ext='magb',istep=i_step)
    endif
    if(savedata%fsplot_plane)then
      call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs),&
      ext='magb',istep=i_step,plane=.true.)
    endif
  endif
endif

end subroutine save_potential_variables
!===============================================================================

subroutine save_displacement_variables(i_step )
! USES
use global 
use postprocess
use free_surface
use math_constants
use math_library,only:rquick_sort,stress_invariant
use set_precision
use nondimensionpar
#if (USE_MPI)
use ghost_library_mpi
#else
use serial_library
#endif 
implicit none

! IO variables
integer,intent(in) :: i_step
! Local variables  
integer :: i_comp,i_node

real(kind=kreal):: dm 
real(kind=kreal) :: dsbar,lode_theta,sigm

! plot displacement
if(savedata%disp)then
  call write_vector_to_file(nnode,DIM_L*nodalu,ext='dis',istep=i_step) 
  ! On the free surface
  if(savedata%fsplot)then
    call write_vector_to_file_freesurf(nnode_fs,DIM_L*nodalu(:,gnode_fs),&
    ext='dis',istep=i_step)
  endif
  if(savedata%fsplot_plane)then
    call write_vector_to_file_freesurf(nnode_fs,DIM_L*nodalu(:,gnode_fs),&
    ext='dis',istep=i_step,plane=.true.)
  endif
  ! Print confirmation
  if(myrank.eq.0.and.verbose_save_var)then
    write(*,'(a,i6)')'[OK] Saved displacement for step ', i_step
    write(*,*)
  endif
endif

! stress and/or principal stress
! First compute the nodal stress
if(savedata%stress .or. savedata%psigma)then
  call compute_nodal_tensor(stress_elmt,stress_nodal)  
  if(nproc.gt.1)then
    call assemble_ghosts_nodal_vectorn(NST,stress_nodal,stress_nodal)
  endif
  ! compute average on the sharing nodes
  do i_comp=1,NST
    stress_nodal(i_comp,:)=stress_nodal(i_comp,:)/real(node_valency,kreal)
  enddo
endif
! Plot nodal stress.
if(savedata%stress)then
  call write_vector_to_file(nnode,DIM_MOD*stress_nodal,&
  ext='sig',istep=i_step)
  ! On the free surface
  if(savedata%fsplot)then
    call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
    ext='sig',istep=i_step)
  endif
  if(savedata%fsplot_plane)then
    call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
    ext='sig',istep=i_step,plane=.true.)
  endif
  ! Print confirmation
  if(myrank.eq.0.and.verbose_save_var)then
    write(*,'(a,i6)')'Saved stress for step ', i_step
    write(*,*)
  endif
endif
! Plot principal nodal stress.
if(savedata%psigma)then
  ! Compute pricipal stress from nodal stress.
  pstress_nodal=ZERO
  do i_node=1,nnode
    call stress_invariant(stress_nodal(:,i_node),sigm,dsbar,lode_theta)
    pstress_nodal(:,i_node)=sigm+TWO_THIRD*dsbar* &
          sin((/ lode_theta-TWO_THIRD*PI,lode_theta,lode_theta+TWO_THIRD*PI /))
    ! put in ascending order (compresson (-) as a major principal stress!)
    pstress_nodal(:,i_node)=rquick_sort(pstress_nodal(:,i_node),3)
  enddo

  call write_vector_to_file(nnode,DIM_MOD*pstress_nodal,&
  ext='psig',istep=i_step)
  ! On the free surface
  if(savedata%fsplot)then
    call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*pstress_nodal(:,gnode_fs),&
    ext='psig',istep=i_step)
  endif
  if(savedata%fsplot_plane)then
    call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*pstress_nodal(:,gnode_fs),&
    ext='psig',istep=i_step,plane=.true.)
  endif
  ! Print confirmation
  if(myrank.eq.0.and.verbose_save_var)then
    write(*,'(a,i6)')'Saved pricipal stress for step ', i_step
    write(*,*)
  endif
endif
!! plot stress
!if(savedata%stress)then
!  call compute_nodal_tensor(stress_elmt,stress_nodal)  
!  if(nproc.gt.1)then
!    call assemble_ghosts_nodal_vectorn(NST,stress_nodal,stress_nodal)
!  endif
!  ! compute average on the sharing nodes
!  do i_comp=1,NST
!    stress_nodal(i_comp,:)=stress_nodal(i_comp,:)/real(node_valency,kreal)
!  enddo
!  call write_vector_to_file(nnode,DIM_MOD*stress_nodal,&
!  ext='sig',istep=i_step)
!  ! On the free surface
!  if(savedata%fsplot)then
!    call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
!    ext='sig',istep=i_step)
!  endif
!  if(savedata%fsplot_plane)then
!    call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
!    ext='sig',istep=i_step,plane=.true.)
!  endif
!  ! Print confirmation
!  if(myrank.eq.0.and.verbose_save_var)then
!    write(*,'(a,i6)')'[OK] Saved stress for step ', i_step
!    write(*,*)
!  endif
!endif

! plot strain
if(savedata%strain)then
  call compute_nodal_tensor(strain_elmt,strain_nodal)  
  if(nproc.gt.1)then
    call assemble_ghosts_nodal_vectorn(NST,strain_nodal,strain_nodal)
  endif
  ! compute average on the sharing nodes
  do i_comp=1,NST
    strain_nodal(i_comp,:)=strain_nodal(i_comp,:)/real(node_valency,kreal)
  enddo
  call write_vector_to_file(nnode,strain_nodal,ext='eps',istep=i_step)
  ! On the free surface
  if(savedata%fsplot)then
    call write_vector_to_file_freesurf(nnode_fs,strain_nodal(:,gnode_fs),&
    ext='eps',istep=i_step)
  endif
  if(savedata%fsplot_plane)then
    call write_vector_to_file_freesurf(nnode_fs,strain_nodal(:,gnode_fs),&
    ext='eps',istep=i_step,plane=.true.)
  endif
  ! Print confirmation
  if(myrank.eq.0.and.verbose_save_var)then
    write(*,'(a,i6)')'[OK] Saved strain for step ', i_step
    write(*,*)
  endif
  if(trim(devel_example).eq.'axial_rod')then
    write(77,*)0.0,strain_nodal(1,2099)                               
    flush(77) 
  endif
endif

end subroutine save_displacement_variables
!===============================================================================

end module save_variables
!===============================================================================
