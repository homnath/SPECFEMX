module save_variables
contains 


subroutine save_pot_variables(nodalphi,nodalg, nodalB, node_valency, i_step)
! USES
use global 
use postprocess
use free_surface
use math_constants
use set_precision
use dimensionless
#if (USE_MPI)
use ghost_library_mpi
#else
use serial_library
#endif 

    ! IO variables: 
    real(kind=kreal),allocatable :: nodalphi(:),nodalg(:,:), nodalB(:,:)
    integer, allocatable :: node_valency(:)
    integer :: i_step

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
      if(myrank.eq.0)then
        write(*,'(a,i6)')'  ✓ Saved gravity potential for step ', i_step
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
      if(myrank.eq.0)then
        write(*,'(a,i6)')'  ✓ Saved magnetic potential for step ', i_step
        write(*,*)
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
      if(myrank.eq.0)then
        write(*,'(a,i6)')'  ✓ Saved grav. acceleration for step ', i_step
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

end subroutine save_pot_variables





subroutine save_displacement_variables(strain_elmt, strain_nodal, &
                                       stress_elmt, stress_nodal, & 
                                       nodalu, node_valency, i_step )

    ! USES
use global 
use postprocess
use free_surface
use math_constants
use set_precision
use dimensionless
#if (USE_MPI)
use ghost_library_mpi
#else
use serial_library
#endif 

  ! IO variables
  ! Local variables  
  real(kind=kreal),allocatable :: nodalu(:,:)
  integer, allocatable :: node_valency(:)
  integer :: i_step, i_comp
  real(kind=kreal),allocatable :: strain_elmt(:,:,:),strain_nodal(:,:),      &
                                  stress_elmt(:,:,:),stress_nodal(:,:)

  real(kind=kreal):: dm 


  ! CODE: 
  
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
    if(myrank.eq.0)then
      write(*,'(a,i6)')'  ✓ Saved displacement for step ', i_step
      write(*,*)
    endif
  endif


  ! plot stress
  if(savedata%stress)then
    call compute_nodal_tensor(stress_elmt,stress_nodal)  
    if(nproc.gt.1)then
      call assemble_ghosts_nodal_vectorn(NST,stress_nodal,stress_nodal)
    endif
    ! compute average on the sharing nodes
    do i_comp=1,NST
      stress_nodal(i_comp,:)=stress_nodal(i_comp,:)/real(node_valency,kreal)
    enddo
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
    if(myrank.eq.0)then
      write(*,'(a,i6)')'  ✓ Saved stress for step ', i_step
      write(*,*)
    endif
  endif


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
    call write_vector_to_file(nnode,strain_nodal,&
    ext='eps',istep=i_step)
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
    if(myrank.eq.0)then
      write(*,'(a,i6)')'  ✓ Saved strain for step ', i_step
      write(*,*)
    endif
    if(trim(devel_example).eq.'axial_rod')then
      write(77,*)0.0,strain_nodal(1,2099)                               
      flush(77) 
    endif
  endif

end subroutine save_displacement_variables





end module save_variables