!===============================specfem3d.F90===================================
! REVISION:
!   HNG, Aug 25,2011; HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
!-------------------------------------------------------------------------------
module write_data
    contains
    !-------------------------------------------------------------------------------
    subroutine write_data_step(i_tstep,node_valency)
    ! Import necessary modules
    use dimensionless
    use global
    use string_library,only : parse_file
    use math_constants
    use conversion_constants
    use math_library
    use element,only:hex8_gnode
    use global_dof
    use free_surface
    use write_ensight
    use postprocess
    implicit none
    
    integer,intent(in) :: i_tstep
    integer,dimension(:),intent(in) :: node_valency ! Nodal valency
    character(len=250) :: myfname=' => write_data_step.f90'
    character(len=500) :: errsrc
    
    ! i,j are dummy vars for interation
    integer :: i,j
    ! istat: status indicator for allocation (can be used in other contexts)
    integer :: istat
    
    ! do-loop indices
    integer :: i_elmt,i_nliter,i_node,i_comp
    integer :: ielmt,imat,idof,iedof!element ID for gdof, node, etc.
    
    real(kind=kreal),dimension(nst),parameter :: unit_voigt=(/one,one,one,ZERO,    &
    ZERO,ZERO /)!Voigt representation for vector
    real(kind=kreal) :: jacw !determinant of Jacobian*gll_weight
    real(kind=kreal) :: dt,dt_vp ! time step
    
    real(kind=kreal) :: sfac ! slip factor
    
    ! KSP convergence reason
    integer :: ksp_convreason
    ! ksp_iter: conjugate gradient iteration,
    ! ksp_tot: total cg interation, nl_iter: nonlinear iteration, nl_tot: total nl
    ! iteration
    integer :: ksp_iter,ksp_tot,nl_iter,nl_tot
    logical :: nl_isconv ! logical variable to check convergence of
    ! nonlinear (NL) iterations
    
    
    !strain_elmt: strain for all elements
    !stress_elmt: stress for each element
    !stress_nodal: nodal stress for all elements in processor
    real(kind=kreal),allocatable :: strain_elmt(:,:,:),strain_nodal(:,:),      &
    stress_elmt(:,:,:),stress_nodal(:,:),evpt(:,:,:)
    !bcnodalv: prescribed BC nodal variables
    !nodalu: nodal displacement for all nodes (not just BC)
    real(kind=kreal),allocatable :: bcnodalv(:,:),nodalu(:,:)
    real(kind=kreal),allocatable :: nodalphi(:),nodalg(:,:)
    ! magnetization
    real(kind=kreal),allocatable :: nodalB(:,:)
    
    ! Frequency
    logical :: isscale_ang_freq=.true.
    real(kind=kreal) :: freq,ang_freq,scale_ang_freq2
    
    logical :: isgravity,ispseudoeq ! gravity load and pseudostatic load
    
    ! Time at current time step
    real(kind=kreal) :: t
    integer :: i_step,istep
    integer :: istep0 !first step
    real(kind=kreal) :: step !current step (t or f)
    
    !Factor for time unit conversion
    real(kind=kreal) :: tunitfac
    real(kind=kreal) :: cpu_tstart,cpu_tend,telap,max_telap,mean_telap
    character(len=20) :: format_str
    
    integer :: i_maxwell
    integer,allocatable :: eid_elas(:),eid_viscoelas(:)
    real(kind=kreal),allocatable :: relaxtime(:,:),muratio(:),tratio(:)
    ! q0: initial state variable for viscoelastic rheology
    real(kind=kreal),allocatable :: visco_q0(:,:,:,:),q0(:,:),elas_e0(:,:,:)
    real(kind=kreal) :: e0(nst) !e0: initial strain
    
    real(kind=kreal) :: vesigma(nst)
    real(kind=kreal) :: trace_vsigma0,trace_strain
    real(kind=kreal) :: esigma0_dev(nst),esigma_dev(nst)
    real(kind=kreal) :: maxresload,maxbodyload
    
    integer :: geq,inum,nequ
    logical,allocatable :: iseq(:)
    integer,allocatable :: gdofu(:) 
    integer :: tot_neq,max_neq,min_neq
    ! number of active ghost partitions for a node
    integer,allocatable :: ngpart_node(:)
    character(len=250) :: errtag ! error message
    integer :: errcode
    errtag=""; errcode=-1
    
    if(ISDISP_DOF)then
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
      endif
      !if(devel_mgll)then
      !  call compute_save_density_perturbation(nodalu,errcode,errtag)
      !endif
    endif
    if(ISPOT_DOF)then
      ! plot gravity potential
      if(savedata%gpot)then
        call write_scalar_to_file(nnode,DIM_GPOT*nodalphi,ext='gpot',istep=i_step) 
        ! On the free surface
        if(savedata%fsplot)then
          call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs), &
          ext='gpot',istep=i_step) 
        endif
        if(savedata%fsplot_plane)then
          call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs), &
          ext='gpot',istep=i_step,plane=.true.) 
        endif
      endif
      if(savedata%mpot)then
        call write_scalar_to_file(nnode,DIM_MPOT*nodalphi,ext='mpot',istep=i_step) 
        ! On the free surface
        if(savedata%fsplot)then
          call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs), &
          ext='mpot',istep=i_step) 
        endif
        if(savedata%fsplot_plane)then
          call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs), &
          ext='mpot',istep=i_step,plane=.true.) 
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
            call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs), &
            ext='grav',istep=i_step)
          endif
          if(savedata%fsplot_plane)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs), &
            ext='grav',istep=i_step,plane=.true.)
          endif
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
            call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs), &
            ext='magb',istep=i_step)
          endif
          if(savedata%fsplot_plane)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs), &
            ext='magb',istep=i_step,plane=.true.)
          endif
        endif
      endif
    endif
    
    return
    end subroutine write_data_step
    !===============================================================================
    