!===============================specfem3d.F90===================================
! REVISION:
!   HNG, Aug 25,2011; HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
!-------------------------------------------------------------------------------
subroutine specfem3d()
! Import necessary modules
use dimensionless
use global
use output_to_user
use string_library,only : parse_file
use math_constants
use conversion_constants
use gll_library
use shape_library
use infinite_element
use math_library
use element,only:hex8_gnode
use dof
use fault
use weakform
use gravity
use preprocess
use matrix_vector
use elastic
use traction
use mtraction
use viscoelastic
use global_dof
use plastic_library
use map_location
use earthquake
use source_function
use cmtsolution,only:source_tshift,source_hdur
#if (USE_MPI)
use mpi_library
use ghost_library_mpi
use math_library_mpi
use sparse
use parsolver
use count_elements  !WE
use relaxation_time !WE
use ghost           !WE
use other_forces    !WE 
#if (USE_COMPLEX)
use parsolver_petsc_complex
#else
use parsolver_petsc
#endif
#else
use serial_library
use math_library_serial
use sparse_serial
use solver
use solver_petsc
#endif
use prestress
use bc
use prepare_solver !WE
use free_surface
use benchmark
use write_ensight
use postprocess
use cleanup
use save_variables 
use nonlinearloop
!use bilinear_form
use sea_level
use ice
use time_loop
use bcs_and_dof
use initialise_arrays

!use relaxation_time
implicit none

character(len=250) :: myfname=' => specfem3d.f90'
character(len=500) :: errsrc


! istat: status indicator for allocation (can be used in other contexts)
integer :: istat, numberproc, myproc

! do-loop indices
integer :: i_dof,i_elmt,i_eq,i_gll,i_mat,i_nliter,i_node,i_comp,j_dof,j_node
integer :: ielmt,imat,idof,iedof!element ID for gdof, node, etc.

real(kind=kreal),dimension(nst),parameter :: unit_voigt=(/one,one,one,ZERO,    &
ZERO,ZERO /)!Voigt representation for vector
real(kind=kreal) :: jacw !determinant of Jacobian*gll_weight
real(kind=kreal) :: dt,dt_vp ! time step

real(kind=kreal) :: sfac ! slip factor

! ksp_iter: conjugate gradient iteration,
! ksp_tot: total cg interation, nl_iter: nonlinear iteration, nl_tot: total nl
! iteration
integer :: ksp_iter,ksp_tot,nl_iter,nl_tot
logical :: nl_isconv ! logical variable to check convergence of
! nonlinear (NL) iterations

real(kind=kreal) :: f

integer :: todalOFNodes
! cmat: elastic matric (Cijkl) in Voigt notation
! estrain: elastic strain
! esigma: elastic stress
! sigma: stress (used in different context than esigma?)
! vsigma: viscose stress

! Dynamic arrays
! num: g_num for particular element.
! node_valency: number of elements that share each node.
integer,allocatable::num(:),node_valency(:)

integer :: nzero_dprecon, spareint
real(kind=kreal),allocatable :: dprecon(:),ndscale(:)

real(kind=kreal),allocatable :: bmat(:,:),coord(:,:),deriv(:,:),    &
jac(:,:)
!kmat: stiffness matrix for each element
!storekmat: stiffness matrix for all elements
real(kind=kreal),allocatable :: kmat(:,:),storekmat(:,:,:)
!storemmat: mass matrix for all elements
real(kind=kreal),allocatable :: storemmat(:,:)

!uerr: used to check convergence
!umax: max of displacement magnitude
!uxmax: max of displacement components
!du: incremental solution
!u: solution (summed over du)
real(kind=kreal),allocatable :: du(:),u(:),olddu(:)
real(kind=kreal) :: uerr

integer :: i

! Source frequency function
real(kind=kreal),allocatable :: eld(:),eload(:),bload(:),   &
vload(:),rhoload(:),resload(:)
!eld: elastic displacement on all nodes of the element
!eload: elastic load on element
!extload: external load
!ubcload: load contributed by displacement BC
!load: like resload. Not currently used
real(kind=kreal),allocatable :: slipload(:),extload(:),bodyload(:),selfload(:), &
viscoload(:),ubcload(:),load(:)

!strain_elmt: strain for all elements
!stress_elmt: stress for each element
!stress_nodal: nodal stress for all elements in processor
real(kind=kreal),allocatable :: strain_elmt(:,:,:),strain_nodal(:,:),      &
stress_elmt(:,:,:),stress_nodal(:,:),evpt(:,:,:)
!bcnodalv: prescribed BC nodal variables
!nodalu: nodal displacement for all nodes (not just BC)
real(kind=kreal),allocatable :: bcnodalv(:,:),nodalu(:,:), nodalustore(:,:) 
real(kind=kreal),allocatable :: nodalphi(:),nodalphistore(:),nodalg(:,:)

! Sea level edit - WE 
real(kind=kreal),allocatable :: nodalsl(:)  ! Nodal theta values
real(kind=kreal),allocatable :: nodalslrate(:)  ! Nodal theta values
real(kind=kreal),allocatable :: nodalice(:) ! Nodal I values
real(kind=kreal),allocatable :: nodalicerate(:) ! Nodal rate of I values
real(kind=kreal),allocatable :: iceload(:)  ! load term due to ice.
real(kind=kreal):: maxnodalsl,minnodalsl
! magnetization
real(kind=kreal),allocatable :: nodalB(:,:)
!,psigma(:,:),psigma0(:,:),taumax(:),nsigma(:)
!bodyload: load computed on all nodes of the element
!selfload: self load computed  on all nodes of the element 
!viscoload: load contributed (on all nodes) by visco elements
!bmat: strain displacement matrix (per node). bmat*displacement vector = strain
!bload: load computed on each dof on a particular element
!vload: like bload, for viscose load
!resload: residual load
!coord: coordinates of geometrical nodes
!deriv: derivative of interpolation functions
!dprecon: diagonal preconditioner. Not used for petsc solver
!ndscale: nondimensionlize scale
!ubcload: load contributed by displacement BC
!jac: Jacobian
integer,allocatable :: egdof(:),egdofu(:)
! placeholder array. holds values of gdof_elmt for a given element.

! Frequency
logical :: isscale_ang_freq=.true.
real(kind=kreal) :: freq,ang_freq,scale_ang_freq2

! Viscoelastic parameters
integer :: mdomain
integer :: ielmt_elas,ielmt_viscoelas,iviscoelas,nelmt_elas,            &
nelmt_viscoelas
integer :: tot_nelmt_elas,max_nelmt_elas,min_nelmt_elas
integer :: tot_nelmt_viscoelas,max_nelmt_viscoelas,min_nelmt_viscoelas
integer :: nmatblk_elas
real(kind=kreal) :: min_relaxtime,max_relaxtime

! Time at current time step
real(kind=kreal) :: t
integer :: i_step,istep
integer :: istep0 !first step
real(kind=kreal) :: step !current step (t or f)

!Factor for time unit conversion
real(kind=kreal) :: tunitfac

integer :: i_maxwell
integer,allocatable :: eid_elas(:),eid_viscoelas(:)
real(kind=kreal),allocatable :: relaxtime(:,:),muratio(:),tratio(:)
! q0: initial state variable for viscoelastic rheology
real(kind=kreal),allocatable :: visco_q0(:,:,:,:),q0(:,:),elas_e0(:,:,:)


real(kind=kreal) :: trace_vsigma0
real(kind=kreal) :: esigma0_dev(nst),esigma_dev(nst)

integer :: geq,inum,nequ
logical,allocatable :: iseq(:)
integer,allocatable :: gdofu(:) 
integer :: tot_neq,max_neq,min_neq
! number of active ghost partitions for a node
integer,allocatable :: ngpart_node(:)
character(len=250) :: errtag ! error message
integer :: errcode
errtag=""; errcode=-1

! ______________________________________________________________________
! ______________________________________________________________________
! Calculate relaxation time for viscoelastic/plastic models
call calc_relaxation_time(relaxtime, tunitfac, muratio, & 
                          tratio, min_relaxtime, max_relaxtime)


!Count number of viscoelastic and elastic elements 
call count_elmts(errcode, nelmt_elas, nelmt_viscoelas,     & 
                 mdomain, tot_nelmt_elas, max_nelmt_elas,  &
                 min_nelmt_elas, tot_nelmt_viscoelas,      &
                 max_nelmt_viscoelas, min_nelmt_viscoelas)


! Split up IDs of elastic and viscoelastic elements into new arrays
! New arrays: eid_elas and eid_viscoelas
call split_elas_visco_eids(nelmt_elas, nelmt_viscoelas, eid_elas,&
                           eid_viscoelas)


! Create new elastic/viscoelastic arrays
allocate(elas_e0(nst,ngll,nelmt_viscoelas), &
         visco_q0(nst,nmaxwell,ngll,nelmt_viscoelas),q0(nst,nmaxwell))


! prepare ghost partitions for the communication
if(nproc.gt.1)then
  call prepare_ghost()
endif
deallocate(g_num0) ! Old connectivity no longer necessary
call sync_process()


! prepare fault - Only for Fault-based simulations 
if( iseqsource .and. (eqsource_type.eq.3 .or. eqsource_type.eq.4) )then
  call prepare_fault(errcode,errtag)
  call sync_process
  call control_error(errcode,errtag,stdout,myrank)
  if(myrank.eq.0)then 
    write(*,*)'Created split fault'
  endif 
endif


! Use BCs to determine global DOF index etc
call sort_gdofs_and_bc(bcnodalv, num, egdof, egdofu, coord, deriv, &
                       eld, eload, bload, vload, rhoload, resload, &
                       jac, bmat, nodalu, nodalg, nodalphi, nodalB,& 
                       tot_neq, max_neq, min_neq, nodalphistore, nodalustore)


! Calculate any prestress 
call calculate_prestress(strain_elmt, strain_nodal,       &
                         stress_elmt, stress_nodal,       & 
                         extload, du, dprecon, storekmat, &
                         errcode, errtag, ksp_iter, istat)
   

! compute node valency and assemble all node_valency across processors
call calculate_valency(node_valency, num)
call assemble_ghosts_nodal_iscalar(node_valency,node_valency)


! Update log with KSP details
call log_KSP_summary()


! Initialise RHS vectors + stiffness matrix/mass matrix 
call initialise_RHS_vectors(load, bodyload, selfload, viscoload, &
                            resload, du, u, kmat, storekmat,     &
                            storemmat, rhoload, ubcload, nodalu, & 
                            visco_q0, elas_e0, extload, iceload, & 
                            nodalustore, nodalphistore)


! Timestepping only needed for plastic/viscoelastic situations
if(isplastic)then
  allocate(olddu(0:neq),evpt(nst,ngll,nelmt))
endif
allocate(ngpart_node(nnode))


! compute stable time step for implicit integration
dt=dstep
nl_tot=0


! Prepare PETSC solver
if(solver_type.eq.petsc_solver)then
  ! Prepare sparsity of the stiffness matrix (working out size etc)
  call prepare_sparse() 
  ! petsc solver
  call petsc_initialize() 
  if(myrank.eq.0)then 
    write(*,*)'- PETSC initialisation        : Successful'
  endif 
  ! Create sparse vector, matrix, and preallocate                                                         
  ! TODO: following call is not necessary for RECYCLE         
  call petsc_create_vector()                                                     
  call petsc_matrix_preallocate_size()                                           
  call petsc_create_matrix()                                                     
  call petsc_create_solver()                           
  if(myrank.eq.0)then 
    write(*,*)'- PETSC matrix/vector creation: Successful'
    write(*,*)
  endif
endif 




! WE - only for a slipping fault with lobe split?? 
! WARNING: TODO
! slip gdof for split PC
if(ISDISP_DOF)then
  allocate(iseq(neq))
  iseq=.false.
  do i_node=1,nnode
    do i_dof=1,nndofu
      geq=gdof(idofu(i_dof),i_node)
      if(geq.gt.0)iseq(geq)=.true.
    enddo
  enddo
  nequ=count(iseq)
  allocate(gdofu(nequ))
  gdofu=-9999
  inum=0
  do i_eq=1,neq
    if(iseq(i_eq))then
      inum=inum+1
      gdofu(inum)=i_eq
    endif
  enddo
  if(nequ.ne.inum)then
    write(*,*) 'ERROR: nequ & inum mismatch!'
    call close_process
  endif
  deallocate(iseq)
endif



! prepare background gravity data
call prepare_gravity()


! Allocate for built-in solver preconditioner stuff? 
if(solver_type.eq.builtin_solver .or. solver_diagscale)then
  allocate(dprecon(0:neq))
endif
if(solver_diagscale)then
  allocate(ndscale(0:neq))
endif


! WE commented out
!if(trim(devel_example).eq.'axial_rod')then
!  open(77,file=trim(file_head)//"_strain.dat",action="write",status="replace")
!endif


! Set initial timestep and output to user: 
!   1 for time domain.
!   0 for frequency domain.
istep0=1
if(steptype.eq.FREQSTEP)then
  istep0=0
  if(myrank.eq.0)then
    print*,'Step type: FREQUENCY'
    print*,'f0, f1, df (Hz):',step0,step1,dstep
  endif
endif


! Compute elastic mass matrix once and for all 
call compute_mass_elastic(storemmat,errcode,errtag)


if(isplastic)then
  ! Compute minimum pseudo-time step for viscoplasticity
  dt_vp=dt_viscoplas(nmatblk,nu_blk,phi_blk,ym_blk)
  ! find global dt_vp
  dt_vp=minscal(dt_vp)
endif


! Compute body loads
if(isbodyload)then
  call compute_bodyload(selfload,selfweight=isselfweight)
endif 

call sync_process()



! Initialise ice: 
if(is_ICE)then
  ! Prepare the ice stuff and set the user-inputted initial condition
  call prepare_ice(nodalice, nodalicerate)
  call set_original_ice_level(nodalice)
endif 


! Initialise Sea Level 
if(is_SL)then 
  call prepare_sea_level(nodalsl, nodalslrate)
  call set_original_sea_level(nodalsl)

  ! Calc ocean func using initial SL and output if desired
  call update_ocean_function(nodalice, nodalsl, errcode, errtag)  

  ! Output summary of initial setup to user
  call sync_process()   
  call summarise_SL_input(nodalsl)                       
endif 



! Write initial (pre-looping) values to istep = 0
if(myrank.eq.0)then
  write(*,*)' Saving initial values to ensight: '
endif   
if(savedata%ice)then 
  call write_ice_to_ensight(nodalice, i_step=istep0-1)
endif 
if(savedata%sl)then
  call write_SL_to_ensight(nodalsl, i_step=istep0-1)
endif 
if(savedata%oceanf)then
  call write_OF_to_ensight(i_step=istep0-1)
endif 
call save_pot_variables(nodalphistore, nodalg, nodalB, node_valency,i_step=istep0-1)
call save_displacement_variables(strain_elmt, strain_nodal, &
                                 stress_elmt, stress_nodal, & 
                                 nodalustore, node_valency, i_step=istep0-1)








!----------------------------------------------------------------------
! ++++++++++++++++ STARTING TIME LOOPING ++++++++++++++++++++++++
! For elastic simulations there is only one timestep. 
if (myrank.eq.0)then 
  write(*,*)
  write(*,*) '----------------------------------------------'
  write(*,*) '------------ Starting time loop! -------------'
endif 



loop_step: do i_step=istep0,nstep

  call sync_process()

  if(myrank.eq.0)then   
    write(*,*)
    write(*,*)
    write(*,'(a,i0,a)')' ~~~~~~~~~~~~~~~~~~~~ TIMESTEP ', i_step, ' ~~~~~~~~~~~~~~~~~~~~'
  endif 

  ! determine time (dt) or freq (df) step and current time/freq
  call calc_time_step(i_step, t, dt, freq, ang_freq, scale_ang_freq2)
  call reset_nodal_arrays_loads(nodalu,ubcload,rhoload,nodalphi,nodalslrate)


  ! Update SL area if necessary 
  if(ISSL_DOF)then
      call update_SL_area(nodalsl, nodalu)
  endif


  ! Calculate the change in ice for this timestep
  ! and save icerate to ensight
  if(is_ICE)then 
      call sync_process()
      nodalicerate = ZERO
      call set_ice_rate(nodalice, nodalicerate, i_step)

      
      ! Calculate change in Ice mass expected
      call calculate_ice_change_volume(nodalicerate)
      call sync_process()
      call summarise_ice_vol_change()

      ! Save the Icerate 
      if(savedata%icerate)then 
        ! Save at the timestep before because this is being used to calculate THIS timestep
        call write_icerate_to_ensight(nodalicerate, i_step=i_step-1)
        
        ! But then for the final setup we need something like: 
        if(i_step.eq.nstep)then 
          call write_icerate_to_ensight(nodalicerate*zero, i_step=i_step)
        endif 
      endif 
  endif 

 
  ! Set the stiffness matrix
  call set_elasto_visco_stiffness_matrix(i_step, dt, storekmat, storemmat,& 
                                         rhoload, isscale_ang_freq, & 
                                         ang_freq, scale_ang_freq2, nelmt_viscoelas, & 
                                         eid_viscoelas, relaxtime)
                                         


  !apply traction boundary conditions for first timestep
  !WE Reads and adds the traction to the extload variable
  if((istraction.or.isfstraction).and.i_step==1)then
    call apply_traction(extload,errcode,errtag, nodalu, i_step)
    call control_error(errcode,errtag,stdout,myrank)
  endif
  ! Other types of forces: 
  if(trim(devel_example).eq.'axial_rod')then
    if(i_step>600)extload=ZERO 
  endif
  if(ismtraction)then
    call compute_magnetic_traction(errcode, errtag, extload)
  endif 
  if(iseqsource.and.eqsource_type.eq.3)then
    call compute_split_node_load(t, i_step, sfac, slipload, extload, &
                                 storekmat, errcode, errtag)
  endif 
  if(iseqsource.and.eqsource_type.lt.3.and.i_step==1)then
    call compute_cmt_load(extload, freq)
  endif 


  ! Calculate ice load: 
  if (is_ICE)then 
    call calc_ice_load(iceload, nodalicerate, nodalu, i_step=i_step)
  endif 


  ! Apply non-zero boundary conditions to the bcnodalv array 
  ! Note this is NOT applying the loading terms (e.g. extload)
  call apply_nonzero_bc(num, egdof, kmat, storekmat, bcnodalv, ubcload,&
                        nodalu, nodalphi)


  ! PREPARE BUILT IN SOLVER
  if(solver_type.eq.builtin_solver)then
    call prep_inbuilt_solver(dprecon, egdof, storekmat, ndscale, &
                             nzero_dprecon, nelmt_elas, eid_elas,&
                             nelmt_viscoelas, eid_viscoelas)
  endif


  ! Reset/initialise the count of ksp and non linear iterations 
  ksp_tot=0; nl_iter=0

  ! Some resetting of the loads 
  extload(0)=ZERO
  iceload(0)=ZERO

  ! only for fault
  load        = selfload + extload + ubcload + rhoload + iceload
  load(0)     = ZERO
  bodyload(0) = ZERO
  
  if(isplastic)then
    evpt     = ZERO
    olddu    = ZERO
    bodyload = ZERO
  endif


  ! RESETTING OF U AND DU 
  du = ZERO
  u  = ZERO

  call run_nonlinear_solver(u, du, olddu, storekmat, ndscale,           &
                            dprecon, resload, isscale_ang_freq,         &
                            ksp_iter, bodyload, load, scale_ang_freq2,  &
                            nl_iter, ksp_tot, uerr, nl_isconv, nodalu,  &
                            nodalphi, nodalslrate, bload, egdofu, dt_vp,& 
                            nelmt_elas, num, eld, eload, bmat, deriv,   &
                            eid_elas,strain_elmt, evpt, f, stress_elmt, &
                            visco_q0, elas_e0, vload,nelmt_viscoelas,   &
                            relaxtime, muratio, i_step, tratio,         &
                            eid_viscoelas, dt, q0)




  ! WE: Update our vectors with a timestep: 
  ! WE: For sea level we solve for the time derivatives of the system so then we use
  ! WE: u(t + dt) = u(t) + dt * f(t) where f is the time derivative
  ! WE: These time derivatives are stored in nodalu, nodalphi, nodalsl 

  if (ISSL_DOF)then
    
    ! TODO: MAKE FLEXIBLE DT FOR UPDATE
    dt = 1.0_kreal

    nodalustore   = nodalustore   +  dt*nodalu
    nodalphistore = nodalphistore +  dt*nodalphi
    nodalsl       = nodalsl       + (dt*nodalslrate)
    nodalice      = nodalice      + (dt*nodalicerate)

    ! Update ocean function and ocean area/volume 
    call update_ocean_function(nodalice, nodalsl, errcode, errtag)  
    call sync_process()

    ! Evaluate the ocean nodes as proportion of overall FS nodes
    totaloceannodes = sumscal(oceannodes);  
    if(myrank.eq.0)then 
      write(*,'(a, i0, a, i0)')'Total ocean nodes: ', totaloceannodes,'/', allnodesfs
      write(*,*)
    endif
  endif 

  ! Print warning if not converging
  if(nl_iter>=NL_MAXITER .and. .not.nl_isconv)then
    if(myrank==0)then
      write(*,*)
      write(*,*)'************************************************'
      write(*,'(a)')'WARNING: NON-CONVERGENCE IN NL ITERATIONS!'
      write(*,'(a, g0.4)')' --> Desired tolerance : ', NL_TOL
      write(*,'(a, g0.4)')' --> Achieved tolerance: ',uerr
      write(*,*)'************************************************'

      write(*,*)
    endif
  endif



  ! Save displacement variables to Ensight
  if(ISDISP_DOF)then
    call save_displacement_variables(strain_elmt, strain_nodal, &
                                    stress_elmt, stress_nodal, & 
                                    nodalustore, node_valency, i_step=i_step)
  endif


  ! Save potential variables to Ensight
  if(ISPOT_DOF)then
    call save_pot_variables(nodalphistore, nodalg, nodalB, node_valency,i_step=i_step)
  endif

  ! Save the ice and water to Ensight
  if(ISSL_DOF)then 
    ! Output min/max SL to stdout
    call write_min_max_SL(nodalsl)

    if(savedata%sl)then
      call write_SL_to_ensight(nodalsl, i_step)
    endif 
    if(savedata%oceanf)then
      call write_OF_to_ensight(i_step=i_step)
    endif   
    if(savedata%ice)then 
      call write_ice_to_ensight(nodalice, i_step)
    endif 
  endif 

  ! Update number of non-linear iterations
  nl_tot=nl_tot+nl_iter

enddo loop_step ! i_step time/frequency stepping loop
!----------------------------------------------------------------------
! ++++++++++++++++ END OF TIME LOOPING CODE ++++++++++++++++++++++++






! ----------------------------- CLEANUP --------------------------------
if(myrank==0)then 
  write(*,*); write(*,*)'Completed timesteps. Cleaning up... '
endif 

if(savedata%strain)then
  close(77)
  deallocate(strain_elmt,strain_nodal)
endif
! cleanup solver
if(solver_type.eq.petsc_solver)then
  call petsc_destroy_vector()                                                      
  call petsc_destroy_matrix()                                                      
  call petsc_destroy_solver()                                                      
  call petsc_finalize()
endif

call cleanup_fault()
deallocate(egdof,egdofu)
if(allocated(gdofu))deallocate(gdofu)
deallocate(extload,load,resload,rhoload,ubcload)
deallocate(du,u)
deallocate(nodalu,bcnodalv)
if(ISPOT_DOF)then
  deallocate(nodalphi)
endif

if(isplastic)then
  deallocate(olddu,evpt)
endif
if(solver_type.eq.builtin_solver .and.solver_diagscale)then
  deallocate(dprecon,ndscale)
endif
deallocate(mat_id,mat_domain,gam_blk,ym_blk,coh_blk,nu_blk,phi_blk,psi_blk,srf)
if(allocated(imat_to_imatve))deallocate(imat_to_imatve)
if(allocated(imatve_to_imat))deallocate(imatve_to_imat)
deallocate(g_coord,g_num)
deallocate(node_valency)
deallocate(bmat,deriv,eld,num)
deallocate(kmat)
if(allocated(infinite_iface))deallocate(infinite_iface)
if(allocated(infinite_face_idir))deallocate(infinite_face_idir)
call cleanup_ghost()

return
end subroutine specfem3d
!===============================================================================

