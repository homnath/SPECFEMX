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
use bcs_and_dof
use initialise_arrays

!use relaxation_time
implicit none

character(len=250) :: myfname=' => specfem3d.f90'
character(len=500) :: errsrc


! istat: status indicator for allocation (can be used in other contexts)
integer :: istat

! do-loop indices
integer :: i_dof,i_elmt,i_eq,i_gll,i_mat,i_nliter,i_node,i_comp,j_dof,j_node
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

real(kind=kreal) :: G,K
real(kind=kreal) :: cmat(nst,nst),estrain(nst),dev_strain(nst),  &
esigma(nst),sigma(nst),effsigma(nst),vsigma(nst)
real(kind=kreal) :: devp(nst),erate(nst),evp(nst),      &
flow(nst,nst),m1(nst,nst),m2(nst,nst),m3(nst,nst)
real(kind=kreal) :: dq1,dq2,dq3,dsbar,f,fmax,lode_theta,sigm


! cmat: elastic matric (Cijkl) in Voigt notation
! estrain: elastic strain
! esigma: elastic stress
! sigma: stress (used in different context than esigma?)
! vsigma: viscose stress

! Dynamic arrays
! num: g_num for particular element.
! node_valency: number of elements that share each node.
integer,allocatable::num(:),node_valency(:)

integer :: nzero_dprecon
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
real(kind=kreal) :: uerr,maxu,maxdu
!du: incremental solution
!u: solution (summed over du)
real(kind=kreal),allocatable :: du(:),u(:),olddu(:)

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
real(kind=kreal),allocatable :: bcnodalv(:,:),nodalu(:,:)
real(kind=kreal),allocatable :: nodalphi(:),nodalg(:,:)

! Sea level edit - WE 
real(kind=kreal),allocatable :: nodalsl(:)



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
integer :: ielmt_elas,ielmt_viscoelas,imatve,iviscoelas,nelmt_elas,            &
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

! ______________________________________________________________________
! ______________________________________________________________________


! Calculate relaxation time for viscoelastic/plastic models
allocate(relaxtime(nmaxwell,nmatblk_viscoelas),muratio(nmaxwell),tratio(nmaxwell))
call calc_relaxation_time(relaxtime, tunitfac, muratio, tratio, min_relaxtime, max_relaxtime)


!Count number of viscoelastic and elastic elements 
call count_elmts(errcode, nelmt_elas, nelmt_viscoelas, & 
mdomain, tot_nelmt_elas, max_nelmt_elas, min_nelmt_elas, tot_nelmt_viscoelas, &
max_nelmt_viscoelas, min_nelmt_viscoelas)


! Split up IDs of elastic and viscoelastic elements into new arrays
! New arrays: eid_elas and eid_viscoelas
call split_elas_visco_eids(nelmt_elas, &
nelmt_viscoelas, eid_elas, eid_viscoelas )


! Create new elastic/viscoelastic arrays
allocate(elas_e0(nst,ngll,nelmt_viscoelas), &
visco_q0(nst,nmaxwell,ngll,nelmt_viscoelas),q0(nst,nmaxwell))


! prepare ghost partitions for the communication
if(nproc.gt.1)then
  call prepare_ghost()
endif
deallocate(g_num0) ! Old connectivity no longer necessary
call sync_process


! prepare fault - Only for Fault-based simulations 
if( iseqsource .and. (eqsource_type.eq.3 .or. eqsource_type.eq.4) )then
  log_msg = trim('preparing split fault...') ;   call write_ifproc0()
  call prepare_fault(errcode,errtag)
  call sync_process
  call control_error(errcode,errtag,stdout,myrank)
  log_msg = 'Complete!...' ;  call write_ifproc0() 
endif


! Use BCs to determine global DOF index etc
call sort_gdofs_and_bc(bcnodalv, num, egdof, egdofu, coord, deriv, &
                       eld, eload, bload, vload, rhoload, resload, &
                       jac, bmat, nodalu, nodalg, nodalphi, nodalB,& 
                       tot_neq, max_neq, min_neq)





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
                            visco_q0, elas_e0, extload)



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
  log_msg = 'petsc_initialize: SUCCESS!' ; call write_ifproc0()
  
  ! Create sparse vector, matrix, and preallocate                                                         
  ! TODO: following call is not necessary for RECYCLE         
    call petsc_create_vector()                                                     
    call petsc_matrix_preallocate_size()                                           
    call petsc_create_matrix()                                                     
    call petsc_create_solver()                           
    log_msg = 'petsc_preallocate_matrix_size: SUCCESS!' ; call write_ifproc0()
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
    log_msg = 'ERROR: nequ & inum mismatch!' ; call write_ifproc0()
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


if(trim(devel_example).eq.'axial_rod')then
  open(77,file=trim(file_head)//"_strain.dat",action="write",status="replace")
endif



! Note that the stepping starts from 
!   1 for time domain.
!   0 for frequency domain.
istep0=1
if(steptype.eq.FREQSTEP)then
  istep0=0
endif


! Angular frequency
if(steptype.eq.FREQSTEP)then
  if(myrank.eq.0)then
    print*,'Step type: FREQUENCY'
    print*,'f0, f1, df (Hz):',step0,step1,dstep
  endif
endif


! Compute ELASTIC mass matrix once and for all 
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


! Initialise Sea Level 
if(is_SL)then 
  call prepare_sea_level(nodalsl)
  call set_original_sea_level()

  
  ! Save initial sea level if flagged: 
  if(savedata%sl0)then 
    call write_SL0_to_ensight()
  endif 

  if(savedata%ice0)then 
    call write_ICE0_to_ensight()
  endif 

  call update_ocean_function(u, errcode, errtag, use_s0=.true.)  ! Calc ocean func
  call calculate_SL_A()                                          ! Calc SL area 

  
  call calc_SL_LHS() 
endif 



!call close_process()




!----------------------------------------------------------------------
! ++++++++++++++++ STARTING TIME LOOPING ++++++++++++++++++++++++

! Starting time/frequency loop.
! For elastic only one timestep 
log_msg = trim(' Starting time loop!') ;   call write_ifproc0()


loop_step: do i_step=istep0,nstep
  !t=dt*real(i_step,kreal)
  

  ! _______DETERMINE TIME/FREQ STEP AND OUTPUT TO USER__________
  step=step0+dstep*real(i_step,kreal)

  if(steptype.eq.TIMESTEP)then
    ! Time step.
    t=step
    dt=dstep
    if(myrank==0)then
      write(logunit,'(a,i0,a,g0.6)')'step: ',i_step,' t: ',t
      flush(logunit)
    endif

  elseif(steptype.eq.FREQSTEP)then
    ! Frequency step.
    freq=step
    if(myrank==0)then
      write(logunit,'(a,i0,a,g0.6)')'step: ',i_step,' f: ',freq
      flush(logunit)
    endif

    if(devel_nondim)then
      ang_freq=TWO*freq*DIM_T
    else
      ang_freq=TWO*PI*freq
    endif
    scale_ang_freq2=ONE/(ang_freq*ang_freq)
  endif
  ! ________________________________________________________________________
  

  
  !  ____________________   INITIALISE VALUES   ______________________
  nodalu=ZERO
  if(ISPOT_DOF)then
    nodalphi=ZERO
  endif

  if(ISSL_DOF)then
    nodalsl = ZERO 
  endif

  ubcload=ZERO
  rhoload=ZERO
  !  ___________________  END INITIALISE VALUES   ____________________



  ! ____________________________________________________________________
  !!!!  CALCULATING ELASTIC/ LIENAR VISCOELASIC STIFFNESS MATRIX  
if(steptype.eq.FREQSTEP)then
    ! compute elastic stiffness matrix for time = 0
    if(i_step==0)then
      call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
    endif
    
    ! Set Petsc stiffness matrix
    if(solver_type.eq.petsc_solver)then
      if (ISSL_DOF)then 
        call set_petsc_stiffness_SL(isscale_ang_freq, storekmat, QSL,& 
        slc_uu, slc_pu, slc_ut, slc_up, slc_pp, slc_pt, storemmat, ang_freq, scale_ang_freq2,    &
        reuse_pc_bool=.false.,freq_bool=.true.)   
      else 
        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
        ang_freq, scale_ang_freq2, reuse_pc_bool=.false.,freq_bool=.true.)   
      endif 
    endif
else ! TIMESTEPPING not freqstepping 
  if(i_step==1)then 
    ! compute elastic stiffness matrix for time = 0
    call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
  
    if(solver_type.eq.petsc_solver)then
      if (ISSL_DOF)then 
        call set_petsc_stiffness_SL(isscale_ang_freq, storekmat, QSL, & 
        slc_uu, slc_pu, slc_ut, slc_up, slc_pp, slc_pt, storemmat, ang_freq, scale_ang_freq2,     &
        reuse_pc_bool=.false.,freq_bool=.false.)   
      else 
        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
        ang_freq, scale_ang_freq2, reuse_pc_bool=.false.,freq_bool=.false.)   
      endif 
    endif
  
  elseif(i_step==2)then
    ! Since we use a uniform dt, following routine has to be called only once 
    ! for a linear viscoelastic model. For nonlinear or nonuniform time steps
    ! it has to be called for every time steps or every changing time step.
    ! This will simply overwrite the storekmat for viscoelastic elements.
    call compute_stiffness_viscoelastic(nelmt_viscoelas,             &   
                                        eid_viscoelas, dt, relaxtime,&
                                        storekmat, errcode, errtag)
  
    if(solver_type.eq.petsc_solver)then
      if (ISSL_DOF)then 
        call set_petsc_stiffness_SL(isscale_ang_freq, storekmat, QSL,& 
        slc_uu, slc_pu, slc_ut, slc_up, slc_pp, slc_pt, storemmat, ang_freq, scale_ang_freq2,    &
        reuse_pc_bool=.true.,freq_bool=.false.)   
      else 
        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
        ang_freq, scale_ang_freq2, reuse_pc_bool=.true.,freq_bool=.false.)   
      endif ! ISDOF
    endif ! Petsc solver
  endif ! istep 
endif ! if(steptype.eq.FREQSTEP)
! ____________  FINISHED CALC. STIFFNESS MATRIX  _________________



  ! apply traction boundary conditions
  ! WARNING: i_step==1 is ONLY for rod example
  if((istraction.or.isfstraction).and.i_step==1)then
    log_msg = trim('applying traction...') ;   call write_ifproc0()
    call apply_traction(extload,errcode,errtag)
    call control_error(errcode,errtag,stdout,myrank)
    if(myrank==0)then
      write(logunit,*)'complete!',maxval(abs(extload))
      flush(logunit)
    endif
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


  ! Apply non-zero boundary conditions to the bcnodalv array 
  call apply_nonzero_bc(num, egdof, kmat, storekmat, bcnodalv, ubcload,&
                        nodalu, nodalphi)


  ! PREPARE BUILT IN SOLVER
  if(solver_type.eq.builtin_solver)then
    call prep_inbuilt_solver(dprecon, egdof, storekmat, ndscale,     &
                             nzero_dprecon, nelmt_elas, eid_elas,    &
                             nelmt_viscoelas, eid_viscoelas)
  endif


  ! Reset/initialise the count of ksp and non linear iterations 
  ksp_tot=0; nl_iter=0


  extload(0)=ZERO
  ! only for fault
  load=selfload+extload+ubcload+rhoload
  
  load(0)     = ZERO
  bodyload(0) = ZERO
  du          = ZERO
  u           = ZERO
  

  if(isplastic)then
    evpt     = ZERO
    olddu    = ZERO
    bodyload = ZERO
  endif





  ! ===================== RUN NON LINEAR ITERATIONS ====================
  !bodyload=ZERO; bodyload(0)=ZERO
  ! nonlinear iteration loop
nonlinear: do i_nliter=1,NL_MAXITER
    fmax=0 ! failure indicator for plastic simulation.
    nl_iter=nl_iter+1
   
    if(isplastic)then
      resload=load+bodyload
    else
      resload=load-bodyload
    endif

    resload(0)=ZERO
    maxresload=maxscal(maxval(abs(resload)))
    maxbodyload=maxscal(maxval(abs(bodyload)))
    
    if(myrank==0)then
      write(logunit,'(a,i0,1x,e12.5,1x,e12.5)')' Residual NL: ',i_nliter,&
      maxresload,maxbodyload
      flush(logunit)
    endif

    
    ! starting timer
    call cpu_time(cpu_tstart)

    ! Run solver for this timestep 
    call run_solver(resload, dprecon, ndscale, storekmat, du, &
                    scale_ang_freq2, ksp_iter, errcode, ksp_convreason,&
                    errtag, isscale_ang_freq)

    ! Log the time taken 
    call write_cpu_timer(format_str,cpu_tstart,cpu_tend,telap)


    ! Update ksp number and write in log file 
    ksp_tot=ksp_tot+ksp_iter
    du(0)=ZERO
    maxdu=maxscal(maxval(abs(du)))
    call log_ksp_iteration(maxdu, ksp_iter, ksp_convreason)


    if(isplastic)then
      u=du
    else
      u=u+du
    endif

    ! Get maximum value of u (disp/grav/sl etc)
    maxu=maxscal(maxval(abs(u)))

    ! check convergence
    call check_convergence(uerr, maxu, maxdu, u, & 
    olddu, resload, nl_isconv, i_nliter)

    ! Update nodal vectors following inversion step 
    call sync_process()
    call update_nodal_u_vector(u, nodalu, nodalphi, nodalsl)

    if(myrank.eq.0.and.ISSL_DOF)then
      write(SLlogunit,*)
      write(SLlogunit,*)'Max value: ', maxval(nodalsl)
    endif 

    ! Reset bodyload to ZERO for Viscoelastic iteration.
    ! We need to reconcile platic and viscoelastic iterations.
    if(.not.isplastic)bodyload=ZERO; !viscoload=ZERO


    ! Calculate the stress and strain for elastic/viscoelastic elements
    if(ISDISP_DOF)then
      log_msg = trim('computing elemental stress'); call write_ifproc0() 
      
      ! Compute stress
      ! Elastic elements
      ! This part is repeated for the first step. We should change this 
      ! for efficiency.
      if(isplastic)bload=ZERO

      ! Calculate elastic/plastic stress & strain  
      call calc_stressstrain(egdofu, nl_iter, devp, dt_vp, evp, flow,  &
                             m1, m2, m3, nelmt_elas, nl_isconv, num,   &
                             eld, eload, bload, nodalu, erate,eid_elas,&
                             cmat, estrain, bodyload, sigma, effsigma, &
                             bmat, deriv, jacw, strain_elmt, evpt ,    &
                             stress_elmt, dq1, dq2, dq3, dsbar, f,     &
                             fmax,lode_theta,sigm)

      bodyload(0)=ZERO

      fmax=maxscal(fmax)                                                           
      if(myrank==0)then                                                            
        write(logunit,'(a,f0.6)') &                           
        ' f_max:',fmax
        write(logunit,'(a)')'-------------------------------------------'
        flush(logunit) 
      endif


      ! If all elastic then leave non-linear loop.
      if(allelastic)exit nonlinear

      ! Calculate stress and strain for viscoelastic elements
      call visco_stressstrain(nl_iter, nl_isconv, vesigma, visco_q0,   &
                              elas_e0, e0, i_step, bmat, deriv, eload, &
                              bload, K, G, strain_elmt, i_nliter,      & 
                              stress_elmt, estrain, dev_strain, jacw,  &
                              trace_strain, esigma, vsigma, bodyload,  &
                              vload, imatve, nelmt_viscoelas,relaxtime,& 
                              egdofu, eld, nodalu, muratio, tratio,    &
                              eid_viscoelas, dt, num, q0)

      bodyload(0)=ZERO
      !viscoload(0)=ZERO
    endif !(ISDISP_DOF)


    
    ! time step 0  and i_nliter 0 is entirely elastic
    ! write data for tiem step 0
    if(i_step==1.and.i_nliter==1)then
      if(ISDISP_DOF)then
        call save_displacement_variables(strain_elmt, strain_nodal, &
                                         stress_elmt, stress_nodal, & 
                                         nodalu, node_valency, i_step=0)

        ! Benchmark calculation for elastic result
        if(benchmark_okada .and. ISDISP_DOF)then
          call compute_okada_solution()
        endif
      endif

      if(ISPOT_DOF)then
        call save_pot_variables(nodalphi, nodalg, nodalB, node_valency,&
                                i_step=0)
      endif
      
      if(ISSL_DOF)then 
        write(SLlogunit,*)'Saving SL to free surface'
        write(SLlogunit,*)'Saving the current SL values'
        write(SLlogunit,*)'  --> Min sea level: ', minval(DIM_L*nodalsl)
        write(SLlogunit,*)'  --> Max sea level: ', maxval(DIM_L*nodalsl)

        if(savedata%fsplot)then
          call write_scalar_to_file_freesurf(nnode_fs,  DIM_L*nodalsl, &
          ext='sl',istep=0) 
        endif
        
        if(savedata%fsplot_plane)then
          call write_scalar_to_file_freesurf(nnode_fs,  DIM_L*nodalsl, &
          ext='sl', istep=0,plane=.true.) 
        endif
      endif 



      if(nstep.le.1.and.NL_MAXITER.le.1)then
        exit loop_step 
      endif
    endif ! i_step==1.and.i_nliter==NL_MAXITER


    ! Exit nonlinear loop if converged
    if(nl_isconv)exit nonlinear

  enddo nonlinear ! i_nliter=1,NL_MAXITER
  ! ===================== FINISHED NON-LINEAR LOOP =====================


! Update and save things before the next time step starts
  ! Print warning if not converging
  if(nl_iter>=NL_MAXITER .and. .not.nl_isconv)then
    if(myrank==0)then
      write(logunit,*)'WARNING: nonconvergence in nonlinear iterations!'
      write(logunit,*)'desired tolerance:',NL_TOL,' achieved tolerance:',uerr
      flush(logunit)
    endif
  endif

  ! Update number of non-linear iterations
  nl_tot=nl_tot+nl_iter


  ! --------------- Saving/plotting variables to Ensight ---------------
  if(ISDISP_DOF)then
    call save_displacement_variables(strain_elmt, strain_nodal,    &
                                     stress_elmt, stress_nodal,    & 
                                     nodalu, node_valency, i_step)
  endif

  if(ISPOT_DOF)then
    ! plot gravity potential
    call save_pot_variables(nodalphi, nodalg, nodalB, node_valency,&
                                  i_step)
  endif
  ! ---------- Finished saving/plotting variables to Ensight -----------


  ! Add spacing to log file
  if(myrank==0)then
    write(logunit,*)' ' 
    flush(logunit)
  endif


enddo loop_step ! i_step time/frequency stepping loop
!----------------------------------------------------------------------
! ++++++++++++++++ END OF TIME LOOPING CODE ++++++++++++++++++++++++


write(*,*)'Completed timesteps. Cleaning up... '
write(SLlogunit,*)'Completed timesteps. Cleaning up... '

! ----------------------------- CLEANUP --------------------------------
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

