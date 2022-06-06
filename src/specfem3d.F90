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
use prestress


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
use bc
use free_surface
use benchmark
use write_ensight
use postprocess

!use relaxation_time
implicit none

character(len=250) :: myfname=' => specfem3d.f90'
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

real(kind=kreal) :: G,K
real(kind=kreal) :: cmat(nst,nst),estrain(nst),dev_strain(nst),  &
esigma(nst),sigma(nst),effsigma(nst),vsigma(nst)
real(kind=kreal) :: devp(nst),eps(nst),erate(nst),evp(nst),      &
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

! Source frequency function
real(kind=kreal) :: sff
real(kind=kreal),allocatable :: eld(:),eload(:),bload(:),   &
vload(:),rhoload(:),resload(:)
!eld: elastic displacement on all nodes of the element
!eload: elastic load on element
!extload: external load
!ubcload: load contributed by displacement BC
!load: like resload. Not currently used
real(kind=kreal),allocatable :: slipload(:),extload(:),bodyload(:), &
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
! magnetization
real(kind=kreal),allocatable :: nodalB(:,:)
!,psigma(:,:),psigma0(:,:),taumax(:),nsigma(:)
!bodyload: load computed on all nodes of the element
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

logical :: isgravity,ispseudoeq ! gravity load and pseudostatic load

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
! Old connectivity no longer necessary
deallocate(g_num0)
call sync_process



! prepare fault - UNECESSARY FOR SEA LEVEL 
if( iseqsource .and. (eqsource_type.eq.3 .or. eqsource_type.eq.4) )then
  log_msg = trim('preparing split fault...') ;   call write_ifproc0()
  call prepare_fault(errcode,errtag)
  call sync_process
  call control_error(errcode,errtag,stdout,myrank)
  log_msg = 'Complete!...' ;  call write_ifproc0() 
endif



! Apply displacement boundary conditions
allocate(bcnodalv(nndof,nnode))
bcnodalv=ZERO
log_msg = 'applying bc...'; call write_ifproc0()
allocate(gdof(nndof,nnode),stat=istat)
if (istat/=0)then
  write(logunit,*)'ERROR: cannot allocate memory!'
  flush(logunit)
  stop
endif

allocate(infinite_iface(6,nelmt),infinite_face_idir(6,nelmt))
infinite_iface=.false.
infinite_face_idir=-9999
! activate degrees of freedom
!gdof=0
!gdof(idofu,:)=1
!
!gdof=1
call activate_dof(errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)


! This will ensure that the gdof IDs are same in the finite/infinite interface
! nodes if they lie across different processors.
! This can also be done if we explicitly define the gdof ON/OFF state on those
! inteface nodes across all the processors.
call assemble_ghosts_gdof(nndof,gdof,gdof)
where(gdof>0)gdof=1
call sync_process

! At this point, all godf IDs are consistent across the parallel interfaces
! having the values either 0 or 1.

! Apply Dirichlet boundary conditions
call apply_bc(bcnodalv,errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)


!! Undo the unmatching dipalcement BCs. This may occur in fault implementation
!call sync_process
!call undo_unmatching_displacementBC(bcnodalv)
!call sync_process
call finalize_gdof(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)
log_msg = 'complete!' ; call write_ifproc0()


call modify_ghost_gdof(num, egdof, egdofu, coord, deriv, jac, bmat, &
        eld, eload, bload, vload, nodalu, nodalphi, nodalg, nodalB )


! store elemental global degrees of freedoms from nodal gdof
! this removes the repeated use of reshape later but it has larger size than
! gdof!!!

allocate(gdof_elmt(nedof,nelmt))
gdof_elmt=0
do i_elmt=1,nelmt
  gdof_elmt(:,i_elmt)=reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
enddo




!---------------------------------------
! global indexing
! this process is necessary for petsc implementation and is done in a single
! processor
call sync_process
if(ismpi .and. nproc.gt.1 .and. myrank.eq.0)call gindex()
call sync_process
!----------------------------------------------



! Output details to user 
tot_neq=sumscal(neq); max_neq=maxscal(neq); min_neq=minscal(neq)
if(myrank==0)then
  write(logunit,'(a,i0,a,i0,a,i0)')'degrees of freedoms => total:',tot_neq,&
                                  ' max:',max_neq,' min:',min_neq
  flush(logunit)
endif
log_msg = 'preprocessing...' ; call write_ifproc0()



! SPEAK TO HNG
! Calculate any prestress 
!call calc_prestress(strain_elmt, strain_nodal, stress_elmt, &
!                    stress_nodal, evpt, bodyload, slipload, &
!                    extload, viscoload, ubcload, load, du,  & 
!                    dprecon, kmat, storekmat, ksp_iter,     &
!                    errcode, errtag)
!! ALL OF THIS SHOULD BE IN PRESTRESS.f90
if(savedata%stress.or.isplastic)then
  allocate(stress_elmt(nst,ngll,nelmt),stress_nodal(nst,nnode))
  stress_elmt=ZERO
  endif
  if(savedata%strain)then
  allocate(strain_elmt(nst,ngll,nelmt),strain_nodal(nst,nnode))
  strain_elmt=ZERO
  endif
  
  if(isstress0)then
  if(s0_type==0)then
  ! compute initial stress using SEM itself
  
  allocate(extload(0:neq),du(0:neq),dprecon(0:neq), &
  storekmat(nedof,nedof,nelmt),stat=istat)
  if (istat/=0)then
  write(logunit,*)'ERROR: cannot allocate memory!'
  flush(logunit)
  stop
  endif
  extload=ZERO; isgravity=.true.; ispseudoeq=.false.
  call stiffness_bodyload(nelmt,neq,hex8_gnode,g_num,gdof_elmt,mat_id,gam_blk, &
  storekmat,dprecon,extload,isgravity,ispseudoeq)
  
  log_msg = 'complete...' ; call write_ifproc0()
  log_msg = '--------------------------------------------' ; call write_ifproc0()
  
  ! assemble from ghost partitions
  if(nproc.gt.1)then
  call assemble_ghosts(nndof,neq,dprecon,dprecon)
  endif
  
  dprecon(1:)=one/dprecon(1:); dprecon(0)=ZERO
  
  ! compute displacement due to graviy loading to compute initial stress
  du=ZERO
  call ksp_pcg_solver(neq,nelmt,storekmat,du,extload,   &
  dprecon,gdof_elmt,ksp_iter,errcode,errtag)
  call control_error(errcode,errtag,stdout,myrank)
  
  du(0)=ZERO
  
  call elastic_stress(nelmt,neq,hex8_gnode,g_num,gdof_elmt,du,stress_elmt)
  deallocate(extload,dprecon,du,storekmat)
  elseif(s0_type==1)then
  ! compute initial stress using simple relation for overburden pressure
  call overburden_stress(nelmt,g_num,z_datum,s0_datum,epk0,stress_elmt)
  else
  write(logunit,*)'ERROR: s0_type:',s0_type,' not supported!'
  flush(logunit)
  stop
  endif
  endif
!! ALL OF THIS SHOULD BE IN PRESTRESS.f90


                    




allocate(slipload(0:neq),extload(0:neq),rhoload(0:neq),ubcload(0:neq))

extload=ZERO
rhoload=ZERO

allocate(node_valency(nnode))

! compute node valency only once - needed for averaging value across interfaces
! e.g a point may be shared by 4 nodes if on the corner of 4 elements (valence 4)
! and so the stress is given as the average of these 4 
node_valency=0
do i_elmt=1,nelmt
  ielmt=i_elmt
  num=g_num(:,ielmt)
  node_valency(num)=node_valency(num)+1
enddo
! assemble all node_valceny across the processors
call assemble_ghosts_nodal_iscalar(node_valency,node_valency)


! open summary file
if(myrank==0)then
  write(logunit,'(a)')'KSP_MAXITER, KSP_TOL, NL_MAXITER, NL_TOL'
  write(logunit,'(i0,1x,g0.6,1x,i0,1x,g0.6)')KSP_MAXITER,KSP_RTOL,NL_MAXITER,NL_TOL
  !write(logunit,'(a)')'Number of SRFs'
  !write(logunit,'(i0)')nsrf
  write(logunit,'(a,i0)')'Number of time steps:',nstep
  !write(logunit,'(a)')'STEP, CGITER, NLITER, UXMAX, UMAX'
  flush(logunit)
endif





allocate(load(0:neq),bodyload(0:neq),viscoload(0:neq),             &
resload(0:neq),du(0:neq),u(0:neq),kmat(nedof,nedof),            &
storekmat(nedof,nedof,nelmt),storemmat(nedof,nelmt),stat=istat)
if(istat/=0)then
  write(logunit,*)'ERROR: cannot allocate memory!'
  flush(logunit)
  stop
endif


! Timestepping only needed for plastic/viscoelastic situations
if(isplastic)then
  allocate(olddu(0:neq),evpt(nst,ngll,nelmt))
endif

allocate(ngpart_node(nnode))

! compute stable time step for implicit integration
dt=dstep

nl_tot=0

elas_e0=ZERO
visco_q0=ZERO

nodalu=ZERO
bodyload=ZERO
viscoload=ZERO
slipload=ZERO ! slip load
extload=ZERO ! incremental external load
ubcload=ZERO
load=ZERO
u=ZERO



!! WARNING: this is temporary
!! WARNING: must remove this because it is already withing the time loop
!! earthquake source
!call earthquake_load(neq,extload,errcode,errtag)
!call sync_process
!call control_error(errcode,errtag,stdout,myrank)


! WE Prepare petsc solver - try turning into function 
! Have removed prepare_solver.f90 from Makefile for now 
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
    do i=1,nndofu
      geq=gdof(idofu(i),i_node)
      if(geq.gt.0)iseq(geq)=.true.
    enddo
  enddo
  nequ=count(iseq)
  allocate(gdofu(nequ))
  gdofu=-9999
  inum=0
  do i=1,neq
    if(iseq(i))then
      inum=inum+1
      gdofu(inum)=i
    endif
  enddo
  if(nequ.ne.inum)then
    log_msg = 'ERROR: nequ & inum mismatch!' ; call write_ifproc0()
    call close_process
  endif
  deallocate(iseq)
endif


! prepare background gravity data
! I think this only works for a global model where it does radial integration
! It can use a global model for g0 with local mesh but wont calc gravity
! just from the local mesh (always relies on some glboal model)
call prepare_gravity()


! Built-in preconditioner stuff? 
allocate(storederiv(ndim,ngll,ngll,nelmt),storejw(ngll,nelmt))
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


! This needs to be in its own file! 
! Starting time/frequency loop.
! For elastic only one timestep 
log_msg = trim(' Starting time loop!') ;   call write_ifproc0()


loop_step: do i_step=istep0,nstep
  !t=dt*real(i_step,kreal)
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

  ! Set initial values to 0 
  nodalu=ZERO
  if(ISPOT_DOF)then
    nodalphi=ZERO
  endif
  ubcload=ZERO
  !extload=ZERO
  rhoload=ZERO
  

  if(steptype.eq.FREQSTEP)then
    ! compute elastic stiffness matrix for time = 0
    if(i_step==0)then
      call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
      !if(istraction.and.trcase.eq.TRACTION_INTERNAL)then
      !  call compute_surface_stiffness(storekmat,errcode,errtag)
      !  symmetric_solver=.false.
      !  call control_error(errcode,errtag,stdout,myrank)
      !endif
    endif
    if(solver_type.eq.petsc_solver)then
      call petsc_set_stiffness_matrix_freq(storekmat,storemmat,ang_freq, &
      scale_ang_freq2,isscale_ang_freq)
      log_msg = trim(' petsc_set_stiffness_matrix: SUCCESS!') ;   call write_ifproc0()
      
      call petsc_set_ksp_operator(reuse_pc=.false.)
    endif


  else ! TIMESTEPPING not freq
    if(i_step==1)then 
      ! compute elastic stiffness matrix for time = 0
      call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
      !if(istraction.and.trcase.eq.TRACTION_INTERNAL)then
      !  call compute_surface_stiffness(storekmat,errcode,errtag)
      !  symmetric_solver=.false.
      !  call control_error(errcode,errtag,stdout,myrank)
      !endif

      ! The following is repeated below in the next elseif except for 
      ! the reuse_pc flag being true - make into a function 
      if(solver_type.eq.petsc_solver)then
        call petsc_set_stiffness_matrix(storekmat)
        log_msg = trim(' petsc_set_stiffness_matrix: SUCCESS!') ;   call write_ifproc0()
        call petsc_set_ksp_operator(reuse_pc=.false.)
        call petsc_set_solver()
      endif

    elseif(i_step==2)then
      ! Since we use a uniform dt, following routine has to be called only once 
      ! for a linear viscoelastic model. For nonlinear or nonuniform time steps
      ! it has to be called for every time steps or every changing time step.
      ! This will simply overwrite the storekmat for viscoelastic elements.
      call compute_stiffness_viscoelastic(nelmt_viscoelas,eid_viscoelas,         &
           dt,relaxtime,storekmat,errcode,errtag)

      if(solver_type.eq.petsc_solver)then
        call petsc_set_stiffness_matrix(storekmat)
        log_msg = trim(' petsc_set_stiffness_matrix: SUCCESS!') ;   call write_ifproc0()
        call petsc_set_ksp_operator(reuse_pc=.true.)
        call petsc_set_solver()
      endif
    endif

  endif ! if(steptype.eq.FREQSTEP)

  ! NOW WE ARE AT A POINT WHERE ANY STIFFNESS MATRICES HAVE BEEN CALCULATED 


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

  ! Rod example
  if(trim(devel_example).eq.'axial_rod')then
    if(i_step>600)extload=ZERO 
  endif


  ! Calculate relevant force terms: 
   
  if(ismtraction)then
    call compute_magnetic_traction(errcode, errtag, extload)
  endif 

  if(iseqsource.and.eqsource_type.eq.3)then
    call compute_split_node_load(t, i_step, sfac, slipload, extload, &
                                 storekmat, errcode, errtag)
  endif 
  
  if(iseqsource.and.eqsource_type.lt.3.and.i_step==1)then
    !call compute_moment_tensor(extload, errcode, errtag, &
     !                          freq, sff)

     ! moment-density tensor apparoch: compute equivalent moment-density tensor
    ! from the prescribe slip on the fault
    log_msg = trim(' Earthquake source type: moment-density tensor') ;   call write_ifproc0()

    call earthquake_load(neq,extload,errcode,errtag)
    call sync_process
    call control_error(errcode,errtag,stdout,myrank)

    if(steptype==FREQSTEP)then
      !WARNING: make it general for nsrc
      sff=source_frequency_function_complex(freq,source_hdur(1))
      extload=extload*sff
    endif
  endif 
  


  ! Modify RHS vector for prescribed displacements
  ! WARNING: need to check for nedofu
  do i_elmt=1,nelmt
    ielmt=i_elmt ! all elements
    num=g_num(:,ielmt)
    egdof=gdof_elmt(:,ielmt)
   
    kmat=storekmat(:,:,ielmt)
    iedof=0
    do j=1,nenode
      do i=1,nndof !nndofu
        iedof=iedof+1
        if(bcnodalv(i,num(j))/=ZERO)then
          ubcload(egdof)=ubcload(egdof)-kmat(:,iedof)*bcnodalv(i,num(j))
        endif
      enddo
    enddo
  enddo ! i_elmt
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

  extload(0)=ZERO
  ! set BC nodal displacements to nodalu array
  if(ISDISP_DOF)then
    do i=1,nndofu
      idof=idofu(i)
      do j=1,nnode
        if(bcnodalv(idof,j)/=ZERO)nodalu(i,j)=bcnodalv(idof,j)
      enddo
    enddo
  endif
  ! set BC nodal potential to nodalphi array
  if(ISPOT_DOF)then
    do i=1,nndofphi
      idof=idofphi(i)
      do j=1,nnode
        if(bcnodalv(idof,j)/=ZERO)nodalphi(j)=bcnodalv(idof,j)
      enddo
    enddo
  endif
  ksp_tot=0; nl_iter=0

  ! only for fault
  load=extload+ubcload+rhoload

  du=ZERO; u=ZERO
  
  load(0)=ZERO

  if(isplastic)then
    evpt=ZERO
    olddu=ZERO
    bodyload=ZERO
  endif
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
      write(logunit,'(a,i0,1x,e12.5,1x,e12.5)')' Residual NL: ',i_nliter, &
      maxresload,maxbodyload
      flush(logunit)
    endif

    ! starting timer
    call cpu_time(cpu_tstart)
    if(solver_type.eq.builtin_solver)then
      ! builtin solver
      if(solver_diagscale)then
        ! nondimensionlize load
        resload=ndscale*resload
        call ksp_cg_solver(neq,nelmt,storekmat,du,resload,     &
        gdof_elmt,ksp_iter,errcode,errtag)
        du=ndscale*du
        call control_error(errcode,errtag,stdout,myrank)
      else
        ! pcg solver
        call ksp_pcg_solver(neq,nelmt,storekmat,du,resload,    &
        dprecon,gdof_elmt,ksp_iter,errcode,errtag)
        call control_error(errcode,errtag,stdout,myrank)
      endif
    else
      !petsc solver
      !call petsc_set_stiffness_matrix(storekmat)
      !if(myrank==0)print*,'petsc_set_stiffness_matrix: SUCCESS!'
      if(steptype.eq.FREQSTEP.and.isscale_ang_freq)then
        resload=scale_ang_freq2*resload
      endif

      call petsc_set_vector(resload)
      log_msg = trim(' petsc_set_vector: SUCCESS!') ;   call write_ifproc0()

      call petsc_solve(du(1:),ksp_iter,ksp_convreason)
      log_msg = trim(' petsc_solve: SUCCESS!') ;   call write_ifproc0()

      continue
    endif

    call cpu_time(cpu_tend)
    telap=cpu_tend-cpu_tstart
    if(myrank==0)then
      write(format_str,*)ceiling(log10(real(telap)+1.))+5 ! 1 . and 4 decimals
      format_str='(a,1x,f'//trim(adjustl(format_str))//'.4)'
      write(logunit,fmt=format_str)' Solver Elapsed Time:',telap
    endif

    ksp_tot=ksp_tot+ksp_iter
    du(0)=ZERO
    maxdu=maxscal(maxval(abs(du)))
    if(myrank==0)then
      write(logunit,'(a,i0,1x,a,g0.6)')' KSP iters: ',ksp_iter, &
      'max du: ',maxdu
      if(solver_type.eq.petsc_solver)then
        write(logunit,'(a,i0)')' convergence reason: ',ksp_convreason
      endif
      flush(logunit)
    endif

    if(isplastic)then
      u=du
    else
      u=u+du
    endif

    maxu=maxscal(maxval(abs(u)))
    ! check convergence
    if(isplastic)then
      uerr=maxscal(maxval(abs(u-olddu)))/maxu    
      olddu=u
    else
      uerr=ZERO
      if(maxu.eq.ZERO)then
        uerr=one
      else
        uerr=maxdu/maxu
      endif
    endif
    nl_isconv=uerr.le.NL_TOL
    if(i_nliter>1.and.maxscal(maxval(abs(resload))).le.ZEROtol)nl_isconv=.true.
    if(myrank==0)then
      write(logunit,'(a,g0.6,1x,a,g0.6)')' UErr:',maxdu/maxu,'maxu:',maxu
      flush(logunit)
    endif

    call sync_process()
    ! update total nodal solution vector
    ! time steps are not incremental!!
    ! therefore, NO u(t+1)=u(t)+du
    ! u contains both diaplacement and/or gravity
    ! displacement
    if(ISDISP_DOF)then
      do i=1,nndofu
        idof=idofu(i)
        do i_node=1,nnode
          if(gdof(idof,i_node)/=0)then
            nodalu(i,i_node)=u(gdof(idof,i_node))
          endif
        enddo
      enddo
    endif
    ! gravity
    if(ISPOT_DOF)then
      do i=1,nndofphi
        idof=idofphi(i)
        do i_node=1,nnode
          if(gdof(idof,i_node)/=0)then
            ! \phi is a scalar
            nodalphi(i_node)=u(gdof(idof,i_node))
          endif
        enddo
      enddo
    endif
    bodyload=ZERO; !viscoload=ZERO

    if(ISDISP_DOF)then
      log_msg = trim('computing elemental stress') ;   call write_ifproc0()

      
      ! Compute stress
      ! Elastic elements
      ! This part is repeated for the first step. We should change this for
      ! efficiency.
      if(isplastic)bload=ZERO
      do i_elmt=1,nelmt_elas
        ielmt=eid_elas(i_elmt)
        imat=mat_id(ielmt)
        num=g_num(:,ielmt)
        egdofu=gdof_elmt(edofu,ielmt)
        eld=reshape(nodalu(:,g_num(:,ielmt)),(/nedofu/))
        bload=ZERO
        do i=1,ngll ! Integration loop
          call compute_cmat_elastic(bulkmod_elmt(i,ielmt), &
          shearmod_elmt(i,ielmt),cmat)
          
          deriv=storederiv(:,:,i,ielmt)
          jacw=storejw(i,ielmt)
        
          call compute_bmat_stress(deriv,bmat)
          estrain=matmul(bmat,eld)
          if(isplastic)then
            estrain=estrain-evpt(:,i,ielmt)
          endif
          sigma=matmul(cmat,estrain)

          if(savedata%strain)strain_elmt(:,i,ielmt)=estrain
          if(savedata%stress .and. .not.isplastic)stress_elmt(:,i,ielmt)=sigma

          if(isplastic)then
            effsigma=sigma+stress_elmt(:,i,ielmt)
            call stress_invariant(effsigma,sigm,dsbar,lode_theta)
            ! check whether yield is violated
            call mohcouf(phi_blk(imat),coh_blk(imat),sigm,dsbar,lode_theta,f)
            if(f>fmax)fmax=f

            if(f>=zero)then !.or.(nl_isconv.or.nl_iter==nl_maxiter))then
              call mohcouq(psi_blk(imat),dsbar,lode_theta,dq1,dq2,dq3)
              call formm(effsigma,m1,m2,m3)
              flow=f*(m1*dq1+m2*dq2+m3*dq3)

              erate=matmul(flow,effsigma)
              evp=erate*dt_vp
              evpt(:,i,ielmt)=evpt(:,i,ielmt)+evp
              devp=matmul(cmat,evp)
              
              ! if not converged we need body load for next iteration
              if(.not.nl_isconv .and. nl_iter/=nl_maxiter)then
                !devp(1:3)=devp(1:3)-wpressure(num(i))
                eload=matmul(devp,bmat)
                bload=bload+eload*jacw
              endif
            endif!(f>=zero)
            if(nl_isconv.or.nl_iter==nl_maxiter)then
              devp=sigma
              ! compute von Mises effective plastic strain
              !vmeps(num(i))=vmeps(num(i))+sqrt(two_third*                         &
              !dot_product(evpt(:,i,ielmt),evpt(:,i,ielmt)))
              ! update stresses
              stress_elmt(:,i,ielmt)=effsigma
              !phif_blkr=atan(tnph/srf(i_srf))
              !sf=(sigm*sin(phif_blkr)-cohf_blk*cos(phif_blkr))/&
              !(-dsbar*(cos(lode_theta)/      &
              !sqrt(r3)-sin(lode_theta)*sin(phif_blkr)/r3))
              !if(sf<scf(num(i)))scf(num(i))=sf
            endif

          else
            eload=matmul(sigma,bmat)
            bload=bload+eload*jacw
          endif
        enddo ! i=1,ngll
        if(nl_isconv .or. nl_iter==nl_maxiter)cycle
        bodyload(egdofu)=bodyload(egdofu)+bload
        !print*,maxval(abs(bload)),maxval(abs(bodyload))
      enddo ! i_elmt

      fmax=maxscal(fmax)                                                           
      if(myrank==0)then                                                            
        write(logunit,'(a,i4,a,f0.6,a,f12.6,a,f12.6)') &                           
        ' nl_iter:',nl_iter,' f_max:',fmax,' uerr:',uerr,' umax:',maxdu
        write(logunit,'(a)')'--------------------------------------------'
        flush(logunit) 
      endif

      !---------------------------------------------------------------------------
      
      if(allelastic)exit nonlinear

      ! Viscoelastic elements
      do i_elmt=1,nelmt_viscoelas
        ielmt=eid_viscoelas(i_elmt)
        imat=mat_id(ielmt)
        imatve=imat_to_imatve(imat)

        muratio=muratio_blk(:,imatve)
        tratio=dt/relaxtime(:,imatve)
        num=g_num(:,ielmt)
        egdofu=gdof_elmt(edofu,ielmt)
        eld=reshape(nodalu(:,g_num(:,ielmt)),(/nedofu/))
        bload=ZERO; vload=ZERO
        do i=1,ngll ! integration loop
          K=bulkmod_elmt(i,ielmt) 
          G=shearmod_elmt(i,ielmt)
          !call compute_cmat_maxwell(K,G,tratio,cmat)
          !call compute_cmat_elastic(K,G,cmat)
          deriv=storederiv(:,:,i,ielmt)
          jacw=storejw(i,ielmt)
        
          call compute_bmat_stress(deriv,bmat)
          estrain=matmul(bmat,eld) ! strain at current time step
          if(savedata%strain.and.i_step==1.and.i_nliter==1)then
            ! store elastic strain
            strain_elmt(:,i,ielmt)=estrain
          endif
          trace_strain=estrain(1)+estrain(2)+estrain(3)
          dev_strain(1:3)=(estrain(1:3)-ONE_THIRD*trace_strain)
          dev_strain(4:6)=estrain(4:6)*HALF
          if(savedata%stress.and.i_step==1.and.i_nliter==1)then
            ! store elastic stress
            esigma=TWO*G*dev_strain
            esigma(1:3)=esigma(1:3)+K*trace_strain
            stress_elmt(:,i,ielmt)=esigma
          endif
          !esigma=matmul(cmat,estrain)
          !esigma=TWO*G*dev_strain
          !esigma(1:3)=esigma(1:3)+K*trace_strain

          !----------------------------ZIENCKIEWICZ---------------------------
          e0=elas_e0(:,i,i_elmt)
          q0=visco_q0(:,:,i,i_elmt)
          if(i_step==1)then !.and.i_nliter==1)then
            ! initialize
            e0=dev_strain
            do i_maxwell=1,nmaxwell                                                
              q0(:,i_maxwell)=e0                                                     
            enddo
          endif
          ! compute stress and update state variables
          call compute_stress_genmaxwell(K,G,tratio,muratio,estrain,e0,q0,&
          vesigma,vsigma)
          ! compute load
          !eload=matmul(vsigma,bmat)
          !vload=vload+eload*jacw
          eload=matmul(vesigma,bmat)
          bload=bload+eload*jacw
          ! update
          if(nl_isconv.or.nl_iter==NL_MAXITER)then
            elas_e0(:,i,i_elmt)=e0
            visco_q0(:,:,i,i_elmt)=q0
            if(savedata%strain)strain_elmt(:,i,ielmt)=estrain
            if(savedata%stress)stress_elmt(:,i,ielmt)=vesigma
          endif
          !----------------------------ZIENCKIEWICZ---------------------------
        enddo ! i
        !----compute the total bodyload vector----
        !viscoload(egdofu)=viscoload(egdofu)+vload
        bodyload(egdofu)=bodyload(egdofu)+bload
      enddo ! i_elmt

      bodyload(0)=ZERO
      !viscoload(0)=ZERO
    endif !(ISDISP_DOF)

    ! time step 0  and i_nliter 0 is entirely elastic
    ! write data for tiem step 0
    if(i_step==1.and.i_nliter==1)then
      if(ISDISP_DOF)then
        ! write displacement field
        if(savedata%disp)then
          call write_vector_to_file(nnode,DIM_L*nodalu,&
          ext='dis',istep=0)
          ! On the free surface
          if(savedata%fsplot)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_L*nodalu(:,gnode_fs),&
            ext='dis',istep=0)
          endif
          if(savedata%fsplot_plane)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_L*nodalu(:,gnode_fs),&
            ext='dis',istep=0,plane=.true.)
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
          ext='sig',istep=0)
          ! On the free surface
          if(savedata%fsplot)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
            ext='sig',istep=0)
          endif
          if(savedata%fsplot_plane)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
            ext='sig',istep=0,plane=.true.)
          endif
        endif
        !if(devel_mgll)then
        !  call compute_save_density_perturbation(nodalu,errcode,errtag)
        !endif
        ! save strain at a point
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
          ext='eps',istep=0)
          ! On the free surface
          if(savedata%fsplot)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*strain_nodal(:,gnode_fs),&
            ext='eps',istep=0)
          endif
          if(savedata%fsplot_plane)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*strain_nodal(:,gnode_fs),&
            ext='eps',istep=0,plane=.true.)
          endif
          if(trim(devel_example).eq.'axial_rod')then
            write(77,*)0.0,strain_nodal(1,2099)                               
            flush(77) 
          endif
        endif
        ! Benchmark calculation for elastic result
        if(benchmark_okada .and. ISDISP_DOF)then
          call compute_okada_solution()
        endif
      endif
      if(ISPOT_DOF)then
        ! plot gravity potential
        if(savedata%gpot)then
          call write_scalar_to_file(nnode,DIM_GPOT*nodalphi,ext='gpot',istep=0) 
          ! On the free surface
          if(savedata%fsplot)then
            call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs), &
            ext='gpot',istep=0) 
          endif
          if(savedata%fsplot_plane)then
            call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs), &
            ext='gpot',istep=0,plane=.true.) 
          endif
        endif
        if(savedata%mpot)then
          call write_scalar_to_file(nnode,DIM_MPOT*nodalphi,ext='mpot',istep=0)
          ! On the free surface
          if(savedata%fsplot)then
            call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs), &
            ext='mpot',istep=0) 
          endif
          if(savedata%fsplot_plane)then
            call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs), &
            ext='mpot',istep=0,plane=.true.) 
          endif
        endif
        ! gravitational
        if(savedata%agrav)then
          ! compute acceleration due to gravity
          call compute_gradient_of_scalar(nodalphi,nodalg)
          if(nproc.gt.1)then
            call assemble_ghosts_nodal_vector(nodalg,nodalg)
          endif
          ! compute average on the sharing nodes
          do i_comp=1,ndim
            nodalg(i_comp,:)=nodalg(i_comp,:)/real(node_valency,kreal)
          enddo
          ! plot gravity accelration
          if(savedata%agrav)then
            call write_vector_to_file(nnode,DIM_G*nodalg,ext='grav',istep=0)
            ! On the free surface
            if(savedata%fsplot)then
              call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs), &
              ext='grav',istep=0)
            endif
            if(savedata%fsplot_plane)then
              call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs), &
              ext='grav',istep=0,plane=.true.)
            endif
          endif
        endif
        ! magnetic
        if(savedata%magb)then
          ! compute magnetic field
          call compute_premagnetic_field(nodalphi,nodalB)
          if(nproc.gt.1)then
            call assemble_ghosts_nodal_vector(nodalB,nodalB)
          endif
          ! compute average on the sharing nodes
          do i_comp=1,ndim
            nodalB(i_comp,:)=nodalB(i_comp,:)/real(node_valency,kreal)
          enddo
          ! multiply by \mu_0
          nodalB=MAG_CONS*nodalB
          ! plot magnetic field
          if(savedata%magb)then
            call write_vector_to_file(nnode,DIM_B*nodalB,ext='magb',istep=0)
            ! plot magnetic field on the free surface
            if(savedata%fsplot)then
              call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs), &
              ext='magb',istep=0)
            endif
            if(savedata%fsplot_plane)then
              call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs), &
              ext='magb',istep=0,plane=.true.)
            endif
          endif
        endif
      endif
      
      if(nstep.le.1.and.NL_MAXITER.le.1)then
        exit loop_step 
      endif
    endif ! i_step==1.and.i_nliter==NL_MAXITER
    ! Exit nonlinear loop if converged
    if(nl_isconv)exit nonlinear
  enddo nonlinear ! i_nliter=1,NL_MAXITER

  !!compute nodal strain--------------------------------------------------
  !if(ISDISP_DOF)then
  !  ! strain 
  !  if(savedata%strain)then
  !    if(myrank==0)then
  !      write(logunit,*)'computing nodal strain'
  !      flush(logunit)
  !    endif
  !    call compute_nodal_tensor(strain_elmt,strain_nodal)  
  !    if(nproc.gt.1)then
  !      call assemble_ghosts_nodal_vectorn(NST,strain_nodal,strain_nodal)
  !    endif
  !    ! compute average on the sharing nodes
  !    do i_comp=1,NST
  !      strain_nodal(i_comp,:)=strain_nodal(i_comp,:)/real(node_valency,kreal)
  !    enddo
  !    write(77,*)dt*real(i_step),strain_nodal(1,2099)                               
  !    flush(77)
  !  endif
  !endif
  !-----------------------------------------------------------------------------

  if(nl_iter>=NL_MAXITER .and. .not.nl_isconv)then
    if(myrank==0)then
      write(logunit,*)'WARNING: nonconvergence in nonlinear iterations!'
      write(logunit,*)'desired tolerance:',NL_TOL,' achieved tolerance:',uerr
      flush(logunit)
    endif
  endif
  nl_tot=nl_tot+nl_iter

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
  if(myrank==0)then
    write(logunit,*)' ' 
    flush(logunit)
  endif
enddo loop_step ! i_step time/frequency stepping loop
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
