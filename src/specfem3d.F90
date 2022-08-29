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

! Source frequency function
real(kind=kreal) :: sff
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
! At this point, all gdof IDs are consistent across the parallel interfaces
! having the values either 0 or 1.


! Apply Dirichlet boundary conditions
call apply_bc(bcnodalv,errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)


!! Undo the unmatching dipalcement BCs. This may occur in fault implementation
!call sync_process
!call undo_unmatching_displacementBC(bcnodalv)
!call sync_process
! Finalise the GDOF after BCs have been applied 
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
if(ismpi .and. nproc.gt.1 .and. myrank.eq.0)then 
  call gindex()
  write(logunit,*)'CALLED GINDEX'
endif 
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




! Calculate any prestress 
call calculate_prestress(strain_elmt, strain_nodal, &
                                stress_elmt, stress_nodal, & 
                                extload, du, dprecon, storekmat, &
                                errcode, errtag, ksp_iter, istat)
   
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

! assemble all node_valency across the processors
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


! HERE IS ALLOCATION OF U VECTOR 
allocate(load(0:neq),bodyload(0:neq),selfload(0:neq),viscoload(0:neq), &
resload(0:neq),du(0:neq),u(0:neq),kmat(nedof,nedof),            &
storekmat(nedof,nedof,nelmt), storemmat(nedof,nelmt), stat=istat)

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
selfload=ZERO
viscoload=ZERO
slipload=ZERO ! slip load
extload=ZERO ! incremental external load
ubcload=ZERO
load=ZERO
u=ZERO




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
! I think this only works for a global model where it does radial integration
! It can use a global model for g0 with local mesh but wont calc gravity
! just from the local mesh (always relies on some glboal model)
call prepare_gravity()



! Built-in preconditioner stuff? 
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
! THIS SHOULD BE COMPUTE_mass_elastic_GLOBAL!!!! 
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
call prepare_sea_level()
call calc_SL_LHS(errcode, errtag)


!----------------------------------------------------------------------
! ++++++++++++++++ STARTING TIME LOOPING ++++++++++++++++++++++++

! This needs to be in its own file! 
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
  
  
  ! Set initial values to 0 
  nodalu=ZERO
  if(ISPOT_DOF)then
    nodalphi=ZERO
  endif

  if(ISSL_DOF)then
    nodalsl = ZERO ! Will still then need to allocate nodal SL values?
  endif

  ubcload=ZERO
  !extload=ZERO
  rhoload=ZERO
  


  ! ____________________________________________________________________
  !!!!  CALCULATING ELASTIC/ LIENAR VISCOELASIC STIFFNESS MATRIX  
  if(steptype.eq.FREQSTEP)then
      ! compute elastic stiffness matrix for time = 0
      if(i_step==0)then
        call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
      endif
      
      ! Set Petsc stiffness matrix
      if(solver_type.eq.petsc_solver)then
        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
                                 ang_freq, scale_ang_freq2, &
                                 reuse_pc_bool=.false.,freq_bool=.true.)            
      endif


  else ! TIMESTEPPING not freqstepping 
    if(i_step==1)then 
      ! compute elastic stiffness matrix for time = 0
      call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
    
      if(solver_type.eq.petsc_solver)then
        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
                                 ang_freq, scale_ang_freq2,            &
                                reuse_pc_bool=.false.,freq_bool=.false.) 
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
        call set_petsc_stiffness(isscale_ang_freq, storekmat,storemmat,&  
                                 ang_freq, scale_ang_freq2, &
                                 reuse_pc_bool=.true.,freq_bool=.false.) 
      endif
    endif
  endif ! if(steptype.eq.FREQSTEP)
  ! ____________________________________________________________________

  ! NOW WE ARE AT A POINT WHERE ANY STIFFNESS MATRICES HAVE BEEN 
  ! CALCULATED 

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
  

  ! ____________________________________________________________________________
  ! FINISH OFF COMPUTE_MOMENT_TENSOR subroutine
  if(iseqsource.and.eqsource_type.lt.3.and.i_step==1)then
    !call compute_moment_tensor(extload, errcode, errtag, &
     !                          freq, sff)

     ! moment-density tensor apparoch: compute equivalent moment-density tensor
    ! from the prescribe slip on the fault
    log_msg = trim(' Earthquake source type: moment-density tensor')
    call write_ifproc0()

    call earthquake_load(neq,extload,errcode,errtag)
    call sync_process
    call control_error(errcode,errtag,stdout,myrank)

    if(steptype==FREQSTEP)then
      !WARNING: make it general for nsrc
      sff=source_frequency_function_complex(freq,source_hdur(1))
      extload=extload*sff
    endif
  endif 
  ! ____________________________________________________________________________


  ! Apply non-zero boundary conditions
  call apply_nonzero_bc(num, egdof, kmat, storekmat, bcnodalv, ubcload)


  ! PREPARE BUILT IN SOLVER
  call prep_inbuilt_solver(dprecon, egdof, storekmat, ndscale,     &
                           nzero_dprecon, nelmt_elas, eid_elas,    &
                           nelmt_viscoelas, eid_viscoelas)



  ! ALSO TO DO WITH NON ZERO BOUNDARY CONDITIONS? 
  extload(0)=ZERO
  ! set BC nodal displacements to nodalu array
  if(ISDISP_DOF)then
    do i_dof=1,nndofu
      idof=idofu(i_dof)
      do j_dof=1,nnode
        if(bcnodalv(idof,j_dof)/=ZERO)nodalu(i_dof,j_dof)=bcnodalv(idof,j_dof)
      enddo
    enddo
  endif

  ! set BC nodal potential to nodalphi array
  if(ISPOT_DOF)then
    do i_dof=1,nndofphi
      idof=idofphi(i_dof)
      do j_dof=1,nnode
        if(bcnodalv(idof,j_dof)/=ZERO)nodalphi(j_dof)=bcnodalv(idof,j_dof)
      enddo
    enddo
  endif


  ! ____________________________________________________________
  ! NEED TO ADD IN INITIAL SL SETUP AND ADD TO NODALSL HERE


  ! ____________________________________________________________


  ! Reset/initialise the count of ksp and non linear iterations 
  ksp_tot=0; nl_iter=0

  ! only for fault
  load=selfload+ubcload+rhoload

  du=ZERO; u=ZERO
  
  load(0)=ZERO

  if(isplastic)then
    evpt=ZERO
    olddu=ZERO
    bodyload=ZERO
  endif

  bodyload(0)=ZERO


  ! ===================== RUN NON LINEAR ITERATIONS =====================
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

    ! Run the solver! 
    call run_solver(resload, dprecon, ndscale, storekmat, du, &
                    scale_ang_freq2, ksp_iter, errcode, ksp_convreason, &
                    errtag, isscale_ang_freq)
    
    ! Output time details
    call write_cpu_timer(format_str,cpu_tstart,cpu_tend,telap)


    ! Update total iterations and max values for outputting 
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

    maxu=maxscal(maxval(abs(u))) ! Get max across all processors. 
    
    call check_convergence(uerr, maxu, maxdu, u, & 
                           olddu, resload, nl_isconv, i_nliter)
    call sync_process()


    call update_nodal_u_vector(u, nodalu, nodalphi, nodalsl)
    

    ! Reset bodyload to ZERO for viscoelastic iteration
    ! We need to reconcile plastic and viscoelastic iterations. 
    if(.not.isplastic)bodyload=ZERO; !viscoload=ZERO


    if(ISDISP_DOF)then
      log_msg = trim('computing elemental stress') ;   call write_ifproc0() 
     
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
        write(logunit,'(a,i4,a,f0.6,a,f12.6,a,f12.6)') &                           
        ' nl_iter:',nl_iter,' f_max:',fmax,' uerr:',uerr,' umax:',maxdu
        write(logunit,'(a)')'--------------------------------------------'
        flush(logunit) 
      endif

      ! If only elastic elements then only need 1 iteration and no 
      ! viscoelastic calculations needed so exit entire loop      
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
    endif !(ISDISP_DOF)







    ! CONVERT TO SUBROUTINE - have we already written subroutines that 
    ! can be called for some of these? 
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

  ! =========================END NON LINEAR ITERATIONS LOOP =============================


  ! Warning if there is not convergence 
  if(nl_iter>=NL_MAXITER .and. .not.nl_isconv)then
    if(myrank==0)then
      write(logunit,*)'WARNING: nonconvergence in nonlinear iterations!'
      write(logunit,*)'desired tolerance:',NL_TOL,' achieved tolerance:',uerr
      flush(logunit)
    endif
  endif
  nl_tot=nl_tot+nl_iter


  ! ---------------- Save variables to EnSight: -------------------
  if(ISDISP_DOF)then
    call save_displacement_variables(strain_elmt, strain_nodal, &
                                     stress_elmt, stress_nodal, & 
                                     nodalu, node_valency, i_step)
  endif

  if(ISPOT_DOF)then
    call save_potential_variables(nodalphi,nodalg, nodalB, &
                                  node_valency, i_step)
  endif
  ! -------------------------------------------------------------
  
  ! Seperate loop steps in log file
  if(myrank==0)then 
    write(logunit,*)' '
     flush(logunit)
  endif
  
enddo loop_step ! i_step time/frequency stepping loop





! cleanup solver
call run_cleanup_specfem3d(strain_elmt, strain_nodal,        &
                           evpt, egdof, egdofu, gdofu,       &
                           extload, ubcload, load, du, u,    &
                           olddu, bcnodalv, nodalu, nodalphi,&
                           dprecon, ndscale, num,            &
                           node_valency, bmat, deriv, kmat,  &
                           rhoload, resload, eld)

return
end subroutine specfem3d
!===============================================================================
