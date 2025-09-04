! WE - Jun 14 2022 - Module holds all subroutines needed for the 
!                    nonlinear loop 

module nonlinearloop 

contains
!_______________________________________________________________________________

!#######################################################################
subroutine write_cpu_timer(format_str, cpu_tstart,cpu_tend,telap)
! Writes CPU time to log file for solver elapsed time
! USES
use global 
use set_precision 
implicit none
! IO 
real(kind=kreal) :: cpu_tstart,cpu_tend,telap
character(len=20) :: format_str
! Local 

! Code: 

call cpu_time(cpu_tend)
telap=cpu_tend-cpu_tstart
if(myrank==0)then
  write(format_str,*)ceiling(log10(real(telap)+1.))+5 
  format_str='(a,1x,f'//trim(adjustl(format_str))//'.4)'
  write(logunit,fmt=format_str)' Solver Elapsed Time:',telap
endif

end subroutine write_cpu_timer

!#######################################################################
subroutine run_solver(scale_ang_freq2, ksp_iter, errcode,  &
                      ksp_convreason, errtag, isscale_ang_freq)
! USES
use global 
use set_precision
use shared
use output_to_user

#if (USE_MPI)
use mpi_library
use math_library_mpi
use sparse
use parsolver
#if (USE_PETSC)
#if (USE_COMPLEX)
use parsolver_petsc_complex
#else
use parsolver_petsc
#endif
#endif
#else
use serial_library
use math_library_serial
use sparse_serial
use solver
use solver_petsc
#endif
implicit none

! IO 
real(kind=kreal) :: scale_ang_freq2
integer :: ksp_iter, errcode, ksp_convreason
character(len=250) :: errtag 
logical :: isscale_ang_freq


! Code: 
if(solver_type.eq.builtin_solver)then
  ! builtin solver
  if(solver_diagscale)then
    ! nondimensionlize load
    resload=ndscale*resload
    call ksp_cg_solver(neq,nelmt,storekmat,du,resload,     &
    gdof_elmt,ksp_iter,errcode,errtag)
    du=ndscale*du
    call control_error(errcode,errtag,stdout)
  else
    ! pcg solver
    call ksp_pcg_solver(neq,nelmt,storekmat,du,resload,    &
    dprecon,gdof_elmt,ksp_iter,errcode,errtag)
    call control_error(errcode,errtag,stdout)
  endif
else
   ! petsc solver 
  if(steptype.eq.FREQSTEP.and.isscale_ang_freq)then
    resload=scale_ang_freq2*resload
  endif

  
  call petsc_set_vector(resload)
  log_msg=trim(' petsc_set_vector: SUCCESS!');call write_ifproc0

  call sync_process()

  !call petsc_print_vector()
  !call petsc_print_matrix()

  call petsc_solve(du(1:), ksp_iter, ksp_convreason)
  log_msg = trim(' petsc_solve: SUCCESS!') ; call write_ifproc0

  continue

endif
end subroutine run_solver
!#######################################################################

subroutine check_convergence(i_nliter,maxu,maxdu,nl_isconv)
! USES
use global 
use math_constants
#if (USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif 

implicit none 

! IO 
integer,intent(in) :: i_nliter
real(kind=kreal),intent(in) :: maxu,maxdu
logical,intent(out) :: nl_isconv ! is there nonlinear convergence bool 
! Local
real(kind=kreal) :: uerr

! Code: 
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
if(i_nliter>1.and.maxscal(maxval(abs(resload))).le.ZEROTOL)then 
  nl_isconv=.true.
endif 

if(myrank==0)then
 write(logunit,'(a,g0.6,1x,a,g0.6)')' UErr:',uerr,'maxu:',maxu
 flush(logunit)
endif

end subroutine check_convergence

!#######################################################################

subroutine update_nodal_u_vector(nodalslrate)
  ! USES
  use global
  use free_surface
  use math_constants
  implicit none 

  ! IO 
  real(kind=kreal),allocatable :: nodalslrate(:)
  ! local: 
  integer :: i_dof, i_node, idof, i, num(maxngll2d), i_elmt, i_gll 


  ! Code: 
  ! update total nodal solution vector
  if(ISDISP_DOF)then
    do i_dof=1,nndofu
      idof=idofu(i_dof)
      do i_node=1,nnode
        if(gdof(idof,i_node)/=0)then
          nodalu(i_dof,i_node)=u(gdof(idof,i_node))
        endif
      enddo
    enddo
  endif

  ! gravity
  if(ISPOT_DOF)then
    do i_dof=1,nndofphi
      idof=idofphi(i_dof)
      do i_node=1,nnode
        if(gdof(idof,i_node)/=0)then
          ! \phi is a scalar
          nodalphi(i_node)=u(gdof(idof,i_node))
        endif
      enddo
    enddo
  endif


  ! Sea level
  if(ISSL_DOF)then
    do i_elmt=1,nelmt_fs
      num = gnum_fs(:,i_elmt)
      do i_gll=1,maxngll2d
        i_node = num(i_gll)

        if(gdof(5, i_node)/=0)then
          nodalslrate(rgnum_fs(i_gll, i_elmt)) = u(gdof(5,i_node))
        endif
      enddo !i_gll 
    enddo

  endif

  if(myrank.eq.0)then
    write(logunit,*)'  ✓ Updated nodal u vectors. '
  endif

end subroutine update_nodal_u_vector


!#######################################################################

subroutine calc_stressstrain(nl_iter,nl_isconv,fmax)

  ! USES
  use global
  use local
  use set_precision
  use plastic_library
  use math_constants
  use math_library
  use weakform,only:compute_bmat_stress
  use elastic
  implicit none 

  ! IO 
  logical,intent(in) :: nl_isconv
  integer,intent(in) :: nl_iter
  real(kind=kreal),intent(inout) :: fmax

  ! Local 
  integer :: i_elmt, ielmt, imat, i_gll
  real(kind=kreal) :: dq1, dq2, dq3, dsbar, f, lode_theta, sigm, &
                      jacw, dt_vp, cmat(nst,nst), estrain(nst),        &
                      sigma(nst), effsigma(nst), devp(nst), erate(nst),&
                      evp(nst), m1(nst,nst), m2(nst,nst), m3(nst,nst), & 
                      flow(nst,nst)
  
  
  ! CODE: 
  ! Elastic elements
  ! This part is repeated for the first step. We should change this for
  ! efficiency.

  if(.not.isplastic)bodyload=ZERO
  do i_elmt=1,nelmt_elas
    ielmt=eid_elas(i_elmt)
    imat=mat_id(ielmt)
    num=g_num(:,ielmt)
    egdofu=gdof_elmt(edofu,ielmt)
    eld=reshape(nodalu(:,g_num(:,ielmt)),(/nedofu/))
    bload=ZERO
    do i_gll=1,ngll ! Integration loop
      call compute_cmat_elastic(bulkmod_elmt(i_gll,ielmt), &
      shearmod_elmt(i_gll,ielmt),cmat)
      
      deriv=storederiv(:,:,i_gll,ielmt)
      jacw=storejw(i_gll,ielmt)
    
      call compute_bmat_stress(deriv,bmat)
      estrain=matmul(bmat,eld)
      if(isplastic)then
        estrain=estrain-evpt(:,i_gll,ielmt)
      endif
      sigma=matmul(cmat,estrain)

      if(savedata%strain)strain_elmt(:,i_gll,ielmt)=estrain
      if(savedata%stress .and. .not.isplastic)stress_elmt(:,i_gll,ielmt)=sigma

      if(isplastic.and.isplastic_blk(imat))then
        effsigma=sigma+stress_elmt(:,i_gll,ielmt)
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
          evpt(:,i_gll,ielmt)=evpt(:,i_gll,ielmt)+evp
          devp=matmul(cmat,evp)
          
          ! if not converged we need body load for next iteration
          if(.not.nl_isconv .and. nl_iter/=nl_maxiter)then
            eload=matmul(devp,bmat)
            bload=bload+eload*jacw
          endif
        endif!(f>=zero)
        if(nl_isconv.or.nl_iter==nl_maxiter)then
          devp=sigma
          stress_elmt(:,i_gll,ielmt)=effsigma
        endif

      else ! if(isplastic)then
        eload=matmul(sigma,bmat)
        bload=bload+eload*jacw
      endif
    enddo ! i=1,ngll
    if(nl_isconv .or. nl_iter==nl_maxiter)cycle
    bodyload(egdofu)=bodyload(egdofu)+bload
  enddo ! i_elmt
  
  if(myrank.eq.0)then
    write(logunit,*)'  ✓ Calculated stress and strain '
  endif


end subroutine calc_stressstrain
!-------------------------------------------------------------------------------

! Calculate stress and strain for viscoelastic elements
subroutine visco_stressstrain(nl_iter, nl_isconv, vesigma,  & 
                              e0, i_step,& 
                              K, G, i_nliter,      &
                              estrain,dev_strain, jacw,   &
                              trace_strain,esigma,vsigma,    &
                              imatve, dt, q0)
  ! USES 
  use global
  use local
  use math_constants
  use weakform,only:compute_bmat_stress
  use viscoelastic

  implicit none 

  ! IO 
  real(kind=kreal) :: trace_strain, jacw, dt, K, G, estrain(nst),      & 
                      dev_strain(nst), esigma(nst), vsigma(nst),       &  
                      e0(nst), vesigma(nst)

  real(kind=kreal),allocatable :: q0(:,:)

  integer :: i_step, i_nliter, imatve, nl_iter
  real(kind=kreal) :: muratio(nmaxwell),tratio(nmaxwell)

  logical :: nl_isconv

  ! Local: 
  integer :: i_gll, i_maxwell, ielmt, imat,  i_elmt

  ! Code: 
  do i_elmt=1,nelmt_viscoelas
    ielmt=eid_viscoelas(i_elmt)
    imat=mat_id(ielmt)
    imatve=imat_to_imatve(imat)

    !muratio=muratio_blk(:,imatve)
    !tratio=dt/relaxtime(:,imatve)
    num=g_num(:,ielmt)
    egdofu=gdof_elmt(edofu,ielmt)
    eld=reshape(nodalu(:,g_num(:,ielmt)),(/nedofu/))
    bload=ZERO;! vload=ZERO
    do i_gll=1,ngll ! integration loop
      K=bulkmod_elmt(i_gll,ielmt) 
      G=shearmod_elmt(i_gll,ielmt)
      muratio=muratio_elmt(:,i_gll,i_elmt)
      tratio=dt/relaxtime_elmt(:,i_gll,i_elmt)
 
      deriv=storederiv(:,:,i_gll,ielmt)
      jacw=storejw(i_gll,ielmt)
    
      call compute_bmat_stress(deriv,bmat)
      estrain=matmul(bmat,eld) ! strain at current time step

      ! store elastic strain
      if(savedata%strain.and.i_step==1.and.i_nliter==1)then
        strain_elmt(:,i_gll,ielmt)=estrain
      endif

      trace_strain=estrain(1)+estrain(2)+estrain(3)
      dev_strain(1:3)=(estrain(1:3)-ONE_THIRD*trace_strain)
      dev_strain(4:6)=estrain(4:6)*HALF
      
      ! store elastic stress
      if(savedata%stress.and.i_step==1.and.i_nliter==1)then
        esigma=TWO*G*dev_strain
        esigma(1:3)=esigma(1:3)+K*trace_strain
        stress_elmt(:,i_gll,ielmt)=esigma
      endif

      !----------------------------ZIENCKIEWICZ-------------------------
      e0=elas_e0(:,i_gll,i_elmt)
      q0=visco_q0(:,:,i_gll,i_elmt)
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
        elas_e0(:,i_gll,i_elmt)=e0
        visco_q0(:,:,i_gll,i_elmt)=q0
        if(savedata%strain)strain_elmt(:,i_gll,ielmt)=estrain
        if(savedata%stress)stress_elmt(:,i_gll,ielmt)=vesigma
      endif
      !----------------------------ZIENCKIEWICZ-------------------------
    enddo ! i

    !----compute the total bodyload vector----
    bodyload(egdofu)=bodyload(egdofu)+bload

  enddo ! i_elmt

  if(myrank.eq.0)then
    write(logunit,*)'  ✓ Calculated viscoelastic stress and strain '
  endif

end subroutine visco_stressstrain
!_______________________________________________________________________________

subroutine run_nonlinear_solver(isscale_ang_freq,&
                                ksp_iter,scale_ang_freq2, nl_iter, ksp_tot, uerr,&
                                nl_isconv, nodalslrate, dt_vp, & 
                                f, &
                                i_step, dt, q0)
use global
use local
#if (USE_MPI)
use mpi_library
use ghost_library_mpi
use math_library_mpi
#else
use serial_library
use math_library_serial
use sparse_serial
use solver
use solver_petsc
#endif
use output_to_user
implicit none 

! IO variables: 
real(kind=kreal),allocatable :: nodalslrate(:), q0(:,:)
real(kind=kreal) :: uerr

logical :: isscale_ang_freq,nl_isconv
integer :: ksp_iter, nl_iter,ksp_tot, i_step
real(kind=kreal) :: scale_ang_freq2, dt_vp, f, dt


! Local variables
integer :: ksp_convreason   ! KSP convergence reason
integer :: i_nliter,imatve,i
real(kind=kreal)  :: fmax, jacw, dq1,dq2,dq3,dsbar,lode_theta,sigm, sigma(nst)
real(kind=kreal)  :: maxresload,maxbodyload
real(kind=kreal)  :: cpu_tstart,cpu_tend, telap,max_telap,mean_telap
real(kind=kreal)  :: maxu,maxdu, estrain(nst), esigma(nst)
real(kind=kreal)  :: devp(nst),evp(nst),flow(nst,nst), m1(nst,nst),& 
                     m2(nst,nst),m3(nst,nst), effsigma(nst), cmat(nst,nst),erate(nst)
real(kind=kreal) :: vesigma(nst), vsigma(nst)
real(kind=kreal) :: e0(nst) !e0: initial strain
real(kind=kreal) :: dev_strain(nst)
real(kind=kreal) :: G,K,trace_strain


character(len=20)  :: format_str
character(len=250) :: errtag ! error message
integer :: errcode
errtag=""; errcode=-1

! ===================== RUN NON LINEAR ITERATIONS ====================
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
  
  !if(myrank==0)then
  !  write(logunit,'(a,i0,1x,e12.5,1x,e12.5)')' Residual NL: ',i_nliter,&
  !  maxresload,maxbodyload
  !endif

  ! starting timer
  call cpu_time(cpu_tstart)

  ! Run NL solver for this timestep 
  ! For our purpose all this does is sets the RHS vector to be resload 
  ! And then calls the 'run' command from petsc
  ! This is one single NL iteration
  call run_solver(scale_ang_freq2, ksp_iter, errcode, ksp_convreason,&
                  errtag, isscale_ang_freq)

  ! Log the time taken 
  call write_cpu_timer(format_str,cpu_tstart,cpu_tend,telap)

  ! Update ksp number and write in log file
  ksp_tot=ksp_tot+ksp_iter
  du(0)=ZERO
  maxdu=maxscal(maxval(abs(du)))
  call log_ksp_iteration(maxdu, ksp_iter, ksp_convreason)


  ! Update the u array with du 
  ! time steps are not incremental!!
  ! therefore, NOT u(t+1)=u(t)+du
  ! u contains both diaplacement and/or gravity
  ! displacement
  if(isplastic.or.steptype.eq.FREQSTEP)then
    u=du
  else
    u=u+du
  endif

  ! Get maximum value of u (disp/grav/sl etc)
  maxu=maxscal(maxval(abs(u)))
  ! check convergence
  call check_convergence(i_nliter,maxu, maxdu, nl_isconv)

  ! Update nodal vectors following inversion step 
  call sync_process()

  ! Copy values from u --> nodalu, nodalphi etc... 
  call update_nodal_u_vector(nodalslrate)

  ! Reset bodyload to ZERO for Viscoelastic iteration.
  ! We need to reconcile platic and viscoelastic iterations.
  if(.not.isplastic)then
    if(myrank.eq.0)then
      write(*,*)'  --> set bodyload to 0'
    endif
    bodyload=ZERO; !viscoload=ZERO
  endif 

  ! Calculate the stress and strain for elastic/viscoelastic elements
  if(ISDISP_DOF)then

    ! Compute stress for Elastic elements;this part is repeated for the first step. 
    ! We should change this for efficiency.
    if(isplastic)bload=ZERO
    ! Calculate elastic/plastic stress & strain  
    call calc_stressstrain(nl_iter,nl_isconv,fmax)
    bodyload(0)=ZERO

    ! If all elastic then leave non-linear loop because only one timestep 
    !if(allelastic) then 
    !  write(logunit,*)' ALL ELASTIC --> EXITING NON LINEAR'
    !  exit nonlinear
    !endif 

    ! Calculate stress and strain for viscoelastic elements
    call visco_stressstrain(nl_iter, nl_isconv, vesigma,  &
                            e0, i_step, &
                            K, G, i_nliter,      & 
                            estrain, dev_strain, jacw,  &
                            trace_strain, esigma, vsigma, &
                            imatve, dt, q0)
    bodyload(0)=ZERO
    !viscoload(0)=ZERO
  endif !(ISDISP_DOF) 

  ! time step 0  and i_nliter 0 is entirely elastic
  ! write data for tiem step 0
  !if(i_step==1.and.i_nliter==1)then
  !write(logunit,*)'istep = 1 and i_nliter = 1'

    ! If only allowed on NL iteration and only 1 timestep then exit whole loop
    !if(nstep.le.1.and.NL_MAXITER.le.1)then
    !  write(logunit,*)'nstep.le.1.and.NL_MAXITER.le.1 - exiting loop step '
    !  exit loop_step 
    !endif
  !endif ! i_step==1.and.i_nliter==1
  ! Exit nonlinear loop if converged

  if(nl_isconv) then 
    if(myrank.eq.0)then
      write(*,*)'NL is converged...exiting non-linear '
    endif
    exit nonlinear
  endif 

enddo nonlinear ! i_nliter=1,NL_MAXITER
!======================= FINISHED NON-LINEAR ITERATIONS ========================

return 

end subroutine run_nonlinear_solver


end module nonlinearloop
!===============================================================================
