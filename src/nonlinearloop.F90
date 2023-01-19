! WE - Jun 14 2022 - Module holds all subroutines needed for the 
!                    nonlinear loop 

module nonlinearloop 

contains


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
subroutine run_solver(resload, dprecon, ndscale, storekmat, du, &
                      scale_ang_freq2, ksp_iter, errcode,  &
                      ksp_convreason, errtag, isscale_ang_freq)
    ! USES
    use global 
    use set_precision
    use output_to_user

#if (USE_MPI)
use mpi_library
use math_library_mpi
use sparse
use parsolver
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



    implicit none
    
    ! IO 
    real(kind=kreal),allocatable :: resload(:), ndscale(:),dprecon(:), &
                                    storekmat(:,:,:), du(:)
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
          call control_error(errcode,errtag,stdout,myrank)
        else
          ! pcg solver
          call ksp_pcg_solver(neq,nelmt,storekmat,du,resload,    &
          dprecon,gdof_elmt,ksp_iter,errcode,errtag)
          call control_error(errcode,errtag,stdout,myrank)
        endif
      else
         ! petsc solver 
        if(steptype.eq.FREQSTEP.and.isscale_ang_freq)then
          resload=scale_ang_freq2*resload
        endif
  
        call petsc_set_vector(resload)
        log_msg=trim(' petsc_set_vector: SUCCESS!');call write_ifproc0()
  

        write(*,*)'About to run solver...iteration:'
        !call petsc_print_vector()
        !call petsc_print_matrix()

        call petsc_solve(du(1:), ksp_iter, ksp_convreason)
        log_msg = trim(' petsc_solve: SUCCESS!') ; call write_ifproc0()
  
        continue

      endif
end subroutine run_solver
!#######################################################################

subroutine check_convergence(uerr, maxu, maxdu, u, & 
                             olddu, resload, nl_isconv, i_nliter)
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
  real(kind=kreal) :: uerr,maxu,maxdu
  real(kind=kreal),allocatable :: u(:),olddu(:), resload(:)
  logical :: nl_isconv ! is there nonlinear convergence bool 
  integer :: i_nliter
  ! Local: none

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

subroutine update_nodal_u_vector( u, nodalu, nodalphi, nodalsl)
  ! USES
  use global
  use free_surface
  use math_constants
  implicit none 


  ! IO 
  real(kind=kreal),allocatable :: nodalu(:,:), u(:), nodalphi(:), nodalsl(:)
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
          nodalsl(rgnum_fs(i_gll, i_elmt)) = u(gdof(5,i_node))
        endif
      enddo !i_gll 
    enddo

    if(myrank.eq.0)then
      write(SLlogunit,*)
      write(SLlogunit,*)'Max nodal Sea Level value: ', maxval(nodalsl)
    endif 
  endif




end subroutine update_nodal_u_vector



!#######################################################################

subroutine calc_stressstrain(egdofu, nl_iter, devp, dt_vp, evp, flow,  &
                             m1, m2, m3, nelmt_elas, nl_isconv, num,   &
                             eld, eload, bload, nodalu, erate,eid_elas,&
                             cmat, estrain, bodyload, sigma, effsigma, &
                             bmat, deriv, jacw, strain_elmt, evpt ,    &
                             stress_elmt, dq1, dq2, dq3, dsbar, f,     &
                             fmax,lode_theta,sigm)

  ! USES
  use global
  use set_precision
  use plastic_library
  use math_constants
  use math_library
  use weakform,only:compute_bmat_stress
  use elastic
  implicit none 

  ! IO 
  logical :: nl_isconv
  integer :: nelmt_elas, nl_iter
  integer,allocatable :: egdofu(:), eid_elas(:), num(:)


  real(kind=kreal) :: dq1, dq2, dq3, dsbar, f, fmax, lode_theta, sigm, &
                      jacw, dt_vp, cmat(nst,nst), estrain(nst),        &
                      sigma(nst), effsigma(nst), devp(nst), erate(nst),&
                      evp(nst), m1(nst,nst), m2(nst,nst), m3(nst,nst), & 
                      flow(nst,nst)

  real(kind=kreal),allocatable :: eld(:),eload(:), bload(:), bmat(:,:),&
                                  nodalu(:,:), bodyload(:), deriv(:,:),&
                                  strain_elmt(:,:,:), &
                                  stress_elmt(:,:,:), evpt(:,:,:)

  ! Local 
  integer :: i_elmt, ielmt, imat, i_gll
  
  
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


end subroutine calc_stressstrain



!#######################################################################

! Calculate stress and strain for viscoelastic elements
subroutine visco_stressstrain(nl_iter, nl_isconv, vesigma, visco_q0,   & 
                              elas_e0, e0, i_step, bmat,  deriv, eload,& 
                              bload, K, G, strain_elmt, i_nliter,      &
                              stress_elmt, estrain,dev_strain, jacw,   &
                              trace_strain,esigma,vsigma, bodyload,    &
                              vload, imatve, nelmt_viscoelas,relaxtime,&
                              egdofu, eld, nodalu, muratio, tratio,    &
                              eid_viscoelas,dt, num, q0)
  ! USES 
  use global
  use math_constants
  use weakform,only:compute_bmat_stress
  use viscoelastic

  implicit none 

  ! IO 
  real(kind=kreal) :: trace_strain, jacw, dt, K, G, estrain(nst),      & 
                      dev_strain(nst), esigma(nst), vsigma(nst),       &  
                      e0(nst), vesigma(nst)

  real(kind=kreal),allocatable :: bmat(:,:),deriv(:,:), bodyload(:),   &
                                  visco_q0(:,:,:,:), q0(:,:),          &
                                  elas_e0(:,:,:), relaxtime(:,:),      &
                                  muratio(:), tratio(:), eld(:),       &
                                  nodalu(:,:), eload(:), bload(:),     &
                                  strain_elmt(:,:,:), vload(:),        &
                                  stress_elmt(:,:,:)

  integer :: i_step, i_nliter, imatve, nelmt_viscoelas, nl_iter
  integer,allocatable :: eid_viscoelas(:), num(:), egdofu(:)

  logical :: nl_isconv

  ! Local: 
  integer :: i_gll, i_maxwell, ielmt, imat,  i_elmt

  ! Code: 
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
    do i_gll=1,ngll ! integration loop
      K=bulkmod_elmt(i_gll,ielmt) 
      G=shearmod_elmt(i_gll,ielmt)
 
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



end subroutine visco_stressstrain




end module nonlinearloop