! This module contains preprocessing library routines.
! REVISION:
!   HNG, Jul 07,2011
module matrix_vector
    implicit none
    character(len=250),private :: myfname=" => matrix_vector.f90"
    character(len=500),private :: errsrc
    
    contains
    !-------------------------------------------------------------------------------
    
    ! This subrotine computes the stiffness matrix, diagonal preconditioner, and
    ! body loads (gravity and pseudostatic loads) optionally
    ! TODO: optional precoditioner,optional assembly of stiffness
    subroutine compute_bodyload(bodyload,selfweight,pseudoeq)
    use set_precision
    use nondimensionpar
    use global,only:myrank,NDIM,nst,nelmt,ngll,nedof,nedofu,nedofphi,nenode,ngnode,&
    ngllx,nglly,ngllz,ngll,g_coord,gdof_elmt,g_num,mat_domain,mat_id,agrav,massdens_elmt,&
    bulkmod_elmt,shearmod_elmt,isempty_blk,rho_blk,ym_blk,magnetization_elmt,&
    infinite_iface,infinite_face_idir,pole_coord0,pole_coord1,       &
    pole_type,pole_axis,axis_range,ISDISP_DOF,ISPOT_DOF,storederiv,    &
    POT_TYPE,PGRAVITY,PMAGNETIC,    &
    element_is_infinite,storejw,storeinterpf_infinite,devel_nondim, &
    isdxval,isdyval,isdzval,devel_gaminf,infquad, &
    edofu,edofphi,grav0_nodal,dgrav0_elmt,ISGRAV0,&
    imat_to_imatmag,magnetization_blk,ismagnet_blk,eqkx,eqky,eqkz
    use elastic,only:compute_cmat_elastic
    use math_library,only:determinant,invert
    use weakform,only:compute_bmat_stress
    use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights
    implicit none
    ! only intact elements
    real(kind=kreal),intent(inout) :: bodyload(0:)
    logical,intent(in),optional :: selfweight,pseudoeq
    
    real(kind=kreal) :: zero=0.0_kreal
    integer :: egdof(nedof),num(nenode)
    integer :: i,i_elmt,ielmt,k
    real(kind=kreal) :: density,g ! acceleration due to gravity
    ! jacw=jacobian*weight
    real(kind=kreal) :: jacw
    real(kind=kreal) :: interpf(NGLL)
    real(kind=kreal) :: eload(nedof),eqload(nedof)
    integer :: nip,nipinf,ielmt_infinite
    logical :: isinf ! flag to check if the element if on the infinite domain
    logical :: isfaces(6)
    
    if(.not.present(selfweight).and. .not.present(pseudoeq))then
      write(*,'(a)')'WARNING: both "selfweight" and "pseudoeq" are absent &
      &NO body load computed!'
      stop
    endif
    if(.not.selfweight.and. .not.pseudoeq)then
      write(*,'(a)')'WARNING: both "selfweight" and "pseudoeq" are switched OFF &
      &NO body load computed!'
      stop
    endif
    
    ! set appropriate acceleration due to gravity
    g = agrav
    if(devel_nondim)then
     g = g*NONDIM_ACCEL
    endif
    ! counter for infinite elements
    ielmt_infinite=0
    ! compute body load
    do i_elmt=1,nelmt
      ielmt=i_elmt
      isinf=element_is_infinite(ielmt)
      if(isinf)ielmt_infinite=ielmt_infinite+1
      num=g_num(:,i_elmt)
      egdof=gdof_elmt(:,i_elmt)
      eload=zero; eqload=zero
      do i=1,ngll
        density=massdens_elmt(i,i_elmt)
        if(isinf)then
          ! infinite element
          interpf=storeinterpf_infinite(i,:,ielmt_infinite) !lagrange_gl(i,:)
        else !(isinf)
          ! standard element
          interpf=lagrange_gll(i,:)
        endif
    
        jacw=storejw(i,ielmt)
        
        eload(3:nedof:3)=eload(3:nedof:3)+interpf*jacw*density
      enddo ! i=1,ngll
    
      ! compute body loads
      ! gravity load and add to extload
      ! WARNING: this must be changed to massdens_elmt
      if(present(selfweight))then
        if(selfweight)then
          bodyload(egdof)=bodyload(egdof)-eload*g
        endif
      endif
    
      if(present(pseudoeq))then
        if(pseudoeq)then
          ! compute pseudostatic earthquake loads and add to extload
          eqload(1:nedof:3)=eqkx*eload(3:nedof:3)
          eqload(2:nedof:3)=eqky*eload(3:nedof:3)
          eqload(3:nedof:3)=eqkz*eload(3:nedof:3)
          bodyload(egdof)=bodyload(egdof)+eqload*g
        endif
      endif
    enddo ! i_elmt=1,nelmt
    bodyload(0)=zero
    end subroutine compute_bodyload
    !===============================================================================
    

    subroutine compute_storekmatSL(storekmatSL, kSL)
      ! This subroutine computes the contribution to the stiffness matrix from the SL 
      ! and stores it in storekmatSL
      use math_constants
      use free_surface
      ! IO variables: 
      real(kind=kreal) :: storekmatSL(:,:,:), kSL(:,:)

      ! local variables: 
      integer :: i_elmtfs

    ! Sea level free surface contributions to stiffness matrix:
      ! Initialise
      storekmatSL = zero
      ! Loop for each face on the free surface (note some elements, those with more than one face on the FS)
      ! will be considered multiple times, but for different dofs.
      do i_elmtfs = 1, nelmt_fs
        kSL = zero 
        call calc_SL_stiffness(i_elmtfs, kSL)
        storekmatSL(:,:,id_elem_fs(i_elmtfs)) = storekmatSL(:,:,id_elem_fs(i_elmtfs)) +  kSL
      enddo 

      end subroutine compute_storekmatSL
!===============================================================================

    ! This subrotine computes the stiffness matrix, and
    ! body loads contributed by mass density or magnetiztion.
    ! TODO: optional precoditioner,optional assembly of stiffness
    subroutine compute_stiffness_elastic(errcode,errtag)
    use set_precision
    use global
    use element 
    use elastic,only:compute_cmat_elastic
    use math_constants,only:HALF,ONE,ZERO,FOUR,GRAV_CONS,VACUUM_PERMITTIVITY,PI
    use math_library,only:determinant,invert,issymmetric
    use weakform
    use shape_library
    !use sea_level
    use gll_library
    use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights, prepare_integration
    use infinite_element
    use free_surface, only: nelmt_fs, iface_fs, id_elem_fs

    !use ieee_arithmetic
    implicit none
    integer,intent(out) :: errcode
    character(len=250),intent(out) :: errtag
    real(kind=kreal),parameter :: FOUR_PI_G=FOUR*PI*GRAV_CONS
    real(kind=kreal),parameter :: FOUR_PI_G_INV=ONE/FOUR_PI_G
    integer :: i,i_gll
    integer :: i_elmt,ielmt,imat
    integer :: imatve,mdomain
    integer :: num(nenode),egdof(nedof)
    real(kind=kreal) :: cmat(nst,nst)
    real(kind=kreal) :: detjac !determinant of Jacobian
    real(kind=kreal) :: coord(ngnode,NDIM),jac(NDIM,NDIM),eload(nedof)
    real(kind=kreal) :: dgmat(NDIM,NDIM),eg0(ngll,NDIM),g0(NDIM),dg0(6)
    
    real(kind=kreal) :: kmat(nedof,nedof),xval
    ! jacw=jacobian*weight
    real(kind=kreal) :: jacw
    real(kind=kreal) :: interpf(NGLL),deriv(NDIM,nenode)
    
    ! magnetization
    integer :: imatmag
    real(kind=kreal) :: M(NDIM),Mgll(NDIM,ngll)
    real(kind=kreal) :: divM,dxMx,dyMy,dzMz
   
    ! electrical conductivity
    real(kind=kreal) :: econgll(ngll)

    real(kind=kreal) :: bmatu(NST,NEDOFU),wmat_gradphi(NEDOFU,NDIM),      &
    rmat_gradphi(NDIM,NEDOFPHI),wmat_sPE(NEDOFPHI,NDIM),rmat_sPE(NDIM,NEDOFU), &
    rmat_term1(NDIM,NEDOFU),wmat_term1(NEDOFU,NDIM),                           &
    rmat_term2(NDIM,NEDOFU),wmat_term2(NEDOFU,NDIM),                           &
    rmat_term3(1,NEDOFU),wmat_term3(NEDOFU,1),                                 &
    rmat_term4(1,NEDOFU),wmat_term4(NEDOFU,1)
    
    ! for infinite elements
    integer,parameter :: nginf=8
    integer :: i_face,idir, iloop 
    real(kind=kreal) :: tcoord
    real(kind=kreal) :: gaminf
    real(kind=kreal) :: polex(4,NDIM)
    real(kind=kreal) :: coordinf(nginf,NDIM)
    real(kind=kreal),allocatable :: shape_infinite(:,:),dshape_infinite(:,:,:)
    real(kind=kreal),allocatable :: lagrange_gl(:,:),dlagrange_gl(:,:,:)
    real(kind=kreal),allocatable :: GLw(:)
    integer :: nip,nipinf,ielmt_infinite
    logical :: isinf ! flag to check if the element if on the infinite domain
    logical :: isfaces(6)
    logical :: empty_row 

    integer :: emp, filled 
    
    integer :: i_elmtfs, ios


    errtag="ERROR: unknown!"
    errcode=-1
    errsrc=trim(myfname)//' => compute_stiffness_elastic'
    
    storekmat=zero
    rhoload=zero
   

    ! CALCULATE THE REST OF MATRIX
    if(myrank.eq.0)then 
      write(*,*)' --> Calculating normal stiffness matrix'
    endif 
    !storekmat=zero
    rhoload=zero
    ! Purely elastic elements
    ! Viscoelastic elements are elastic at time = 0
    ! Following loops through nelmt_elas+nelmt_viscoelas
    ielmt_infinite=0
    do i_elmt=1,nelmt
      ielmt=i_elmt
      isinf=element_is_infinite(ielmt)
      if(isinf)ielmt_infinite=ielmt_infinite+1
      num=g_num(:,ielmt)
      imat=mat_id(ielmt)
      mdomain=mat_domain(imat) 
      nip=ngll
    
      ! set magnetization at GLL points
      Mgll=ZERO
      if(POT_TYPE==PMAGNETIC)then
        if(ismagnet_blk(imat))then
          do i_gll=1,ngll 
            Mgll(:,i_gll)=magnetization_elmt(:,i_gll,ielmt)
          enddo
        endif
      endif
      ! set electrical conductivity at GLL points
      econgll=ZERO
      if(POT_TYPE==PELECTRIC)then
        if(iselectric_blk(imat))then
          do i_gll=1,ngll
            econgll(i_gll)=econductivity_elmt(i_gll,ielmt)
          enddo
        endif
      endif

      !H=dgrav0_elmt(:,:,ielmt)
      eg0=transpose(grav0_nodal(:,num))
    
      egdof=gdof_elmt(:,i_elmt)
        
      kmat=zero
      eload=zero
      do i=1,nip
        if(ISDISP_DOF)then
          call compute_cmat_elastic(bulkmod_elmt(i,ielmt),shearmod_elmt(i,ielmt),  &
          cmat)
        endif
    
        if(isinf)then
          ! infinite element
          interpf=storeinterpf_infinite(i,:,ielmt_infinite) !lagrange_gl(i,:)
        else !(isinf)
          ! standard element
          interpf=lagrange_gll(i,:)
        endif
    
        deriv=storederiv(:,:,i,ielmt) 
        jacw=storejw(i,ielmt)
        
        if(ISDISP_DOF)then
          ! compute only for nonempty elements
          if(.not.isempty_blk(imat))then
            call compute_bmat_stress(deriv,bmatu)
            kmat(edofu,edofu)=kmat(edofu,edofu)+matmul(matmul(transpose(bmatu),cmat),bmatu)*jacw
             
            if(ISGRAV0)then
              g0=eg0(i,:)
              ! compute dg0=\nabla g0
              dgmat=matmul(deriv,eg0)
              dg0(1)=dgmat(1,1)
              dg0(2)=dgmat(2,2)
              dg0(3)=dgmat(3,3)
              dg0(4)=dgmat(1,2)
              dg0(5)=dgmat(1,3)
              dg0(6)=dgmat(2,3)
    
              call compute_rmat_term1(massdens_elmt(i,ielmt),interpf,rmat_term1)
              call compute_wmat_term1(g0,dg0,interpf,deriv,wmat_term1)
              call compute_rmat_term2(g0,dg0,interpf,deriv,rmat_term2)
              call compute_wmat_term2(massdens_elmt(i,ielmt),interpf,wmat_term2)
              kmat(edofu,edofu)=kmat(edofu,edofu)-HALF*(matmul(wmat_term1,rmat_term1)+ &
              matmul(wmat_term2,rmat_term2))*jacw
      
              call compute_rmat_term3(massdens_elmt(i,ielmt),g0,interpf,rmat_term3)
              call compute_wmat_term3(deriv,wmat_term3)
              call compute_rmat_term4(massdens_elmt(i,ielmt),deriv,rmat_term4)
              call compute_wmat_term4(g0,interpf,wmat_term4)
              kmat(edofu,edofu)=kmat(edofu,edofu)+HALF*(matmul(wmat_term3,rmat_term3)+ &
              matmul(wmat_term4,rmat_term4))*jacw
            endif
            
            if(ISPOT_DOF)then
              ! w.rho*grad(phi)
              call compute_wmat_gradphi(interpf,wmat_gradphi)
              call compute_rmat_gradphi(massdens_elmt(i,ielmt),deriv,rmat_gradphi)          
              kmat(edofu,edofphi)=kmat(edofu,edofphi)+matmul(wmat_gradphi,rmat_gradphi)*jacw
              ! grad(w).rho*s
              ! This term is a transpose of the previous. Therefore, it will be
              ! computed later
              call compute_wmat_sPE(deriv,wmat_sPE)
              call compute_rmat_sPE(massdens_elmt(i,ielmt),interpf,rmat_sPE)
              kmat(edofphi,edofu)=kmat(edofphi,edofu)+matmul(wmat_sPE,rmat_sPE)*jacw
            endif
          endif
        endif
    
        if(ISPOT_DOF)then
          if(POT_TYPE==PELECTRIC)then
            kmat(edofphi,edofphi)=kmat(edofphi,edofphi)+matmul(transpose(deriv),deriv)*econgll(i)*jacw
          else
            kmat(edofphi,edofphi)=kmat(edofphi,edofphi)+matmul(transpose(deriv),deriv)*jacw
          endif

          if(.not.ISDISP_DOF)then
            ! Magnetic
            if(POT_TYPE==PMAGNETIC)then
              ! Code segment below is not taken into account. Magnetization is 
              ! implemented via magnetic traction. See apply_mtraction.f90.
              if(ismagnet_blk(imat))then
                ! first compute the divergence of M
                ! \nabla.M=dxMx+dyMy+dzMz
                dxMx=dot_product(deriv(1,:),Mgll(1,:))
                dyMy=dot_product(deriv(2,:),Mgll(2,:))
                dzMz=dot_product(deriv(3,:),Mgll(3,:))
                divM=dxMx+dyMy+dzMz
                if(maxval(abs(M)).gt.ZERO)then
                  eload(edofphi)=eload(edofphi)+lagrange_gll(i,:)*divM*jacw
                endif
              endif
            endif
            ! Gravity
            if(POT_TYPE==PGRAVITY)then
              eload(edofphi)=eload(edofphi)+lagrange_gll(i,:)*massdens_elmt(i,ielmt)*jacw
            endif
            ! Charge density
            if(POT_TYPE==PCHARGE)then
              eload(edofphi)=eload(edofphi)+lagrange_gll(i,:)*charge_density_elmt(i,ielmt)*jacw
            endif
          endif
        endif
    
      enddo !nip

      ! kmat terms for grad(w).rho*s
      !kmat(edofphi,edofu)=transpose(kmat(edofu,edofphi))
      if(ISDISP_DOF.and.ISPOT_DOF)then
        if(devel_nondim)then
          ! Note: PI*G is nondimensionalized
          kmat(edofphi,edofphi)=0.25_kreal*kmat(edofphi,edofphi)
        else
          kmat(edofphi,edofphi)=FOUR_PI_G_INV*kmat(edofphi,edofphi)
        endif
      endif


      if(.not.issymmetric(kmat))then
        write(*,*)'ERROR: matrix is unsymmetric!'
        stop
      endif
  
      storekmat(:,:,ielmt)=kmat
      if(.not.ISDISP_DOF .and. ISPOT_DOF)then
        if(POT_TYPE==PGRAVITY .or. POT_TYPE==PCHARGE)then
          rhoload(egdof)=rhoload(egdof)+eload
        endif
      endif
    enddo ! i_elmt

   
    ! rhoload is computed if only the ISPOT_DOF is TRUE
    ! multiply rhoload by 4*PI*G
    if(.not.ISDISP_DOF.and.ISPOT_DOF)then
      if(POT_TYPE==PGRAVITY)then
        if(.not.devel_nondim)then
          rhoload=FOUR_PI_G*rhoload
        else
          rhoload=FOUR*rhoload
          ! Note: PI*G is nondimensionalized
        endif
       elseif(POT_TYPE==PCHARGE)then
        if(.not.devel_nondim)then
          rhoload=(-ONE/VACUUM_PERMITTIVITY)*rhoload
        else
          rhoload=FOUR*rhoload
          ! Note: PI*G is nondimensionalized
        endif
       endif
      rhoload(0)=ZERO
    endif
    
    end subroutine compute_stiffness_elastic
    !===============================================================================
    
    ! This subroutine computes the free surafce contribution
    ! on the stiffness matrix.
    subroutine compute_surface_stiffness(errcode,errtag)
    use global
    use math_constants
    use element,only:hexface,hexface_sign
    use integration,only:dshape_quad4_xy,dshape_quad4_yz,dshape_quad4_zx,          &
                         gll_weights_xy,gll_weights_yz,gll_weights_zx,             &
                         lagrange_gll_xy,lagrange_gll_yz,lagrange_gll_zx
    use dof,only:set_face_vecdof,set_face_scaldof
    use weakform,only:compute_rmat_sn
    implicit none
    integer,intent(out) :: errcode
    character(len=250),intent(out) :: errtag
    
    integer :: i_face,i_gll,ios
    integer :: ielmt,iface,nface
    integer :: num(nenode)
    integer :: nfdofu,nfdofphi,nfdof,nfgll
    !face nodal dof, face gll points
    integer :: tractype,count_trac
    logical :: trac_stat
    real(kind=kreal) :: coord(NDIM,4)
    real(kind=kreal) :: detjac
    real(kind=kreal),dimension(NDIM) :: face_normal,dx_dxi,dx_deta
    
    integer,allocatable :: edof(:),imapuf(:),imapphif(:)
    real(kind=kreal),dimension(:,:,:),allocatable :: dshape_quad4
    real(kind=kreal),dimension(:),allocatable :: gll_weights
    real(kind=kreal),dimension(:,:),allocatable :: lagrange_gll
    real(kind=kreal),dimension(:,:,:),allocatable :: dlagrange_gll
    real(kind=kreal),allocatable :: interpf(:)                                       
    real(kind=kreal),allocatable :: wmat_sn(:,:),rmat_sn(:,:),kmat(:,:),rhofgll(:)
    
    character(len=80) :: fname
    character(len=80) :: data_path
    
    errtag="ERROR: unknown!"
    errcode=-1
    errsrc=trim(myfname)//' => compute_surface_stiffness'
    
    ! set data path
    if(ismpi.and.nproc.gt.1)then
      data_path=trim(part_path)
    else
      data_path=trim(inp_path)
    endif
    
    fname=trim(data_path)//trim(trfile)//trim(ptail_inp)
    open(unit=11,file=trim(fname),status='old',action='read',iostat=ios)
    if (ios /= 0)then
      write(errtag,'(a)')'ERROR: input file "'//trim(fname)//'" cannot be opened!'
      return
    endif
    
    allocate(dshape_quad4(2,4,maxngll2d))
    allocate(gll_weights(maxngll2d),lagrange_gll(maxngll2d,maxngll2d),             &
    dlagrange_gll(2,maxngll2d,maxngll2d))
    
    trac_stat=.true. ! necessary for empty trfile
    count_trac=0
    traction: do
      read(11,*,iostat=ios)tractype
      if(ios/=0)exit traction
      count_trac=count_trac+1
      trac_stat=.false.
    
      read(11,*) ! skip this line, we do not need
      read(11,*)nface

      do i_face=1,nface
        read(11,*)ielmt,iface
        if(iface==1 .or. iface==3)then
          ! ZX faces
          nfgll=ngllzx
          lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_zx
          gll_weights(1:nfgll)=gll_weights_zx
          dshape_quad4(:,:,1:nfgll)=dshape_quad4_zx
        elseif(iface==2 .or. iface==4)then
          ! YZ faces
          nfgll=ngllyz
          lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_yz
          gll_weights(1:nfgll)=gll_weights_yz
          dshape_quad4(:,:,1:nfgll)=dshape_quad4_yz
        elseif(iface==5 .or. iface==6)then
          ! XY faces
          nfgll=ngllxy
          lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
          gll_weights(1:nfgll)=gll_weights_xy
          dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
        else
          write(errtag,'(a)')'ERROR: wrong face ID for traction!'
          exit traction
        endif
        nfdof=nfgll*nndof
        nfdofu=nfgll*NNDOFU                                                            
        nfdofphi=nfgll*NNDOFPHI                       
        allocate(edof(nfgll))
        allocate(rhofgll(nfgll))
        allocate(interpf(nfgll))
        allocate(kmat(nfdof,nfdof))                                                   
        allocate(imapuf(nfdofu),imapphif(nfdofphi))                                        
        allocate(wmat_sn(nfdofphi,1),rmat_sn(1,nfdofu))                                  
        call set_face_vecdof(nfgll,3,1,imapuf) !ndofphi=1
        call set_face_scaldof(nfgll,4,0,imapphif) !idofphi=4
      
        print*,'******************************* '
        print*,'Element, face: ', ielmt,iface
        print*,'nfdof:',nfdof
        print*,'nfdofu:',nfdofu
        print*,'nfdofphi:',nfdofphi
        print*,'imapuf:',imapuf
        print*,'imapphif:',imapphif
        print*,'iface:',iface
        print*,'edof:',hexface(iface)%edof
        print*,' '

        edof=hexface(iface)%edof
        num=g_num(:,ielmt)
        coord=g_coord(:,num(hexface(iface)%gnode))
        rhofgll=massdens_elmt(hexface(iface)%node,ielmt)
        kmat=zero
        ! compute numerical integration
        do i_gll=1,nfgll
    
          ! compute two vectors dx_dxi and dx_deta
          dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
          dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))
    
          ! Normal = (dx_dxi x dx_deta)
          face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
          face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
          face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
    
          detjac=sqrt(dot_product(face_normal,face_normal))
          face_normal=hexface_sign(iface)*face_normal/detjac
          print*,face_normal
          interpf=lagrange_gll(i_gll,:)
          wmat_sn(:,1)=interpf
          ! -w_\phi \rho s.n                                                                      
          call compute_rmat_sn(nfgll,face_normal,interpf,rmat_sn)                           
          kmat(imapphif,imapuf)=kmat(imapphif,imapuf)-matmul(wmat_sn,rmat_sn)* &
                                rhofgll(i_gll)*detjac*gll_weights(i_gll)
        enddo ! i_gll
        !kmat(imapuf,imapphif)=transpose(kmat(imapphif,imapuf))
        storekmat(edof,edof,ielmt)= storekmat(edof,edof,ielmt)+kmat!*1e0
        print*,minval(abs(kmat)),maxval(abs(kmat)),minval(rhofgll),maxval(rhofgll)
        deallocate(edof)                                                   
        deallocate(rhofgll)                                                   
        deallocate(interpf)
        deallocate(kmat)                                                   
        deallocate(imapuf,imapphif)                                        
        deallocate(wmat_sn,rmat_sn)                                  
      enddo ! i_face
      trac_stat=.true.
    enddo traction
    
    close(11)
    deallocate(dshape_quad4)
    deallocate(gll_weights,lagrange_gll,dlagrange_gll)
    !deallocate(kmat,wmat_sn,rmat_sn)
    !deallocate(imapuf,imapphif)
    if(.not.trac_stat)then
      write(errtag,'(a)')'ERROR: all tractions cannot be read!'
      return
    endif
    
    errcode=0
    
    return
    end subroutine compute_surface_stiffness
    !===============================================================================
    
    ! This subrotine computes the stiffness matrix for the viscoelastic elements.
    ! WARNING: it has to be modified for gravity perturbation.
    subroutine compute_stiffness_viscoelastic(dt,errcode,errtag)
    use set_precision
    use global,only:myrank,NDIM,nst,ngll,nedof,nedofu,nedofphi,nenode,ngnode,      &
    ngllx,nglly,ngllz,ngll,g_coord,gdof_elmt,g_num,mat_domain,mat_id,massdens_elmt,&
    bulkmod_elmt,shearmod_elmt,rho_blk,ym_blk,imat_to_imatve, &
    ISDISP_DOF,ISPOT_DOF,storederiv,        &
    POT_TYPE,PGRAVITY,PMAGNETIC,    &
    storejw,devel_nondim,isdxval,isdyval,isdzval, &
    edofu,edofphi,grav0_nodal,dgrav0_elmt,ISGRAV0, &
    muratio_blk,visco_model,VISCO_MAXWELL,VISCO_ZENER,VISCO_GENMAXWELL,nmaxwell
    use global,only:nelmt_viscoelas,eid_viscoelas,relaxtime,storekmat
    use element,only:hex8_gnode
    use viscoelastic,only:compute_cmat_maxwell,compute_cmat_zener, &
    compute_cmat_genmaxwell
    use math_constants,only:HALF,ONE,ZERO,FOUR,GRAV_CONS,PI
    use math_library,only:determinant,invert,issymmetric
    use weakform
    use shape_library
    use gll_library
    use infinite_element
    implicit none
    real(kind=kreal),intent(in) :: dt
    integer,intent(out) :: errcode
    character(len=250),intent(out) :: errtag
    integer :: i
    integer :: i_elmt,ielmt,imat
    integer :: imatve,mdomain
    integer :: num(nenode)
    real(kind=kreal) :: cmat(nst,nst)
    real(kind=kreal) :: jacw !jacobian*weight
    real(kind=kreal) :: eload(nedof)
    real(kind=kreal) :: dgmat(NDIM,NDIM),eg0(ngll,NDIM),g0(NDIM),dg0(6)
    
    real(kind=kreal) :: kmatu(nedofu,nedofu)
    
    real(kind=kreal) :: interpf(NGLL),deriv(NDIM,nenode)
    real(kind=kreal) :: Mvec(NDIM)
    real(kind=kreal) :: bmatu(NST,NEDOFU),wmat_gradphi(NEDOFU,NDIM),      &
    rmat_gradphi(NDIM,NEDOFPHI),wmat_sPE(NEDOFPHI,NDIM),rmat_sPE(NDIM,NEDOFU), &
    rmat_term1(NDIM,NEDOFU),wmat_term1(NEDOFU,NDIM),                           &
    rmat_term2(NDIM,NEDOFU),wmat_term2(NEDOFU,NDIM),                           &
    rmat_term3(1,NEDOFU),wmat_term3(NEDOFU,1),                                 &
    rmat_term4(1,NEDOFU),wmat_term4(NEDOFU,1)
    
    real(kind=kreal) :: muratio(nmaxwell),tratio(nmaxwell)
    
    errtag="ERROR: unknown!"
    errcode=-1
    errsrc=trim(myfname)//' => compute_stiffness_viscoelastic'
    
    ! Viscoelastic part of storekmat will be replaced. 
    ! kmat(edofu,edofu) is the only viscoelastic portion.
    ! kmat(edofphi,edofphi) and other offdiagonal terms must remain the same.
    
    ! If statement here isn't necessary since viscoelasticity matters only if 
    ! if there are displacement DOFs.
    if(.not.ISDISP_DOF)return
    
    ! Viscoelastic elements
    do i_elmt=1,nelmt_viscoelas
      ielmt=eid_viscoelas(i_elmt)
      imat=mat_id(ielmt)
      imatve=imat_to_imatve(imat)
      muratio=muratio_blk(:,imatve)
      ! For a more gneneral case tratio can also be variable within an element
      tratio=dt/relaxtime(:,imatve)
      
      kmatu=zero
      do i=1,ngll
        deriv=storederiv(:,:,i,ielmt)
        jacw=storejw(i,ielmt)
       
        call compute_cmat_genmaxwell(bulkmod_elmt(i,ielmt),shearmod_elmt(i,ielmt),&
        tratio,muratio,cmat)
    
        call compute_bmat_stress(deriv,bmatu)
        kmatu=kmatu+matmul(matmul(transpose(bmatu),cmat),bmatu)*jacw
      enddo !i
      storekmat(edofu,edofu,ielmt)=kmatu
    
    enddo ! i_elmt
    
    end subroutine compute_stiffness_viscoelastic
    !===============================================================================
    
    ! This subrotine computes the stiffness matrix, diagonal preconditioner, and
    ! body loads (gravity and pseudostatic loads) optionally
    ! TODO: optional precoditioner,optional assembly of stiffness
    subroutine stiffness_bodyload(nelmt,neq,gnod,g_num,gdof_elmt,mat_id,gam,       &
    storkm,dprecon,extload,gravity,pseudoeq)
    use set_precision
    use global,only:NDIM,nst,ngll,nedof,nenode,ngnode,g_coord,eqkx,eqky,eqkz,      &
    nmatblk,bulkmod_elmt,shearmod_elmt
    use elastic,only:compute_cmat_elastic
    use math_library,only:determinant,invert
    use weakform,only:compute_bmat_stress
    use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights
    implicit none
    integer,intent(in) :: nelmt,neq,gnod(8) ! nelmt (only intact elements)
    integer,intent(in) :: g_num(nenode,nelmt),gdof_elmt(nedof,nelmt),mat_id(nelmt)
    ! only intact elements
    real(kind=kreal),intent(in) :: gam(nmatblk)
    real(kind=kreal),intent(out) :: storkm(nedof,nedof,nelmt),dprecon(0:neq)
    real(kind=kreal),intent(inout),optional :: extload(0:neq)
    logical,intent(in),optional :: gravity,pseudoeq
    
    real(kind=kreal) :: detjac,zero=0.0_kreal
    real(kind=kreal) :: cmat(nst,nst),coord(ngnode,NDIM),jac(NDIM,NDIM),           &
    deriv(NDIM,nenode),bmatu(nst,nedof),eld(nedof),eqload(nedof),km(nedof,nedof)
    integer :: egdof(nedof),num(nenode)
    integer :: i,i_elmt,k
    
    if(present(extload).and.(.not.present(gravity) .or. .not.present(pseudoeq)))then
      write(*,'(a)')'ERROR: both "gravity" and "pseudoeq" must be defined for &
      &"extload"!'
      stop
    endif
    
    ! compute stiffness matrices
    storkm=zero; dprecon=zero
    do i_elmt=1,nelmt
      num=g_num(:,i_elmt)
      coord=transpose(g_coord(:,num(gnod)))
      egdof=gdof_elmt(:,i_elmt)
      km=zero; eld=zero; eqload=zero
      do i=1,ngll
        call compute_cmat_elastic(bulkmod_elmt(i,i_elmt),shearmod_elmt(i,i_elmt),  &
        cmat)
        ! compute Jacobian at GLL point using 20 noded element
        jac=matmul(dshape_hex8(:,:,i),coord)
        detjac=determinant(jac)
        call invert(jac)
    
        deriv=matmul(jac,dlagrange_gll(:,i,:))
        call compute_bmat_stress(deriv,bmatu)
        km=km+matmul(matmul(transpose(bmatu),cmat),bmatu)*detjac*gll_weights(i)
        eld(3:nedof:3)=eld(3:nedof:3)+lagrange_gll(i,:)*detjac*gll_weights(i)
        !eld(2:nedof-1:3)=eld(2:nedof-1:3)+fun(:)*detjac*weights(i)
      enddo ! i=1,ngll
      storkm(:,:,i_elmt)=km
      do k=1,nedof
        dprecon(egdof(k))=dprecon(egdof(k))+km(k,k)
      enddo
    
      if(.not.present(extload))cycle
      ! compute body loads
      ! gravity load and add to extload
      ! WARNING: this must be changed to massdens_elmt
      if(gravity)extload(egdof)=extload(egdof)-eld*gam(mat_id(i_elmt))
      if(pseudoeq)then
        ! compute pseudostatic earthquake loads and add to extload
        eqload(1:nedof:3)=eqkx*eld(3:nedof:3)
        eqload(2:nedof:3)=eqky*eld(3:nedof:3)
        eqload(3:nedof:3)=eqkz*eld(3:nedof:3)
        extload(egdof)=extload(egdof)+eqload*gam(mat_id(i_elmt)) ! KN
      endif
    enddo ! i_elmt=1,nelmt
    !write(*,*)'complete!'
    dprecon(0)=zero
    if(present(extload))extload(0)=zero
    end subroutine stiffness_bodyload
    !===============================================================================
    
    ! This subrotine computes the mass matrix.
    subroutine compute_mass_elastic(storemmat,errcode,errtag)
    use set_precision
    use global,only:myrank,NDIM,nst,nelmt,ngll,nenode,ngnode,&
    ngllx,nglly,ngllz,ngll,g_coord,gdof_elmt,g_num,massdens_elmt,&
    storejw,devel_nondim, &
    isdxval,isdyval,isdzval,devel_gaminf,infquad, &
    edofu,edofphi,grav0_nodal,dgrav0_elmt,ISGRAV0,&
    imat_to_imatmag,magnetization_blk,ismagnet_blk
    use element,only:hex8_gnode,map2exodus_hex8
    use math_constants,only:HALF,ONE,ZERO,FOUR,GRAV_CONS,PI
    use math_library,only:determinant,invert,issymmetric
    use weakform
    use shape_library
    use gll_library
    use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights,       &
    prepare_integration
    !use ieee_arithmetic
    implicit none
    real(kind=kreal),intent(out) :: storemmat(:,:)
    integer,intent(out) :: errcode
    character(len=250),intent(out) :: errtag
    integer :: i,i_gll
    integer :: i_elmt,ielmt,imat,ignode
    integer :: num(nenode)
    real(kind=kreal) :: detjac !determinant of Jacobian
    real(kind=kreal) :: coord(ngnode,NDIM),jac(NDIM,NDIM)
    
    real(kind=kreal) :: xval
    real(kind=kreal) :: interpf(NGLL),deriv(NDIM,nenode)
    
    ! jacw=jacobian*weight
    real(kind=kreal) :: jacw
    integer :: nip
    
    errtag="ERROR: unknown!"
    errcode=-1
    errsrc=trim(myfname)//' => compute_mass_elastic'
    
    storemmat=zero
    ! Elastic elements
    ! Following loops through nelmt
    do i_elmt=1,nelmt
      ielmt=i_elmt
      nip=ngll
      
      do i=1,nip
          ignode=num(i)
          jacw=storejw(i,ielmt)
          ! The mass matrix element for X, Y, and Z degrees of freedom at 
          ! a GLL point is same.
          ! Therefore, we store only the one value per GLL point.
    
          ! NOTE: there is no quadrature summation in the following statement
          ! because the summation results in the Kronecker's delta giving 
          ! a simple GLL point-wise expression for the mass matrix element. 
          storemmat(i,i_elmt)=storemmat(i,i_elmt)+massdens_elmt(i,ielmt)*jacw
      enddo
    enddo ! i_elmt
    
    end subroutine compute_mass_elastic
    !===============================================================================
    
    ! This subrotine computes the mass matrix.
    subroutine compute_mass_elastic_global(storemmat_global,errcode,errtag)
    use set_precision
    use global,only:myrank,NDIM,nst,nelmt,ngll,nenode,ngnode,&
    ngllx,nglly,ngllz,ngll,g_coord,gdof_elmt,g_num,massdens_elmt,&
    storejw,devel_nondim, &
    isdxval,isdyval,isdzval,devel_gaminf,infquad, &
    edofu,edofphi,grav0_nodal,dgrav0_elmt,ISGRAV0,&
    imat_to_imatmag,magnetization_blk,ismagnet_blk
    use element,only:hex8_gnode,map2exodus_hex8
    use math_constants,only:HALF,ONE,ZERO,FOUR,GRAV_CONS,PI
    use math_library,only:determinant,invert,issymmetric
    use weakform
    use shape_library
    use gll_library
    use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights,       &
    prepare_integration
    !use ieee_arithmetic
    implicit none
    real(kind=kreal),intent(out) :: storemmat_global(:)
    integer,intent(out) :: errcode
    character(len=250),intent(out) :: errtag
    integer :: i,i_gll
    integer :: i_elmt,ielmt,imat,ignode
    integer :: num(nenode)
    real(kind=kreal) :: detjac !determinant of Jacobian
    real(kind=kreal) :: coord(ngnode,NDIM),jac(NDIM,NDIM)
    
    real(kind=kreal) :: xval
    real(kind=kreal) :: interpf(NGLL),deriv(NDIM,nenode)
    
    ! jacw=jacobian*weight
    real(kind=kreal) :: jacw
    integer :: nip
    
    errtag="ERROR: unknown!"
    errcode=-1
    errsrc=trim(myfname)//' => compute_mass_elastic'
    
    storemmat_global=zero
    ! Elastic elements
    ! Following loops through nelmt
    do i_elmt=1,nelmt
      ielmt=i_elmt
      num=g_num(:,ielmt)
      coord=transpose(g_coord(:,num(hex8_gnode)))
      nip=ngll
    
      do i=1,nip
          ignode=num(i)
          jacw=storejw(i,ielmt)
        ! NOTE: there is no quadrature summation in the following statement
        ! because the summation results in the Kronecker's delta giving 
        ! a simple GLL point-wise expression for the mass matrix element. 
         
          storemmat_global(ignode)=storemmat_global(ignode)+massdens_elmt(i,ielmt)*jacw
      enddo
    enddo ! i_elmt
    
    end subroutine compute_mass_elastic_global
    !===============================================================================
    
    ! This subrotine computes the stiffness matrix, and
    ! body loads contributed by mass density or magnetiztion.
    ! TODO: optional precoditioner,optional assembly of stiffness
    subroutine compute_stiffness_elasticOLD(storekmat,rhoload,errcode,errtag)
    use set_precision
    use global,only:myrank,NDIM,nst,nelmt,ngll,nedof,nedofu,nedofphi,nenode,ngnode,&
    ngllx,nglly,ngllz,ngll,g_coord,gdof_elmt,g_num,mat_domain,mat_id,massdens_elmt,&
    bulkmod_elmt,shearmod_elmt,isempty_blk,rho_blk,ym_blk,magnetization_elmt,&
    infinite_iface,infinite_face_idir,pole_coord0,pole_coord1,       &
    pole_type,pole_axis,axis_range,ISDISP_DOF,ISPOT_DOF,storederiv,    &
    POT_TYPE,PGRAVITY,PMAGNETIC,    &
    storejw,devel_nondim, &
    isdxval,isdyval,isdzval,devel_gaminf,infquad, &
    edofu,edofphi,grav0_nodal,dgrav0_elmt,ISGRAV0,&
    imat_to_imatmag,magnetization_blk,ismagnet_blk
    use element,only:hex8_gnode,map2exodus_hex8
    use elastic,only:compute_cmat_elastic
    use math_constants,only:HALF,ONE,ZERO,FOUR,GRAV_CONS,PI
    use math_library,only:determinant,invert,issymmetric
    use weakform
    use shape_library
    use gll_library
    use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights,       &
    prepare_integration
    use infinite_element
    !use ieee_arithmetic
    implicit none
    real(kind=kreal),intent(out) :: storekmat(:,:,:)
    real(kind=kreal),intent(out) :: rhoload(0:)
    integer,intent(out) :: errcode
    character(len=250),intent(out) :: errtag
    real(kind=kreal),parameter :: FOUR_PI_G=FOUR*PI*GRAV_CONS
    real(kind=kreal),parameter :: FOUR_PI_G_INV=ONE/FOUR_PI_G
    integer :: i,i_gll
    integer :: i_elmt,ielmt,imat
    integer :: imatve,mdomain
    integer :: num(nenode),egdof(nedof)
    real(kind=kreal) :: cmat(nst,nst)
    real(kind=kreal) :: detjac !determinant of Jacobian
    real(kind=kreal) :: coord(ngnode,NDIM),jac(NDIM,NDIM),eload(nedof)
    real(kind=kreal) :: dgmat(NDIM,NDIM),eg0(ngll,NDIM),g0(NDIM),dg0(6)
    
    real(kind=kreal) :: kmat(nedof,nedof),xval
    real(kind=kreal) :: interpf(NGLL),deriv(NDIM,nenode)
    
    ! magnetization
    integer :: imatmag
    real(kind=kreal) :: M(NDIM),Mgll(NDIM,ngll)
    real(kind=kreal) :: divM,dxMx,dyMy,dzMz
    
    real(kind=kreal) :: bmatu(NST,NEDOFU),wmat_gradphi(NEDOFU,NDIM),      &
    rmat_gradphi(NDIM,NEDOFPHI),wmat_sPE(NEDOFPHI,NDIM),rmat_sPE(NDIM,NEDOFU), &
    rmat_term1(NDIM,NEDOFU),wmat_term1(NEDOFU,NDIM),                           &
    rmat_term2(NDIM,NEDOFU),wmat_term2(NEDOFU,NDIM),                           &
    rmat_term3(1,NEDOFU),wmat_term3(NEDOFU,1),                                 &
    rmat_term4(1,NEDOFU),wmat_term4(NEDOFU,1)
    
    real(kind=kreal) :: tratio
    
    ! for infinite elements
    integer,parameter :: nginf=8
    integer :: i_face,idir
    real(kind=kreal) :: tcoord
    real(kind=kreal) :: gaminf
    real(kind=kreal) :: polex(4,NDIM)
    real(kind=kreal) :: coordinf(nginf,NDIM)
    ! jacw=jacobian*weight
    real(kind=kreal) :: jacw
    real(kind=kreal),allocatable :: shape_infinite(:,:),dshape_infinite(:,:,:)
    real(kind=kreal),allocatable :: lagrange_gl(:,:),dlagrange_gl(:,:,:)
    real(kind=kreal),allocatable :: GLw(:)
    integer :: nip,nipinf
    logical :: isinf ! flag to check if the element if on the infinite domain
    logical :: isfaces(6)
    
    errtag="ERROR: unknown!"
    errcode=-1
    errsrc=trim(myfname)//' => compute_stiffness_elastic'
    
    ! use ng for Gauss quadrature and ngll fro GLL-Radau quadrature
    gaminf=2.00_kreal !1.99_kreal ! 2.0: X1 in the mid-position. gaminf must be > 1.0
    if(devel_gaminf.gt.ONE .and. devel_gaminf.lt.FOUR)gaminf=devel_gaminf
    nipinf=ngll
    
    allocate(shape_infinite(nipinf,nginf),dshape_infinite(NDIM,nipinf,nginf))
    allocate(lagrange_gl(nipinf,ngll),dlagrange_gl(NDIM,nipinf,ngll))
    
    allocate(GLw(nipinf))
    
    storekmat=zero
    rhoload=zero
    ! Purely elastic elements
    ! Viscoelastic elements are elastic at time = 0
    ! Following loops through nelmt_elas+nelmt_viscoelas
    do i_elmt=1,nelmt
      ielmt=i_elmt
      num=g_num(:,ielmt)
      imat=mat_id(ielmt)
      mdomain=mat_domain(imat) 
      coord=transpose(g_coord(:,num(hex8_gnode)))
      nip=ngll
    
      ! set magnetization at GLL points
      Mgll=ZERO
      if(POT_TYPE==PMAGNETIC)then
        if(ismagnet_blk(imat))then
          ! ONLY FOR BLOCK PROPERTIES
          !imatmag=imat_to_imatmag(imat)
          !M=magnetization_blk(:,imatmag)
          !do i_gll=1,ngll 
          !  Mgll(:,i_gll)=M
          !enddo
          do i_gll=1,ngll 
            Mgll(:,i_gll)=magnetization_elmt(:,i_gll,ielmt)
          enddo
        endif
      endif
      !H=dgrav0_elmt(:,:,ielmt)
      eg0=transpose(grav0_nodal(:,num))
    
      isfaces=infinite_iface(:,ielmt)
      isinf=any(isfaces)
      if(count(isfaces).gt.1)isinf=.false.
      ! set coordinates
      if(isinf)then
        !coordinf=transpose(g_coord(:,num(gnodinf)))
        coordinf=transpose(g_coord(:,num(hex8_gnode)))
    
        ! Loop through infinite faces
        do i_face=1,6
          if(.not.isfaces(i_face))cycle
          idir=infinite_face_idir(i_face,ielmt)
          if(i_face==1)then
            ! ymin face
            ! Set X coordinate of the pole
            if(trim(pole_type)=='plane')then
              polex(1,:)=coordinf(3,:)!coordinf(1,:)
              polex(2,:)=coordinf(4,:)!coordinf(4,:)
              polex(3,:)=coordinf(7,:)!coordinf(5,:)
              polex(4,:)=coordinf(8,:)!coordinf(8,:)
              polex(:,idir)=pole_coord0(idir)
            elseif(trim(pole_type)=='axis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              polex(1,pole_axis)=coordinf(3,pole_axis)
              polex(2,pole_axis)=coordinf(4,pole_axis)
              polex(3,pole_axis)=coordinf(7,pole_axis)
              polex(4,pole_axis)=coordinf(8,pole_axis)
            elseif(trim(pole_type)=='pointaxis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              tcoord=coordinf(3,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(1,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(1,:)=pole_coord1
              endif
              tcoord=coordinf(4,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(2,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(2,:)=pole_coord1
              endif
              tcoord=coordinf(7,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(3,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(3,:)=pole_coord1
              endif
              tcoord=coordinf(8,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(4,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(4,:)=pole_coord1
              endif
            else
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
            endif
            ! X2
            coordinf(2,:)=polex(1,:)+gaminf*(coordinf(3,:)-polex(1,:))
            coordinf(1,:)=polex(2,:)+gaminf*(coordinf(4,:)-polex(2,:))
            coordinf(6,:)=polex(3,:)+gaminf*(coordinf(7,:)-polex(3,:))
            coordinf(5,:)=polex(4,:)+gaminf*(coordinf(8,:)-polex(4,:))
            ! Set near face coordinates to the pole coordinates
            ! X0
            coordinf(3,:)=polex(1,:)
            coordinf(4,:)=polex(2,:)
            coordinf(7,:)=polex(3,:)
            coordinf(8,:)=polex(4,:)
          
          elseif(i_face==2)then
            ! xmax face
            ! Set X coordinate of the pole
            if(trim(pole_type)=='plane')then
              polex(1,:)=coordinf(1,:)!coordinf(1,:)
              polex(2,:)=coordinf(4,:)!coordinf(4,:)
              polex(3,:)=coordinf(5,:)!coordinf(5,:)
              polex(4,:)=coordinf(8,:)!coordinf(8,:)
              polex(:,idir)=pole_coord0(idir)
            elseif(trim(pole_type)=='axis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              polex(1,pole_axis)=coordinf(1,pole_axis)
              polex(2,pole_axis)=coordinf(4,pole_axis)
              polex(3,pole_axis)=coordinf(5,pole_axis)
              polex(4,pole_axis)=coordinf(8,pole_axis)
            elseif(trim(pole_type)=='pointaxis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              tcoord=coordinf(1,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(1,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(1,:)=pole_coord1
              endif
              tcoord=coordinf(4,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(2,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(2,:)=pole_coord1
              endif
              tcoord=coordinf(5,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(3,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(3,:)=pole_coord1
              endif
              tcoord=coordinf(8,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(4,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(4,:)=pole_coord1
              endif
            else
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
            endif
            ! X2
            coordinf(2,:)=polex(1,:)+gaminf*(coordinf(1,:)-polex(1,:))
            coordinf(3,:)=polex(2,:)+gaminf*(coordinf(4,:)-polex(2,:))
            coordinf(6,:)=polex(3,:)+gaminf*(coordinf(5,:)-polex(3,:))
            coordinf(7,:)=polex(4,:)+gaminf*(coordinf(8,:)-polex(4,:))
            ! Set near face coordinates to the pole coordinates
            ! X0
            coordinf(1,:)=polex(1,:)
            coordinf(4,:)=polex(2,:)
            coordinf(5,:)=polex(3,:)
            coordinf(8,:)=polex(4,:)
          
          elseif(i_face==3)then
            ! ymax face
            ! Set Y coordinate of the pole 
            if(trim(pole_type)=='plane')then
              polex(1,:)=coordinf(1,:)
              polex(2,:)=coordinf(2,:)
              polex(3,:)=coordinf(5,:)
              polex(4,:)=coordinf(6,:)
              polex(:,idir)=pole_coord0(idir)
            elseif(trim(pole_type)=='axis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              polex(1,pole_axis)=coordinf(1,pole_axis)
              polex(2,pole_axis)=coordinf(2,pole_axis)
              polex(3,pole_axis)=coordinf(5,pole_axis)
              polex(4,pole_axis)=coordinf(6,pole_axis)
            elseif(trim(pole_type)=='pointaxis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              tcoord=coordinf(1,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(1,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(1,:)=pole_coord1
              endif
              tcoord=coordinf(2,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(2,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(2,:)=pole_coord1
              endif
              tcoord=coordinf(5,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(3,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(3,:)=pole_coord1
              endif
              tcoord=coordinf(6,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(4,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(4,:)=pole_coord1
              endif
            else
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
            endif
            ! X2
            coordinf(4,:)=polex(1,:)+gaminf*(coordinf(1,:)-polex(1,:))
            coordinf(3,:)=polex(2,:)+gaminf*(coordinf(2,:)-polex(2,:))
            coordinf(8,:)=polex(3,:)+gaminf*(coordinf(5,:)-polex(3,:))
            coordinf(7,:)=polex(4,:)+gaminf*(coordinf(6,:)-polex(4,:))
            ! Set near face coordinates to the pole coordinates
            ! X0
            coordinf(1,:)=polex(1,:)
            coordinf(2,:)=polex(2,:)
            coordinf(5,:)=polex(3,:)
            coordinf(6,:)=polex(4,:)
    
          elseif(i_face==4)then
            ! xmin face
            ! Set Y coordinate of the pole 
            if(trim(pole_type)=='plane')then
              polex(1,:)=coordinf(2,:)
              polex(2,:)=coordinf(3,:)
              polex(3,:)=coordinf(6,:)
              polex(4,:)=coordinf(7,:)
              polex(:,idir)=pole_coord0(idir)
            elseif(trim(pole_type)=='axis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              polex(1,pole_axis)=coordinf(2,pole_axis)
              polex(2,pole_axis)=coordinf(3,pole_axis)
              polex(3,pole_axis)=coordinf(6,pole_axis)
              polex(4,pole_axis)=coordinf(7,pole_axis)
            elseif(trim(pole_type)=='pointaxis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              tcoord=coordinf(2,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(1,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(1,:)=pole_coord1
              endif
              tcoord=coordinf(3,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(2,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(2,:)=pole_coord1
              endif
              tcoord=coordinf(6,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(3,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(3,:)=pole_coord1
              endif
              tcoord=coordinf(7,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(4,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(4,:)=pole_coord1
              endif
            else
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
            endif
            ! X2
            coordinf(1,:)=polex(1,:)+gaminf*(coordinf(2,:)-polex(1,:))
            coordinf(4,:)=polex(2,:)+gaminf*(coordinf(3,:)-polex(2,:))
            coordinf(5,:)=polex(3,:)+gaminf*(coordinf(6,:)-polex(3,:))
            coordinf(8,:)=polex(4,:)+gaminf*(coordinf(7,:)-polex(4,:))
            ! Set near face coordinates to the pole coordinates
            ! X0
            coordinf(2,:)=polex(1,:)
            coordinf(3,:)=polex(2,:)
            coordinf(6,:)=polex(3,:)
            coordinf(7,:)=polex(4,:)
    
          elseif(i_face==5)then
            ! zmin face
            ! Set Y coordinate of the pole 
            if(trim(pole_type)=='plane')then
              polex(1,:)=coordinf(5,:)
              polex(2,:)=coordinf(6,:)
              polex(3,:)=coordinf(8,:)
              polex(4,:)=coordinf(7,:)
              polex(:,idir)=pole_coord0(idir)
            elseif(trim(pole_type)=='axis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              polex(1,pole_axis)=coordinf(5,pole_axis)
              polex(2,pole_axis)=coordinf(6,pole_axis)
              polex(3,pole_axis)=coordinf(8,pole_axis)
              polex(4,pole_axis)=coordinf(7,pole_axis)
            elseif(trim(pole_type)=='pointaxis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              tcoord=coordinf(5,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(1,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(1,:)=pole_coord1
              endif
              tcoord=coordinf(6,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(2,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(2,:)=pole_coord1
              endif
              tcoord=coordinf(8,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(3,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(3,:)=pole_coord1
              endif
              tcoord=coordinf(7,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(4,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(4,:)=pole_coord1
              endif
            else
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
            endif
            ! X2
            coordinf(1,:)=polex(1,:)+gaminf*(coordinf(5,:)-polex(1,:))
            coordinf(2,:)=polex(2,:)+gaminf*(coordinf(6,:)-polex(2,:))
            coordinf(4,:)=polex(3,:)+gaminf*(coordinf(8,:)-polex(3,:))
            coordinf(3,:)=polex(4,:)+gaminf*(coordinf(7,:)-polex(4,:))
            ! Set near face coordinates to the pole coordinates
            ! X0
            coordinf(5,:)=polex(1,:)
            coordinf(6,:)=polex(2,:)
            coordinf(8,:)=polex(3,:)
            coordinf(7,:)=polex(4,:)
    
          elseif(i_face==6)then
            ! zmax face
            ! Set Z coordinate of the pole 
            if(trim(pole_type)=='plane')then
              polex(1,:)=coordinf(1,:)
              polex(2,:)=coordinf(2,:)
              polex(3,:)=coordinf(3,:)
              polex(4,:)=coordinf(4,:)
              polex(:,idir)=pole_coord0(idir)
            elseif(trim(pole_type)=='axis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              polex(1,pole_axis)=coordinf(1,pole_axis)
              polex(2,pole_axis)=coordinf(2,pole_axis)
              polex(3,pole_axis)=coordinf(3,pole_axis)
              polex(4,pole_axis)=coordinf(4,pole_axis)
            elseif(trim(pole_type)=='pointaxis')then
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
              tcoord=coordinf(1,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(1,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(1,:)=pole_coord1
              endif
              tcoord=coordinf(2,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(2,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(2,:)=pole_coord1
              endif
              tcoord=coordinf(3,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(3,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(3,:)=pole_coord1
              endif
              tcoord=coordinf(4,pole_axis)
              if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
                polex(4,pole_axis)=tcoord
              endif
              if(tcoord.ge.axis_range(2))then
                polex(4,:)=pole_coord1
              endif
            else
              polex(1,:)=pole_coord0
              polex(2,:)=pole_coord0
              polex(3,:)=pole_coord0
              polex(4,:)=pole_coord0
            endif
            ! X2
            coordinf(5,:)=polex(1,:)+gaminf*(coordinf(1,:)-polex(1,:))
            coordinf(6,:)=polex(2,:)+gaminf*(coordinf(2,:)-polex(2,:))
            coordinf(7,:)=polex(3,:)+gaminf*(coordinf(3,:)-polex(3,:))
            coordinf(8,:)=polex(4,:)+gaminf*(coordinf(4,:)-polex(4,:))
            ! Set near face coordinates to the pole coordinates
            ! X0
            coordinf(1,:)=polex(1,:)
            coordinf(2,:)=polex(2,:)
            coordinf(3,:)=polex(3,:)
            coordinf(4,:)=polex(4,:)
    
          endif
    
        enddo ! i_face
    
        ! convert gnode ordering to indicial order to match the infinite shape
        ! functions ordering
        ! if we use map2exodus_hex8 on exodus order it becomes indicial order
        coordinf=coordinf(map2exodus_hex8,:)
        nip=nipinf !ngll ! we use GLL or GLL-Radau quadrature
        
        ! Radau or Gauss quadrature
        call shape_function_infiniteGLHEX8ZW(infquad,ngllx,nglly,ngllz,    &
        ngll,nip,isfaces,shape_infinite,dshape_infinite,lagrange_gl,       &
        dlagrange_gl,GLw)
        
      endif
    
      egdof=gdof_elmt(:,i_elmt)
        
      kmat=zero
      eload=zero
      do i=1,nip
        if(ISDISP_DOF)then
          call compute_cmat_elastic(bulkmod_elmt(i,ielmt),shearmod_elmt(i,ielmt),  &
          cmat)
        endif
    
        if(isinf)then
          ! infinite element
          interpf=lagrange_gl(i,:)
          
          jac=matmul(dshape_infinite(:,i,:),coordinf)
          detjac=determinant(jac)
          if(detjac.le.zero.and.myrank==0)then
            write(*,*)'ERROR: zero or negative jacobian in infinite element!'
            write(*,*)'HINT: check "pole_type" and "infquad"!'
            write(*,*)'HINT: make sure that coordinates units are consistent!'
            write(*,*)myrank,i_elmt,i,nip,detjac
            write(*,*)isfaces
            write(*,*)infquad
            write(*,*)coordinf
            stop
          endif
          call invert(jac)
          deriv=matmul(jac,dlagrange_gl(:,i,:))
    
          ! set derivative constraint
          if(isdxval)deriv(1,:)=ZERO
          if(isdyval)deriv(2,:)=ZERO
          if(isdzval)deriv(3,:)=ZERO
          
          deriv=storederiv(:,:,i,ielmt) 
          !jacw=detjac*GLw(i)
          jacw=storejw(i,ielmt)
          
          if(ISDISP_DOF)then
            ! compute only for nonempty elements
            if(.not.isempty_blk(imat))then
              call compute_bmat_stress(deriv,bmatu)
              kmat(edofu,edofu)=kmat(edofu,edofu)+matmul(matmul(transpose(bmatu),cmat),bmatu)*jacw
    
              if(ISGRAV0)then
                g0=eg0(i,:)
                ! compute dg0=\nabla g0
                dgmat=matmul(deriv,eg0)
                dg0(1)=dgmat(1,1)
                dg0(2)=dgmat(2,2)
                dg0(3)=dgmat(3,3)
                dg0(4)=dgmat(1,2)
                dg0(5)=dgmat(1,3)
                dg0(6)=dgmat(2,3)
    
                call compute_rmat_term1(massdens_elmt(i,ielmt),interpf,rmat_term1)
                call compute_wmat_term1(g0,dg0,interpf,deriv,wmat_term1)
                call compute_rmat_term2(g0,dg0,interpf,deriv,rmat_term2)
                call compute_wmat_term2(massdens_elmt(i,ielmt),interpf,wmat_term2)
                kmat(edofu,edofu)=kmat(edofu,edofu)-HALF*(matmul(wmat_term1,rmat_term1)+ &
                matmul(wmat_term2,rmat_term2))*jacw
      
                call compute_rmat_term3(massdens_elmt(i,ielmt),g0,interpf,rmat_term3)
                call compute_wmat_term3(deriv,wmat_term3)
                call compute_rmat_term4(massdens_elmt(i,ielmt),deriv,rmat_term4)
                call compute_wmat_term4(g0,interpf,wmat_term4)
                kmat(edofu,edofu)=kmat(edofu,edofu)+HALF*(matmul(wmat_term3,rmat_term3)+ &
                matmul(wmat_term4,rmat_term4))*jacw
              endif
              
              if(ISPOT_DOF)then
                ! w.rho*grad(phi)
                call compute_wmat_gradphi(interpf,wmat_gradphi)
                call compute_rmat_gradphi(massdens_elmt(i,ielmt),deriv,rmat_gradphi)          
                kmat(edofu,edofphi)=kmat(edofu,edofphi)+matmul(wmat_gradphi,rmat_gradphi)*jacw
                ! grad(w).rho*s
                ! this term is transpose of the previous. Therefore it will be
                ! computed later
                call compute_wmat_sPE(deriv,wmat_sPE)
                call compute_rmat_sPE(massdens_elmt(i,ielmt),interpf,rmat_sPE)
                kmat(edofphi,edofu)=kmat(edofphi,edofu)+matmul(wmat_sPE,rmat_sPE)*jacw
              endif
            endif
          endif
    
          if(ISPOT_DOF)then
            kmat(edofphi,edofphi)=kmat(edofphi,edofphi)+matmul(transpose(deriv),deriv)*jacw
            if(.not.ISDISP_DOF)then
              if(POT_TYPE==PMAGNETIC)then
                if(ismagnet_blk(imat))then
                  ! first compute divegence of M
                  ! \nabla.M=dxMx+dyMy+dzMz
                  dxMx=dot_product(deriv(1,:),Mgll(1,:))
                  dyMy=dot_product(deriv(2,:),Mgll(2,:))
                  dzMz=dot_product(deriv(3,:),Mgll(3,:))
                  divM=dxMx+dyMy+dzMz
                  if(maxval(abs(M)).gt.ZERO)then
                    eload(edofphi)=eload(edofphi)+lagrange_gl(i,:)*divM*jacw
                  endif
                endif
              else
                eload(edofphi)=eload(edofphi)+lagrange_gl(i,:)*massdens_elmt(i,ielmt)*jacw
              endif
            endif
          endif
    
        else ! (isinf) 
          ! standard element
          interpf=lagrange_gll(i,:)
        
          jac=matmul(dshape_hex8(:,:,i),coord)
          detjac=determinant(jac)
          call invert(jac)
          deriv=matmul(jac,dlagrange_gll(:,i,:))
          ! set derivative constraint
          if(isdxval)deriv(1,:)=ZERO
          if(isdyval)deriv(2,:)=ZERO
          if(isdzval)deriv(3,:)=ZERO
          
          deriv=storederiv(:,:,i,ielmt) 
          !jacw=detjac*gll_weights(i)
          jacw=storejw(i,ielmt)
         
          if(ISDISP_DOF)then
            ! compute only for nonempty elements
            if(.not.isempty_blk(imat))then
              call compute_bmat_stress(deriv,bmatu)
              kmat(edofu,edofu)=kmat(edofu,edofu)+matmul(matmul(transpose(bmatu),cmat),bmatu)*jacw
               
              if(ISGRAV0)then
                g0=eg0(i,:)
                ! compute dg0=\nabla g0
                dgmat=matmul(deriv,eg0)
                dg0(1)=dgmat(1,1)
                dg0(2)=dgmat(2,2)
                dg0(3)=dgmat(3,3)
                dg0(4)=dgmat(1,2)
                dg0(5)=dgmat(1,3)
                dg0(6)=dgmat(2,3)
    
                call compute_rmat_term1(massdens_elmt(i,ielmt),interpf,rmat_term1)
                call compute_wmat_term1(g0,dg0,interpf,deriv,wmat_term1)
                call compute_rmat_term2(g0,dg0,interpf,deriv,rmat_term2)
                call compute_wmat_term2(massdens_elmt(i,ielmt),interpf,wmat_term2)
                kmat(edofu,edofu)=kmat(edofu,edofu)-HALF*(matmul(wmat_term1,rmat_term1)+ &
                matmul(wmat_term2,rmat_term2))*jacw
      
                call compute_rmat_term3(massdens_elmt(i,ielmt),g0,interpf,rmat_term3)
                call compute_wmat_term3(deriv,wmat_term3)
                call compute_rmat_term4(massdens_elmt(i,ielmt),deriv,rmat_term4)
                call compute_wmat_term4(g0,interpf,wmat_term4)
                kmat(edofu,edofu)=kmat(edofu,edofu)+HALF*(matmul(wmat_term3,rmat_term3)+ &
                matmul(wmat_term4,rmat_term4))*jacw
              endif
              if(ISPOT_DOF)then
                ! w.rho*grad(phi)
                call compute_wmat_gradphi(interpf,wmat_gradphi)
                call compute_rmat_gradphi(massdens_elmt(i,ielmt),deriv,rmat_gradphi)          
                kmat(edofu,edofphi)=kmat(edofu,edofphi)+matmul(wmat_gradphi,rmat_gradphi)*jacw
                ! grad(w).rho*s
                ! This term is a transpose of the previous. Therefore, it will be
                ! computed later
                call compute_wmat_sPE(deriv,wmat_sPE)
                call compute_rmat_sPE(massdens_elmt(i,ielmt),interpf,rmat_sPE)
                kmat(edofphi,edofu)=kmat(edofphi,edofu)+matmul(wmat_sPE,rmat_sPE)*jacw
              endif
            endif
          endif
    
          if(ISPOT_DOF)then
            kmat(edofphi,edofphi)=kmat(edofphi,edofphi)+matmul(transpose(deriv),deriv)*jacw
            if(.not.ISDISP_DOF)then
              if(POT_TYPE==PMAGNETIC)then
                if(ismagnet_blk(imat))then
                  ! first compute the divergence of M
                  ! \nabla.M=dxMx+dyMy+dzMz
                  dxMx=dot_product(deriv(1,:),Mgll(1,:))
                  dyMy=dot_product(deriv(2,:),Mgll(2,:))
                  dzMz=dot_product(deriv(3,:),Mgll(3,:))
                  divM=dxMx+dyMy+dzMz
                  if(maxval(abs(M)).gt.ZERO)then
                    eload(edofphi)=eload(edofphi)+lagrange_gll(i,:)*divM*jacw
                  endif
                endif
              else !if(POT_TYPE==PMAGNETIC)
                eload(edofphi)=eload(edofphi)+lagrange_gll(i,:)*massdens_elmt(i,ielmt)*jacw
              endif
            endif
          endif
    
        endif
      enddo
      ! kmat terms for grad(w).rho*s
      !kmat(edofphi,edofu)=transpose(kmat(edofu,edofphi))
      if(ISDISP_DOF.and.ISPOT_DOF)then
        if(devel_nondim)then
          ! Note: PI*G is nondimensionalized
          kmat(edofphi,edofphi)=0.25_kreal*kmat(edofphi,edofphi)
        else
          kmat(edofphi,edofphi)=FOUR_PI_G_INV*kmat(edofphi,edofphi)
        endif
      endif
      if(.not.issymmetric(kmat))then
        write(*,*)'ERROR: matrix is unsymmetric!'
        stop
      endif
      !do i=1,nedof
      !  do j=1,nedof
      !    xval=kmat(i,j)                                                  
      !    if(ieee_is_nan(xval).or. .not.ieee_is_finite(xval))then                    
      !      write(*,*)'ERROR: stiffness matrix has nonfinite value/s!',myrank,ielmt,isinf,&
      !      mat_id(ielmt),xval,minval(abs(kmat)),maxval(abs(kmat))         
      !      stop                                                                     
      !    endif     
      !  enddo
      !enddo
      storekmat(:,:,ielmt)=kmat
      if(.not.ISDISP_DOF .and. ISPOT_DOF)then
        rhoload(egdof)=rhoload(egdof)+eload
      endif
    enddo ! i_elmt
    
    ! rhoload is computed if only the ISPOT_DOF is TRUE
    ! multiply rhoload by 4*PI*G
    if(.not.ISDISP_DOF.and.ISPOT_DOF)then
      if(.not.devel_nondim)then
        if(POT_TYPE==PMAGNETIC)then
          !rhoload=rhoload
        else
          rhoload=FOUR_PI_G*rhoload
        endif
      else
        if(POT_TYPE==PMAGNETIC)then
          !rhoload=rhoload
        else
          rhoload=FOUR*rhoload
          ! Note: PI*G is nondimensionalized
        endif
      endif
      rhoload(0)=ZERO
    endif
    deallocate(shape_infinite,dshape_infinite)
    deallocate(lagrange_gl,dlagrange_gl)
    deallocate(GLw)
    
    end subroutine compute_stiffness_elasticOLD
    !===============================================================================
    

! ####################### FUNCS FOR GETTING FS DETAILS #########################
subroutine get_fs_details(i_elmt, iface, nfgll, gw, dsq4)
  use global 
  use set_precision
  use integration 
  use free_surface
  implicit none 

  ! IO variables: 
  integer           :: i_elmt, nfgll ,iface 
  real(kind=kreal)  :: gw(:), dsq4(:,:,:)


  ! Face number (ie between 1 and 6) and get related properties
  iface = iface_fs(i_elmt)    
  if(iface==1 .or. iface==3)then
      nfgll             = ngllzx
      gw(1:nfgll)       = gll_weights_zx
      dsq4(:,:,1:nfgll) = dshape_quad4_zx

    elseif(iface==2 .or. iface==4)then
      nfgll             = ngllyz
      gw(1:nfgll)       = gll_weights_yz
      dsq4(:,:,1:nfgll) = dshape_quad4_yz

    elseif(iface==5 .or. iface==6)then
      nfgll             = ngllzx
      gw(1:nfgll)       = gll_weights_xy
      dsq4(:,:,1:nfgll) = dshape_quad4_xy
    else
      !write(errtag,'(a)')'ERROR: wrong face ID for traction!'
      return
  endif
end subroutine get_fs_details
!===============================================================================

subroutine calc_SL_stiffness(i_elmtfs,kmatSL)
    ! Uses 
    use global
    use set_precision
    use element
    use free_surface
    use integration
    use math_constants
    implicit none  
    
    ! IO 
    integer                        :: i_elmtfs        ! loops
    real(kind=kreal)               :: kmatSL(nedof,nedof)

    ! Local
    real(kind=kreal)               :: detjac2d        ! 2d jacobian
    integer                        :: iface           ! face ID for elmt 
    integer                        :: i_elmt          ! face ID for elmt 
    integer                        :: nfgll           ! ngll on 2D face
    real(kind=kreal), allocatable  :: gw(:)           ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4), face_normal(3),& 
                                    dx_dxi(NDIM), dx_deta(NDIM)
    integer :: num4(4), gid, phi_ind, u_ind, j,k, abg, xyg, gid_abg, gid_xyg, i_dim, i
    real(kind=kreal) ::   pi_2d_abg, pi_2d_xyg, area_inv, &
                        ival, iival, rho_over_g ,face_normal_len, cos_theta

    real(kind=kreal) :: g0abg, grav_abgj, Cabg, utfj, g0xyg, Cxyg, v1,rho_Ag, vertical(3), unit_normal(3)
    integer :: num(nenode)


    integer :: abgdof(5), xygdof(5), iloop, gidloc(maxngll2d)
    integer :: face_nodes(maxngll2d),  ggdof_elmt_fs(5,maxngll2d), & 
                dof_u(NDIM, maxngll2d), dof_u_tmp(NDIM,ngll), dof_phi(maxngll2d), dof_sl(maxngll2d)

    integer :: fgdof(nndof*maxngll2d)    , nfdof        

    ! Code
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))
    

    ! Inverse area
    area_inv = ONE/SLarea


      ! Get details of element's face that lies on free surface
    call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4)
    num4   = gnum4_fs(:, i_elmtfs)
    coord  = g_coord(:,num4)
    

    ! Node values within the element (1-27) - returns 9 points
    face_nodes = hexface(iface)%node

    ! The global IDs of the nodes? array of length maxngll2d
    gidloc =  gnum_fs(:, i_elmtfs)

    ! The GLOBAL DOFs values of the entire element  (5,maxngll2d)
    ggdof_elmt_fs(:,:) = ggdof(:, gnum_fs(:, i_elmtfs))
    !nfdof = nfgll * NNDOF

    ! Get the element-scale degrees of freedom  
    dof_u_tmp = reshape(edofu, (/NDIM, ngll/)) 
    dof_u     = dof_u_tmp(:,face_nodes)  !(NDIM, maxngll2d)
    dof_phi   = edofphi(face_nodes)      !(maxngll2d)
    dof_sl    = edofsl(face_nodes)       !(maxngll2d)



    do abg = 1, nfgll ! ABG
      gid_abg = gidloc(abg)            ! Global ID of node ABG
      g0abg   = ABS(g0_nodal(gid_abg)) ! g0 abg - in SL equation, this scalar g is positive

      ! Get Jacobian_2d x weights for ABG 
      dx_dxi  = matmul(coord,dshape4(1,:,abg))
      dx_deta = matmul(coord,dshape4(2,:,abg))
      face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
      face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
      face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

      face_normal_len =  sqrt(dot_product(face_normal,face_normal))
      pi_2d_abg     = gw(abg) * face_normal_len ! Weights*jacw

      ! Get the values here because they are repeated lots 
      Cabg       = oceanf(i_elmtfs, abg)  ! Ocean func abg


      ! We need the dot product of the vertical with the normal to the surface 
      ! The sign is not relevant becaause the water is always pushing down from the top surface
      !vertical    = zero
      !vertical(3) = one

      ! Normalise the length of the face normal to get unit normal to the free surface
      !unit_normal = face_normal/face_normal_len
      !cos_theta   = ABS(dot_product(unit_normal,vertical) ) 

      ! Factor of rho/g outside of integral 
      rho_over_g =  (rho_water/g0abg)       ! rho/g
      rho_Ag     =  (rho_over_g/SLarea)     ! rho/(g*Area)

      ! DIAGONAL theta_tilde theta_dot
      kmatSL(dof_sl(abg), dof_sl(abg)) = kmatSL(dof_sl(abg), dof_sl(abg)) &
                                       - (theta_tf * pi_2d_abg * g0abg * rho_water)

      ! COUPLING TERMS: 
      ! theta_tilde Phi_dot 
      kmatSL(dof_sl(abg), dof_phi(abg)) = kmatSL(dof_sl(abg), dof_phi(abg))  &
                                        - (g0abg * pi_2d_abg * theta_tf  * rho_over_g)

      ! phi_tilde Phi_dot 
      kmatSL(dof_phi(abg), dof_phi(abg)) = kmatSL(dof_phi(abg), dof_phi(abg)) &
                                         - (phi_tf * Cabg * pi_2d_abg  * rho_over_g)
      
 
      do j=1,NDIM
          ! Grav0_nodal is defined as a negative vector (i.e. -9.8), but grad Phi is positive
          ! g = - nabla Phi 
          grav_abgj = -grav0_nodal(j, gid_abg)

          ! u_tilde Phi_dot 
          kmatSL(dof_u(j, abg), dof_phi(abg)) = kmatSL(dof_u(j, abg), dof_phi(abg))& 
                                              - (Cabg * pi_2d_abg *  u_tf(j)*grav_abgj*rho_over_g) 

          ! theta_tilde u_dot 
          kmatSL(dof_sl(abg), dof_u(j, abg))  = kmatSL(dof_sl(abg), dof_u(j, abg))& 
                                              - (pi_2d_abg*g0abg*theta_tf*grav_abgj*rho_over_g)

          ! phi_tilde u_dot 
          kmatSL(dof_phi(abg), dof_u(j,abg)) = kmatSL(dof_phi(abg), dof_u(j,abg)) & 
                                             - (Cabg * pi_2d_abg * phi_tf * grav_abgj * rho_over_g)

          do k=1,NDIM
            ! u_tilde u_dot  
            ! NOTE THE NEGATIVE in grav0_nodal is needed for same reason as above 
            kmatSL(dof_u(k,abg),dof_u(j,abg)) = kmatSL(dof_u(k,abg),dof_u(j,abg))  & 
                                              - (Cabg * pi_2d_abg * grav_abgj *    & 
                                                u_tf(k) * (-grav0_nodal(k, gid_abg)) * rho_over_g)
          enddo !k

      enddo  ! j 


      ! Diagonal+non-diagonal components of Kmat coupling SL 
      ! Note here that ABG is the index of the variable
      ! and XYG is the test function so when we assemble the 
      ! matrix, ABG should be the column index 

      do xyg = 1, nfgll !XYG
        gid_xyg = gidloc(xyg)  
        g0xyg   =  ABS(g0_nodal(gid_xyg))      ! g0 abg 

        ! Get Jacobian_2d x weights for XYG 
        dx_dxi  = matmul(coord,dshape4(1,:,xyg))
        dx_deta = matmul(coord,dshape4(2,:,xyg))
        face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
        face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
        face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
        pi_2d_xyg  =  gw(xyg)*sqrt(dot_product(face_normal,face_normal)) ! Weights*jacw
        
        Cxyg       =  oceanf(i_elmtfs, xyg)  ! Ocean func abg

        ! SL_tilde, Phi_dot coupling 
        kmatSL(dof_sl(xyg), dof_phi(abg)) = kmatSL(dof_sl(xyg), dof_phi(abg)) + & 
                                            (g0xyg * theta_tf * pi_2d_abg * Cabg * pi_2d_xyg * rho_Ag)

        ! Phi_tilde, Phi_dot coupling 
        kmatSL(dof_phi(xyg),dof_phi(abg)) = kmatSL(dof_phi(xyg),dof_phi(abg)) + &
                                            (Cxyg * phi_tf * pi_2d_abg * Cabg * pi_2d_xyg * rho_Ag)

        do j=1,NDIM
          v1 = pi_2d_abg * Cabg * (-grav0_nodal(j, gid_abg)) * pi_2d_xyg 

          ! U_tilde, Phi_dot coupling  
          kmatSL(dof_u(j, xyg), dof_phi(abg)) = kmatSL(dof_u(j, xyg), dof_phi(abg)) & 
                                              + (pi_2d_abg * Cabg * pi_2d_xyg * Cxyg * u_tf(j) & 
                                                 *  (-grav0_nodal(j, gid_xyg)) * rho_Ag)

          ! SL_tilde, U_dot coupling
          kmatSL(dof_sl(xyg), dof_u(j, abg)) = kmatSL(dof_sl(xyg), dof_u(j, abg)) + (v1 * g0xyg * theta_tf * rho_Ag) 

          ! Phi_tilde, U_dot coupling  
          kmatSL(dof_phi(xyg), dof_u(j,abg)) = kmatSL(dof_phi(xyg), dof_u(j,abg)) + (v1 * Cxyg * phi_tf * rho_Ag)

          do k = 1, NDIM 
            ! U_tilde, U_dot coupling  
            kmatSL(dof_u(k,xyg), dof_u(j,abg)) = kmatSL(dof_u(k,xyg), dof_u(j,abg))  & 
                                                 + (v1 * Cxyg * u_tf(k) *  (-grav0_nodal(k, gid_xyg)) * rho_Ag)!*cos_theta
          enddo ! k 
        enddo ! j
      enddo! xyg

    enddo! abg
    deallocate(gw)
    deallocate(dshape4)
end subroutine calc_SL_stiffness
!===============================================================================
  

  
  subroutine check_KSL(contains_empty, start_index)
    ! Checks that the local SL stiffness matrix for an element does not have any empty rows 
    use set_precision
    use global 
    ! IO variables: 
    real(kind=kreal) :: row(nedof)
    logical :: contains_empty
    integer :: start_index
    ! Local: 
    integer :: i,j, rowctr


    contains_empty = .false. 
    
    do i = start_index, 4*ngll + maxngll2d !row loop  
      rowctr = 0 

      do j= 1, nedof
        if (ksl(i,j).ne.0) then 
          rowctr = rowctr + 1 
        endif 
      enddo 

      if (rowctr.eq.0)then 
        !write(*,*)'   EMPTY ROW, ID: ', i
        contains_empty = .true.
      endif 

    enddo ! row loop 

  end subroutine check_KSL


    end module matrix_vector
    !===============================================================================
