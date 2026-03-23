! This module contains preprocessing library routines.
! REVISION:
!   HNG, Jul 07,2011
module preprocess
implicit none
character(len=250),private :: myfname=" => preprocess.f90"
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

! This subrotine computes and stores elemental derivative and integration
! information.
subroutine precompute_derivative_integration(errcode,errtag)
use set_precision
use global,only:myrank,NDIM,nst,nelmt,ngll,nedof,nedofu,nedofphi,nenode,ngnode,&
ngllx,nglly,ngllz,ngll,g_coord,gdof_elmt,g_num,mat_domain,mat_id,massdens_elmt,&
bulkmod_elmt,shearmod_elmt,isempty_blk,rho_blk,ym_blk,magnetization_elmt,&
infinite_iface,infinite_face_idir,pole_coord0,pole_coord1,       &
pole_type,pole_axis,axis_range,ISDISP_DOF,ISPOT_DOF,    &
POT_TYPE,PGRAVITY,PMAGNETIC,    &
element_is_infinite,storederiv,storejw,storeinterpf_infinite,devel_nondim, &
isdxval,isdyval,isdzval,devel_gaminf,infquad, &
edofu,edofphi,grav0_nodal,dgrav0_elmt,ISGRAV0,&
imat_to_imatmag,magnetization_blk,ismagnet_blk,nelmt_infinite
use element,only:hex8_gnode,map2exodus_hex8
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

integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: i,i_gll
integer :: i_elmt,ielmt,imat
integer :: imatve,mdomain
integer :: num(nenode),egdof(nedof)
real(kind=kreal) :: cmat(nst,nst)
real(kind=kreal) :: detjac !determinant of Jacobian
real(kind=kreal) :: coord(ngnode,NDIM),jac(NDIM,NDIM),eload(nedof)
real(kind=kreal) :: dgmat(NDIM,NDIM),eg0(ngll,NDIM),g0(NDIM),dg0(6)

real(kind=kreal) :: interpf(NGLL),deriv(NDIM,nenode)

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
integer :: nip,nipinf,ielmt_infinite
logical :: isinf ! flag to check if the element if on the infinite domain
logical :: isfaces(6)

errtag="ERROR: unknown!"
errcode=-1
errsrc=trim(myfname)//' => precompute_derivative_integration'

! use ng for Gauss quadrature and ngll fro GLL-Radau quadrature
gaminf=2.00_kreal !1.99_kreal ! 2.0: X1 in the mid-position. gaminf must be > 1.0
if(devel_gaminf.gt.ONE .and. devel_gaminf.lt.FOUR)gaminf=devel_gaminf
nipinf=ngll

allocate(shape_infinite(nipinf,nginf),dshape_infinite(NDIM,nipinf,nginf))
allocate(lagrange_gl(nipinf,ngll),dlagrange_gl(NDIM,nipinf,ngll))

allocate(GLw(nipinf))
element_is_infinite=.false.
storederiv=zero
if(allocated(storeinterpf_infinite))storeinterpf_infinite=zero
! Purely elastic elements
! Viscoelastic elements are elastic at time = 0
! Following loops through nelmt_elas+nelmt_viscoelas
ielmt_infinite=0
do i_elmt=1,nelmt
  ielmt=i_elmt
  num=g_num(:,ielmt)
  imat=mat_id(ielmt)
  mdomain=mat_domain(imat) 
  coord=transpose(g_coord(:,num(hex8_gnode)))
  nip=ngll

  isfaces=infinite_iface(:,ielmt)
  isinf=any(isfaces)
  if(count(isfaces).gt.1)isinf=.false.
  element_is_infinite(ielmt)=isinf
  ! set coordinates
  if(isinf)then
    ielmt_infinite=ielmt_infinite+1
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
    storeinterpf_infinite(:,:,ielmt_infinite)=lagrange_gl
    
  endif

  egdof=gdof_elmt(:,i_elmt)
    
  do i=1,nip
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
      
      storederiv(:,:,i,ielmt)=deriv 
      jacw=detjac*GLw(i)
      storejw(i,ielmt)=jacw

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
      
      storederiv(:,:,i,ielmt)=deriv 
      jacw=detjac*gll_weights(i)
      storejw(i,ielmt)=jacw
    endif ! (isinf)
  enddo
enddo ! i_elmt

deallocate(shape_infinite,dshape_infinite)
deallocate(lagrange_gl,dlagrange_gl)
deallocate(GLw)

end subroutine precompute_derivative_integration
!===============================================================================

end module preprocess
!===============================================================================
