! this module contains the routine for Mohr-Coulomb plasticity, Viscoplastic
! algorithm, and strength reduction techniques
! REVISION
!   HNG, Jul 07,2011; HNG, Apr 09,2010
module plastic_library
use set_precision
use math_constants
use conversion_constants,only:DEG2RAD,RAD2DEG

contains
! this function computes the stable pseudo-time step
! for viscoplastic algorithm (Cormeau 1975, Smith and Griffiths 2003)
function dt_viscoplas(nmatblk,nuf,phif,ymf,ismat) result(dt_min)
implicit none
integer,intent(in) :: nmatblk
! friction angle,Poisson's ratio, strength reduction factor
real(kind=kreal),intent(in) :: ymf(nmatblk),phif(nmatblk),nuf(nmatblk) ! phif in degrees
logical,optional,intent(in) :: ismat(nmatblk)
real(kind=kreal) :: dt_min
real(kind=kreal) :: dt,snphi
integer :: i_mat
logical :: ismat_on(nmatblk)
ismat_on=.true.
if(present(ismat))ismat_on=ismat
! compute minimum pseudo-time step for viscoplasticity
dt_min=inftol
do i_mat=1,nmatblk
  if(.not.ismat_on(i_mat))cycle
  snphi=sin(phif(i_mat)*deg2rad)
  dt=FOUR*(ONE+nuf(i_mat))*(ONE-TWO*nuf(i_mat))/(ymf(i_mat)*(ONE-TWO*nuf(i_mat)+ &
  snphi**2))
  if(dt<dt_min)dt_min=dt
enddo
end function dt_viscoplas
!=======================================================

subroutine strength_reduction(srf,phinu,nmatblk,coh,nu,phi,psi,cohf,nuf,phif,psif,&
istat)
implicit none
real(kind=kreal),intent(in) :: srf
logical :: phinu
integer,intent(in) :: nmatblk
real(kind=kreal),intent(in) :: coh(nmatblk),nu(nmatblk),phi(nmatblk),psi(nmatblk)
! phi,psi in degrees
real(kind=kreal),intent(out) :: cohf(nmatblk),nuf(nmatblk),phif(nmatblk),psif(nmatblk)
! phif,psif in degrees
! status whether the material properties has changed
integer,intent(out) :: istat

integer :: i_mat
real(kind=kreal) :: beta_nuphi,omtnu,snphi,snphif,tnphi,tnpsi

cohf=coh; nuf=nu; phif=phi; psif=psi
istat=0

! strength reduction
if(srf/=one)then
  do i_mat=1,nmatblk
    tnphi=tan(phi(i_mat)*deg2rad)
    phif(i_mat)=atan(tnphi/srf)*rad2deg
    tnpsi=tan(psi(i_mat)*deg2rad)
    psif(i_mat)=atan(tnpsi/srf)*rad2deg
    cohf(i_mat)=coh(i_mat)/srf
  enddo
endif

if(.not.phinu)return
! correction for phi-nu inequality sin(phi)>=1-2*nu
! Reference: Zheng et al 2005, IJNME
! currently only the Poisson's ratio is corrected
do i_mat=1,nmatblk
  omtnu=one-two*nu(i_mat)
  snphif=sin(phif(i_mat)*deg2rad)
  if(snphif<omtnu)then
    snphi=sin(phi(i_mat)*deg2rad)
    beta_nuphi=snphi/omtnu
    if(beta_nuphi<one)beta_nuphi=one+zerotol
    nuf(i_mat)=half*(one-snphif/beta_nuphi)
    istat=1 ! material properties has changed
    !print*,phi,phif
    !print*,nu,nuf
    !print*,snphif,one-two*nuf
    !print*,beta_nuphi
    !stop
  endif
enddo
return
end subroutine strength_reduction
!=======================================================

! this subroutine calculates the value of the yield function
! for a mohr-coulomb material (phi in degrees, theta in radians).
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine mohcouf(phi,c,sigm,dsbar,theta,f)
implicit none
real(kind=kreal),intent(in)::phi,c,sigm,dsbar,theta
real(kind=kreal),intent(out)::f
real(kind=kreal)::phir,snph,csph,csth,snth,r3=3.0_kreal
phir=phi*deg2rad
snph=sin(phir)
csph=cos(phir)
csth=cos(theta)
snth=sin(theta)
f=snph*sigm+dsbar*(csth/sqrt(r3)-snth*snph/r3)-c*csph
return
end subroutine mohcouf
!=======================================================

! this subroutine forms the derivatives of a mohr-coulomb potential
! function with respect to the three stress invariants
! (psi in degrees, theta in radians).
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine mohcouq(psi,dsbar,theta,dq1,dq2,dq3)
implicit none
 real(kind=kreal),intent(in)::psi,dsbar,theta
 real(kind=kreal),intent(out)::dq1,dq2,dq3
 real(kind=kreal)::psir,snth,snps,sq3,c1,csth,cs3th,tn3th,tnth,pt49=0.49_kreal,&
 pt5=0.5_kreal,r3=3.0_kreal

 psir=psi*deg2rad
 snth=sin(theta)
 snps=sin(psir)
 sq3=sqrt(r3)
 dq1=snps

 if(abs(snth).gt.pt49)then
   c1=one
   if(snth.lt.zero)c1=-one
   dq2=(sq3*pt5-c1*snps*pt5/sq3)*sq3*pt5/dsbar
   dq3=zero
 else
   csth=cos(theta)
   cs3th=cos(r3*theta)
   tn3th=tan(r3*theta)
   tnth=snth/csth
   dq2=sq3*csth/dsbar*((one+tnth*tn3th)+snps*(tn3th-tnth)/sq3)*pt5
   dq3=pt5*r3*(sq3*snth+snps*csth)/(cs3th*dsbar*dsbar)
 end if
return
end subroutine mohcouq
!=======================================================

! this subroutine forms the derivatives of the invariants with respect to
! stress in 2- or 3-d.
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
subroutine formm(stress,m1,m2,m3)
 implicit none
 real(kind=kreal),intent(in)::stress(:)
 real(kind=kreal),intent(out)::m1(:,:),m2(:,:),m3(:,:)
 real(kind=kreal)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,  &
   r3=3.0_kreal,r6=6.0_kreal,r9=9.0_kreal
 integer::nst,i,j
 nst=ubound(stress,1)
 select case(nst)
 case(4)
   sx=stress(1); sy=stress(2); sz=stress(4)
   txy=stress(3)

   dx=(two*sx-sy-sz)/r3; dy=(two*sy-sz-sx)/r3; dz=(two*sz-sx-sy)/r3
   sigm=(sx+sy+sz)/r3
   m1=zero; m2=zero; m3=zero
   m1(1,1:2)=one
   m1(2,1:2)=one
   m1(4,1:2)=one
   m1(1,4)=one
   m1(4,4)=one
   m1(2,4)=one
   m1=m1/r9/sigm
   m2(1,1)=two/r3
   m2(2,2)=two/r3
   m2(4,4)= two/r3
   m2(2,4)=-one/r3
   m2(4,2)=-one/r3
   m2(1,2)=-one/r3
   m2(2,1)=-one/r3
   m2(1,4)=-one/r3
   m2(4,1)=-one/r3
   m2(3,3)=two
   m3(3,3)=-dz
   m3(1:2,3)=txy/r3
   m3(3,1:2)=txy/r3
   m3(3,4)=-two*txy/r3
   m3(4,3)=-two*txy/r3
   m3(1,1)=dx/r3
   m3(2,4)=dx/r3
   m3(4,2)=dx/r3
   m3(2,2)=dy/r3
   m3(1,4)=dy/r3
   m3(4,1)=dy/r3
   m3(4,4)=dz/r3
   m3(1,2)=dz/r3
   m3(2,1)=dz/r3
 case(6)
   sx=stress(1);  sy=stress(2);  sz=stress(3)
   txy=stress(4); tyz=stress(5); tzx=stress(6)
   sigm=(sx+sy+sz)/r3
   dx=sx-sigm; dy=sy-sigm; dz=sz-sigm
   m1=zero; m2=zero
   m1(1:3,1:3)=one/(r3*sigm)
   do i=1,3
     m2(i,i)=two
     m2(i+3,i+3)=r6
   enddo
   m2(1,2)=-one
   m2(1,3)=-one
   m2(2,3)=-one
   m3(1,1)=dx
   m3(1,2)=dz
   m3(1,3)=dy
   m3(1,4)=txy
   m3(1,5)=-two*tyz
   m3(1,6)=tzx
   m3(2,2)=dy
   m3(2,3)=dx
   m3(2,4)=txy
   m3(2,5)=tyz
   m3(2,6)=-two*tzx
   m3(3,3)=dz
   m3(3,4)=-two*txy
   m3(3,5)=tyz
   m3(3,6)=tzx
   m3(4,4)=-r3*dz
   m3(4,5)=r3*tzx
   m3(4,6)=r3*tyz
   m3(5,5)=-r3*dx
   m3(5,6)=r3*txy
   m3(6,6)=-r3*dy
   do i=1,6
     do j=i+1,6
       m1(j,i)=m1(i,j)
       m2(j,i)=m2(i,j)
       m3(j,i)=m3(i,j)
     enddo
   enddo
   m1=m1/r3; m2=m2/r3; m3=m3/r3
 case default
   write(*,*)"ERROR: nst size not recognized in formm!"
   stop
 end select
return
end subroutine formm
!=======================================================

! this subroutine calculates the value of the yield function
! for a mohr-coulomb material (phi in degrees, theta in radians).
! this routine was copied and modified from
! Smith and Griffiths (2004): Programming the finite element method
SUBROUTINE mcdpl(phi,psi,dee,stress,pl)
!
! This subroutine forms the plastic stress/strain matrix
! for a Mohr-Coulomb material (phi,psi in degrees).
!
 IMPLICIT NONE
 REAL(kind=kreal),INTENT(IN)::stress(:),dee(:,:),phi,psi  
 REAL(kind=kreal),INTENT(OUT)::pl(:,:)
 REAL(kind=kreal),ALLOCATABLE::dfds(:),dqds(:),ddqds(:),dfdsd(:)
 REAL(kind=kreal)::t1,t2,t3,t4,t5,t6,t8,t10,t12,t13,t14,t15,t16,t17,t18,t19,t20, &
   t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,   &
   t38,t39,t40,t41,t42,t43,t44,t45,t46,t48,t49,t50,t51,t53,t54,t55,t56,   &
   t60,t61,t63,t64,t68,t69,t70,t71,t73,t74,t77,t79,t80,t82,t83,t84,t85,   &
   t86,t89,t92,t93,t94,t97,t98,t101,t103,t106,t110,t111,t113,t122,t129,   &
   t133,t140,t145,t152,t166,t186,t206,pm,pi,phir,snph,snth,sq3,sx,sy,sz,  &
   txy,tyz,tzx,zero=0.0_kreal,pt49=0.49_kreal,one=1.0_kreal,two=2.0_kreal,        &
   d3=3.0_kreal,d4=4.0_kreal,d6=6.0_kreal,d180=180.0_kreal,snps,psir
 REAL(kind=kreal)::denom
 INTEGER::i,j,ih
 
 ih=SIZE(stress)
 ALLOCATE(dfds(ih),dqds(ih),ddqds(ih),dfdsd(ih))
 pi=ACOS(-one) 
 phir=phi*pi/d180
 snph=SIN(phir)
 psir=psi*pi/d180
 snps=SIN(psir)
 sq3=SQRT(d3)
 SELECT CASE(ih)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(3) 
   sz=stress(4)   
   t3=one/d3
   t4=(sx+sy+sz)*t3
   t8=sz-t4
   t10=txy**2
   t16=(sx-sy)**2
   t18=(sy-sz)**2
   t20=(sz-sx)**2
   t25=SQRT((t16+t18+t20)/d6+t10)
   t26=t25**2
   t30=d3*sq3*((sx-t4)*(sy-t4)*t8-t8*t10)/two/t26/t25
   IF(t30>one)t30=one
   IF(t30<-one)t30=-one
   t31=ASIN(t30)
   t33=SIN(t31*t3)
   snth=-t33
   IF(ABS(snth).GT.pt49)THEN
     pm=-one
     if(snth.LT.zero)pm=one
     t2=snph/d3
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t14=SQRT((t4+t6+t8)*t10+t12)
     t17=one/t14/sq3
     t18=one/two
     t19=t17*t18
     t21=d3+pm*snph
     t23=two*sy
     t24=two*sz
     t31=two*sx
     dfds(1)=t2+t19*t21*(d4*sx-t23-t24)*t10/two
     dfds(2)=t2+t19*t21*(-t31+d4*sy-t24)*t10/two
     dfds(3)=t17*t18*t21*txy
     dfds(4)=t2+t19*t21*(-t23+d4*sz-t31)*t10/two
     t2=snps/d3
     t21=d3+pm*snps
     dqds(1)=t2+t19*t21*(d4*sx-t23-t24)*t10/two
     dqds(2)=t2+t19*t21*(-t31+d4*sy-t24)*t10/two
     dqds(3)=t17*t18*t21*txy
     dqds(4)=t2+t19*t21*(-t23+d4*sz-t31)*t10/two
   ELSE
     t1=one/d3
     t2=snph*t1
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=(t4+t6+t8)*t10+t12
     t14=SQRT(t13)
     t16=d3*sq3
     t18=(sx+sy+sz)*t1
     t19=sx-t18
     t20=sy-t18
     t21=t19*t20
     t22=sz-t18
     t25=t21*t22-t22*t12
     t26=one/two
     t28=t14**2
     t30=one/t28/t14
     t31=t16*t25*t26*t30
     IF(t31>one)t31=one
     IF(t31<-one)t31=-one
     t33=ASIN(t31)
     t34=t33*t1
     t35=COS(t34)
     t36=SIN(t34)
     t38=one/sq3
     t41=one/t14*(t35+t36*snph*t38)
     t43=two*sy
     t44=two*sz
     t46=(d4*sx-t43-t44)*t10
     t49=one-t1
     t53=t19*t1*t22
     t54=t21*t1
     t55=t1*t12
     t60=t16*t25
     t61=t28**2
     t64=t26/t61/t14
     t68=t16*(t49*t20*t22-t53-t54+t55)*t26*t30-d3/two*t60*t64*t46
     t70=d3**2
     t71=sq3**2
     t73=t25**2
     t74=two**2
     t77=t13**2
     t83=SQRT(one-t70*t71*t73/t74/t77/t13)
     t84=one/t83
     t85=t84*t1
     t89=t2*t38
     t94=two*sx
     t97=(-t94+d4*sy-t44)*t10
     t101=t1*t20*t22
     t111=t16*(-t101+t19*t49*t22-t54+t55)*t26*t30-d3/two*t60*t64*t97
     t129=-two*t16*t22*txy*t26*t30-d3*t60*t64*txy
     t140=(-t43+d4*sz-t94)*t10
     t152=t16*(-t101-t53+t21*t49-t49*t12)*t26*t30-d3/two*t60*t64*t140
     dfds(1)=t2+t41*t46/two+t14*(-t36*t68*t85+t35*t68*t84*t89)
     dfds(2)=t2+t41*t97/two+t14*(-t36*t111*t85+t35*t111*t84*t89)
     dfds(3)=t41*txy+t14*(-t36*t129*t85+t35*t129*t84*t89)
     dfds(4)=t2+t41*t140/two+t14*(-t36*t152*t85+t35*t152*t84*t89)
     t2=snps*t1
     t41=one/t14*(t35+t36*snps*t38)
     t89=t2*t38
     dqds(1)=t2+t41*t46/two+t14*(-t36*t68*t85+t35*t68*t84*t89)
     dqds(2)=t2+t41*t97/two+t14*(-t36*t111*t85+t35*t111*t84*t89)
     dqds(3)=t41*txy+t14*(-t36*t129*t85+t35*t129*t84*t89)
     dqds(4)=t2+t41*t140/two+t14*(-t36*t152*t85+t35*t152*t84*t89)
   END IF
 CASE(6)
   sx=stress(1)
   sy=stress(2)
   sz=stress(3)   
   txy=stress(4) 
   tyz=stress(5) 
   tzx=stress(6) 
   t3=one/d3
   t4=(sx+sy+sz)*t3
   t5=sx-t4
   t6=sy-t4
   t8=sz-t4
   t10=tyz**2
   t12=tzx**2
   t14=txy**2
   t23=(sx-sy)**2
   t25=(sy-sz)**2
   t27=(sz-sx)**2
   t32=SQRT((t23+t25+t27)/d6+t14+t10+t12)
   t33=t32**2
   t37=d3*sq3*(t5*t6*t8-t5*t10-t6*t12-t8*t14+two*txy*tyz*tzx)/two/t33/t32
   IF(t37>one)t37=one
   IF(t37<-one)t37=-one
   t38=ASIN(t37)
   t40=SIN(t38*t3)
   snth=-t40
   IF(ABS(snth).GT.pt49)THEN
     pm=-one
     IF(snth.LT.zero)pm=one  
     t2=snph/d3
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=tyz**2
     t14=tzx**2
     t16=SQRT((t4+t6+t8)*t10+t12+t13+t14)
     t19=one/t16/sq3
     t20=one/two
     t21=t19*t20
     t23=d3+pm*snph
     t25=two*sy
     t26=two*sz
     t33=two*sx
     t48=t20*t23
     dfds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
     dfds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
     dfds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
     dfds(4)=t19*t48*txy
     dfds(5)=t19*t48*tyz
     dfds(6)=t19*t48*tzx
     t2=snps/d3
     t23=d3+pm*snps
     t48=t20*t23
     dqds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
     dqds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
     dqds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
     dqds(4)=t19*t48*txy
     dqds(5)=t19*t48*tyz
     dqds(6)=t19*t48*tzx
   ELSE
     t1=one/d3
     t2=snph*t1
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=tyz**2
     t14=tzx**2
     t15=(t4+t6+t8)*t10+t12+t13+t14
     t16=SQRT(t15)
     t18=d3*sq3
     t20=(sx+sy+sz)*t1
     t21=sx-t20
     t22=sy-t20
     t23=t21*t22
     t24=sz-t20
     t29=two*txy
     t32=t23*t24-t21*t13-t22*t14-t24*t12+t29*tyz*tzx
     t33=one/two
     t35=t16**2
     t37=one/t35/t16
     t39=t18*t32*t33*t37
     IF(t39>one)t39=one
     IF(t39<-one)t39=-one
     t40=ASIN(t39)
     t41=t40*t1
     t42=COS(t41)
     t43=SIN(t41)
     t45=one/sq3
     t48=one/t16*(t42+t43*snph*t45)
     t50=two*sy
     t51=two*sz
     t53=(d4*sx-t50-t51)*t10
     t56=one-t1
     t60=t21*t1*t24
     t61=t23*t1
     t63=t1*t14
     t64=t1*t12
     t69=t18*t32
     t70=t35**2
     t73=t33/t70/t16
     t77=t18*(t56*t22*t24-t60-t61-t56*t13+t63+t64)*t33*t37-               &
       d3/two*t69*t73*t53
     t79=d3**2
     t80=sq3**2
     t82=t32**2
     t83=two**2
     t86=t15**2
     t92=SQRT(one-t79*t80*t82/t83/t86/t15)
     t93=one/t92
     t94=t93*t1
     t98=t2*t45
     t103=two*sx
     t106=(-t103+d4*sy-t51)*t10
     t110=t1*t22*t24
     t113=t1*t13
     t122=t18*(-t110+t21*t56*t24-t61+t113-t56*t14+t64)*t33*t37-           &
       d3/two*t69*t73*t106
     t133=(-t50+d4*sz-t103)*t10
     t145=t18*(-t110-t60+t23*t56+t113+t63-t56*t12)*t33*t37-               &
       d3/two*t69*t73*t133
     t166=t18*(-two*t24*txy+two*tyz*tzx)*t33*t37-d3*t69*t73*txy
     t186=t18*(-two*t21*tyz+t29*tzx)*t33*t37-d3*t69*t73*tyz
     t206=t18*(-two*t22*tzx+t29*tyz)*t33*t37-d3*t69*t73*tzx
     dfds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
     dfds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
     dfds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
     dfds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
     dfds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
     dfds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
     t2=snps*t1
     t48=one/t16*(t42+t43*snps*t45)
     t98=t2*t45
     dqds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
     dqds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
     dqds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
     dqds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
     dqds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
     dqds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
   END IF
 END SELECT
 ddqds=MATMUL(dee,dqds)
 dfdsd=MATMUL(dee,dfds) 
 denom=DOT_PRODUCT(dfdsd,dqds)
 DO i=1,ih
   DO j=1,ih
     pl(i,j)=ddqds(i)*dfdsd(j)/denom
   END DO
 END DO
 DEALLOCATE(dfds,dqds,ddqds,dfdsd)
RETURN
END SUBROUTINE mcdpl
!===============================================================================

end module plastic_library
!===============================================================================

