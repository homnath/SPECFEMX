module initialise_arrays 
use shared,only:check_allocate
implicit none
contains 
!_______________________________________________________________________________

subroutine initialise_global_arrays()
use global
use math_constants,only:ZERO
implicit none 

! Local variables
integer :: istat 
allocate(gdof(nndof,nnode),gdof_elmt(nedof,nelmt),stat=istat)
if (istat/=0)then
    write(*,*)'ERROR: cannot allocate memory!'
    stop
endif
allocate(bcnodalv(nndof,nnode),stat=istat)
if (istat/=0)then
    write(*,*)'ERROR: cannot allocate memory!'
    stop
endif
allocate(nodalu(nndofu,nnode),stat=istat)
if (istat/=0)then
    write(*,*)'ERROR: cannot allocate memory!'
    stop
endif
if(ISDISP_DOF)then
  if(savedata%stress.or.isplastic)then
      allocate(stress_elmt(nst,ngll,nelmt),stress_nodal(nst,nnode))
      stress_elmt=ZERO
  endif
  if(savedata%strain)then
      allocate(strain_elmt(nst,ngll,nelmt),strain_nodal(nst,nnode))
      strain_elmt=ZERO
  endif
endif

if(ISSL_DOF)then
    allocate(nodalustore(nndofu,nnode), stat=istat)
    if(istat/=0)then
        write(*,*)'ERROR: cannot allocate memory!'
        stop
    endif
endif 

if(ISPOT_DOF)then
    allocate(nodalphi(nnode),nodalg(ndim,nnode),nodalphistore(nnode), &
    nodalB(ndim,nnode),stat=istat)
    if(istat/=0)then
    write(*,*)'ERROR: cannot allocate memory!'
    stop
    endif
endif

return 
end subroutine initialise_global_arrays
!-------------------------------------------------------------------------------

subroutine initialise_local_arrays()
use global,only:ndim,ngnode,nenode,nedof,nedofu,nedofphi,nst,ISDISP_DOF,ISPOT_DOF
use local
implicit none 

! Local variables
integer :: istat 

allocate(num(nenode),coord(ngnode,ndim),jac(ndim,ndim),deriv(ndim,nenode),     &
bmat(nst,nedofu),kmat(nedof,nedof),eld(nedofu),bload(nedofu),vload(nedofu),   &
eload(nedofu),stat=istat)
allocate(egdof(nedof),stat=istat)
if(ISDISP_DOF)then
  allocate(egdofu(nedofu),stat=istat)
endif
if(ISPOT_DOF)then
  allocate(egdofphi(nedofphi),stat=istat)
endif
if(istat/=0)then
  write(*,*)'ERROR: cannot allocate memory!'
  stop
endif

return 
end subroutine initialise_local_arrays
!-------------------------------------------------------------------------------

subroutine initialise_equation_arrays()

use global 
use math_constants
use set_precision
use shared,only:check_allocate
implicit none 

integer :: istat
character(len=500) :: errsrc

errsrc='initialise_equation_arrays'

! HERE IS ALLOCATION OF U VECTOR 
allocate(storekmat(nedof,nedof,nelmt),storemmat(nedof,nelmt),stat=istat)
allocate(load(0:neq),bodyload(0:neq),selfload(0:neq),viscoload(0:neq), &
resload(0:neq),du(0:neq),u(0:neq),olddu(0:neq),                   &
slipload(0:neq),extload(0:neq),rhoload(0:neq),ubcload(0:neq),          &
iceload(0:neq), stat=istat)

call check_allocate(istat,errsrc)
if(istat/=0)then
  write(logunit,*)'ERROR: cannot allocate memory!',istat
  flush(logunit)
  stop
endif

! Initialise more loads and u vector
elas_e0   = ZERO
visco_q0  = ZERO
nodalu    = ZERO
bodyload  = ZERO
selfload  = ZERO
viscoload = ZERO
slipload  = ZERO ! slip load
ubcload   = ZERO
load      = ZERO
u         = ZERO
extload   = ZERO
rhoload   = ZERO

if(ISSL_DOF)then
  nodalustore    = ZERO 
  nodalphistore  = ZERO
endif

end subroutine initialise_equation_arrays
!-------------------------------------------------------------------------------

end module initialise_arrays
!===============================================================================
