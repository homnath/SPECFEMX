module initialise_arrays 

    implicit none
    contains 



subroutine initialise_RHS_vectors(load, bodyload, selfload, viscoload, &
                                  resload, du, u, kmat, storekmat,     &
                                  storemmat, rhoload, ubcload, nodalu, & 
                                  visco_q0, elas_e0, extload, iceload)

use global 
use math_constants
use set_precision
implicit none 

real(kind=kreal),allocatable :: slipload(:), extload(:), bodyload(:),  &
                                selfload(:), viscoload(:),ubcload(:),  &
                                load(:), resload(:), kmat(:,:), du(:), & 
                                u(:), storekmat(:,:,:), storemmat(:,:),&
                                rhoload(:), nodalu(:,:), iceload(:),   & 
                                visco_q0(:,:,:,:), elas_e0(:,:,:)     
                                

integer :: istat


! HERE IS ALLOCATION OF U VECTOR 
allocate(load(0:neq),bodyload(0:neq),selfload(0:neq),viscoload(0:neq), &
resload(0:neq),du(0:neq),u(0:neq),kmat(nedof,nedof),                   &
storekmat(nedof,nedof,nelmt), storemmat(nedof,nelmt),                  &
slipload(0:neq),extload(0:neq),rhoload(0:neq),ubcload(0:neq),          &
iceload(0:neq), stat=istat)


if(istat/=0)then
  write(logunit,*)'ERROR: cannot allocate memory!'
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

end subroutine initialise_RHS_vectors

end module initialise_arrays