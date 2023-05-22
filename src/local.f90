! DESCRIPTION
!  These modules contain local variables. 
! DEVELOPER
!  Hom Nath Gharti, Queen's University
! REVISION
!  HNG, May 19,2023
! TODO


! This module local parameters/variables
module local
use set_precision
implicit none
! num: g_num for particular element.

integer         , allocatable :: num(:),egdof(:)
integer         , allocatable :: egdofu(:)
integer         , allocatable :: egdofphi(:)

!kmat: stiffness matrix for each element
real(kind=kreal), allocatable :: coord(:,:),deriv(:,:),jac(:,:), &
                                 bmat(:,:),kmat(:,:), &
                                 eld(:),eload(:),bload(:),vload(:)
end module local
!===============================================================================
