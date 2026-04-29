! This module computes source time/frequency function.
module source_function
use set_precision
use math_constants
contains
!-------------------------------------------------------------------------------
complex(kind=kreal) function source_frequency_function_complex(freq,hdur)

implicit none
! Type of source function
integer,parameter :: SFTYPE=0

real(kind=kreal),intent(in) :: freq,hdur
complex(kind=kreal) :: omegath

if(SFTYPE==0)then
  ! Heaviside function
  source_frequency_function_complex = cmplx(ONE,ZERO)
elseif(SFTYPE==1)then
  ! Only imaginary component
  omegath=cmplx(ZERO,TWO*freq*hdur)
  ! Source_frequency_function_complex = sin(omegath)/omegath 
  source_frequency_function_complex = exp(omegath) 
else
  write(*,*)'ERROR: invalid SFTYPE for source_frequency_function_complex!'
  stop
endif

end function source_frequency_function_complex
!===============================================================================

real(kind=kreal) function source_frequency_function(freq,hdur)

implicit none
! Type of source function
integer,parameter :: SFTYPE=0

real(kind=kreal),intent(in) :: freq,hdur
real(kind=kreal) :: omegath

if(SFTYPE==0)then
  ! Heaviside function
  source_frequency_function = ONE
elseif(SFTYPE==1)then
  omegath=TWO*freq*hdur
  source_frequency_function = sin(omegath)/omegath 
else
  write(*,*)'ERROR: invalid SFTYPE for source_frequency_function!'
  stop
endif

end function source_frequency_function
!===============================================================================

real(kind=kreal) function source_time_function(t,hdur)

implicit none

real(kind=kreal),intent(in) :: t,hdur


! Quasi Heaviside
source_time_function = 0.5d0*(1.0d0 + erf(t/hdur))
! comp_source_time_function = dexp(-(t/hdur)**2)/(dsqrt(PI)*hdur)

end function source_time_function
!===============================================================================

real(kind=kreal) function source_time_function_rickr(t,f0)

implicit none
real(kind=kreal),intent(in) :: t,f0

! Ricker
source_time_function_rickr = (ONE-TWO*PI*PI*f0*f0*t*t)*exp(-PI*PI*f0*f0*t*t)

!!! Another source time function they have called 'Ricker' in some old papers,
!!! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006)
!!! in order to benchmark those simulations, the following formula is needed.
! comp_source_time_function_rickr = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

end function source_time_function_rickr
!===============================================================================
end module source_function
!===============================================================================
