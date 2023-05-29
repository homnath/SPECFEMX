! contains scripts originally in SPECFEMX for nondimensionalisation of problem 
! Variables are non-dimensionalised to avoid rounding/Sig Fig errors when using large numbers like Earth radius etc 

! Last modified:  8th June 2022 (WE) 

! (Non)dimensionalize
module nondimension
use set_precision
contains 
!_______________________________________________________________________________

! Set nondimensionlization reference variables.
subroutine set_nondimension_refs()

! USES 
use global 

#if(USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif

implicit none

!-------------------------------------------------------------------------------
! set dimensionalize parameters
! minimum, maximum density
! density can be negative for gravity anomaly calculation
! massdens_elmt is not allocated for the magnetic anomaly computation

    if(myrank.eq.0)then 
      write(*,*)'---------- Nondimensionalisation parameters: ----------'
    endif 

    if(allocated(massdens_elmt))then
      if(infbc)then
        mindensity=minscal(minval(massdens_elmt(:,elmt_finite)))
        maxdensity=maxscal(maxval(massdens_elmt(:,elmt_finite)))
      else
        mindensity=minscal(minval(massdens_elmt))
        maxdensity=maxscal(maxval(massdens_elmt))
      endif

      ! Always use positive value for nondimensionalizing
      ! It may be that water is the largest density value
      maxdensity=max(abs(mindensity),abs(maxdensity))

      if(myrank.eq.0)then 
        write(*, '(a,g0.6)')'* Min. density from model  : ', mindensity
        write(*, '(a,g0.6)')'* Max. density from model  : ', maxdensity
      endif 

      
      if(ISSL_DOF)then
        maxdensity=max(maxdensity, rho_water_dim)
        if(myrank==0)then
          write(*, '(a,g0.6)')'* Water density            : ', rho_water_dim
          write(*, '(a,g0.6)')'  ----> Using max density  : ', maxdensity
          write(*,*)
        endif 
      endif 
    endif


    ! minimum, maximum bulk modulus
    ! bulkmod_elmt is not allocated for the magnetic anomaly computation
    if(allocated(bulkmod_elmt))then
      if(infbc)then
        minbulkmod=minscal(minval(bulkmod_elmt(:,elmt_finite)))
        maxbulkmod=maxscal(maxval(bulkmod_elmt(:,elmt_finite)))
      else
        minbulkmod=minscal(minval(bulkmod_elmt))
        maxbulkmod=maxscal(maxval(bulkmod_elmt))
      endif
      if(myrank==0)then
        write(*,'(a,g0.6,1x,g0.6)')'* Bulk modulus range (N/m2): ',minbulkmod,maxbulkmod
      endif
    endif
    ! minimum, maximum shear modulus
    ! shearmod_elmt is not allocated for the magnetic anomaly computation
    if(allocated(shearmod_elmt))then
      if(infbc)then
        minshearmod=minscal(minval(shearmod_elmt(:,elmt_finite)))
        maxshearmod=maxscal(maxval(shearmod_elmt(:,elmt_finite)))
      else
        minshearmod=minscal(minval(shearmod_elmt))
        maxshearmod=maxscal(maxval(shearmod_elmt))
      endif
      if(myrank==0)then
        write(logunit,'(a,g0.6,1x,g0.6)')'* Shear modulus range (N/m2): ',minshearmod,maxshearmod
      endif
    endif

    ! Add empty line for clarity 
    if(myrank==0)then 
      write(*,*)
    endif
    
    return 
end subroutine set_nondimension_refs 
!-------------------------------------------------------------------------------

! Calculate dimension/nondimension parameters.
subroutine calc_nondimension_pars

  ! USES
  use global
  use math_library_mpi
  use nondimensionpar
  implicit none 
  ! IO variables
  ! Local variables  



  if(.not.devel_nondim)then
    ! DO NOT nondimensionalize
    if(myrank==0)then
      write(*,*)'*****   Nondimensionalise: NO  *****'
    endif
    DIM_DENSITY=ONE
    NONDIM_DENSITY=ONE
  
    DIM_L=ONE
    NONDIM_L=ONE
  
    NONDIM_T=ONE
    DIM_T=ONE
  
    DIM_VEL=ONE
    DIM_ACCEL=ONE
  
    NONDIM_ACCEL=ONE

    DIM_M=ONE
  
    DIM_MOD = ONE
    NONDIM_MOD=ONE
  
    DIM_MTENS=ONE
    NONDIM_MTENS=ONE
  
    DIM_GPOT=ONE
    DIM_G=ONE
  
    DIM_MPOT=ONE
    DIM_B=ONE

    DIM_EPOT=ONE

    DIM_ICELOAD = ONE
  else
    ! nondimensionalize
    DIM_DENSITY=maxdensity                               
    NONDIM_DENSITY=ONE/DIM_DENSITY                               
  
    DIM_L=absmaxcoord
    NONDIM_L=ONE/DIM_L
  
    NONDIM_T=sqrt(PI*GRAV_CONS*maxdensity)
    DIM_T=ONE/NONDIM_T

    DIM_VEL=DIM_L*NONDIM_T
    NONDIM_VEL=ONE/DIM_VEL

    DIM_ACCEL=DIM_VEL*NONDIM_T
    NONDIM_ACCEL=ONE/DIM_ACCEL

    DIM_M=maxdensity*DIM_L*DIM_L*DIM_L
  
    DIM_MOD = DIM_M*NONDIM_L*NONDIM_T*NONDIM_T
    NONDIM_MOD=ONE/DIM_MOD                             
  
    DIM_MTENS=DIM_DENSITY*(DIM_L**5)*NONDIM_T*NONDIM_T
    NONDIM_MTENS=ONE/DIM_MTENS  
  
    DIM_GPOT=PI*GRAV_CONS*maxdensity*DIM_L*DIM_L
    DIM_G=PI*GRAV_CONS*maxdensity*DIM_L
  
    ! F = ma so mass * acceleration?
    DIM_ICELOAD = DIM_M*DIM_ACCEL 



    if(myrank.eq.0)then
      write(*,*)'*****   Nondimensionalise: YES  *****'
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL MASS         : ', DIM_M
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL DENSITY      : ', DIM_DENSITY
      write(*,'(a,g0.6,1x,g0.6)')'* NON-DIMENSONAL DENSITY      : ', NONDIM_DENSITY
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL LENGTH       : ', DIM_L
      write(*,'(a,g0.6,1x,g0.6)')'* NON-DIMENSONAL LENGTH       : ', NONDIM_DENSITY
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL TIME         : ', DIM_T
      write(*,'(a,g0.6,1x,g0.6)')'* NON-DIMENSONAL TIME         : ', NONDIM_T
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL VELOCITY     : ', DIM_VEL
      write(*,'(a,g0.6,1x,g0.6)')'* NON-DIMENSONAL VELOCITY     : ', NONDIM_VEL
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL ACCELERATION : ', DIM_ACCEL
      write(*,'(a,g0.6,1x,g0.6)')'* NON-DIMENSONAL ACCELERATION : ', NONDIM_ACCEL
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL ELASTIC MOD  : ', DIM_MOD
      write(*,'(a,g0.6,1x,g0.6)')'* NON-DIMENSONAL ELASTIC MOD  : ', NONDIM_MOD
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL GRAVITY POT. : ', DIM_GPOT
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSIONAL GRAVITY     : ', DIM_G
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL MOMENT TENSOR: ', DIM_MTENS
      write(*,'(a,g0.6,1x,g0.6)')'* NON-DIMENSONAL MOMENT TENSOR: ', NONDIM_MTENS
      write(*,*)
      write(*,'(a,g0.6,1x,g0.6)')'*     DIMENSONAL ICE LOAD     : ', DIM_ICELOAD
      write(*,*)
    endif
  endif

  end subroutine calc_nondimension_pars
!-------------------------------------------------------------------------------

  subroutine apply_nondimension()
    ! USES
    use global !, only: g_coord, NONDIM_L, ISDISP_DOF, massdens_elmt, &
                !      NONDIM_DENSITY, bulkmod_elmt, NONDIM_MOD, shearmod_elmt, &
                !      pole_coord0, pole_coord1, NONDIM_L, axis_range
    use nondimensionpar
    ! IO variables
    ! Local variables
    ! Code: 

    ! Nondimensionlize
    g_coord=g_coord*NONDIM_L

    if(ISDISP_DOF)then
      massdens_elmt = massdens_elmt * NONDIM_DENSITY
      bulkmod_elmt  = bulkmod_elmt  * NONDIM_MOD
      shearmod_elmt = shearmod_elmt * NONDIM_MOD

      ym_blk        = ym_blk  * NONDIM_MOD
      coh_blk       = coh_blk * NONDIM_MOD
      rho_blk       = rho_blk * NONDIM_DENSITY
      gam_blk       = gam_blk * (NONDIM_DENSITY*NONDIM_ACCEL)
    endif

    if(IS_SL)then
      rho_water    = rho_water_dim*NONDIM_DENSITY
      if(myrank.eq.0)then
        write(*, '(a,g0.6)')'  ----> Using water density   : ', rho_water
      endif
    endif 

    if(IS_ICE)then
      rho_ice    = rho_ice_dim*NONDIM_DENSITY
      if(myrank.eq.0)then
        write(*, '(a,g0.6)')'  ----> Using ice density     : ', rho_ice
        write(*,*)
      endif
    endif 
    
    pole_coord0=pole_coord0*NONDIM_L
    pole_coord1=pole_coord1*NONDIM_L
    axis_range=axis_range*NONDIM_L

  end subroutine apply_nondimension
!-------------------------------------------------------------------------------
end module nondimension
!===============================================================================
