! contains scripts originally in SPECFEMX for nondimensionalisation of problem 
! Variables are non-dimensionalised to avoid rounding/Sig Fig errors when using large numbers like Earth radius etc 

! Last modified:  8th June 2022 (WE) 


module nondimensionalisation 
contains 


subroutine set_nondimensional_params()

! USES 
use global 

#if(USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif

implicit none
! IO variables 

! Local variables 


!-------------------------------------------------------------------------------
! set dimensionalize parameters
! minimum, maximum density
! density can be negative for gravity anomaly calculation
! massdens_elmt is not allocated for the magnetic anomaly computation
    if(allocated(massdens_elmt))then
        if(infbc)then
          mindensity=minscal(minval(massdens_elmt(:,elmt_finite)))
          maxdensity=maxscal(maxval(massdens_elmt(:,elmt_finite)))
        else
          mindensity=minscal(minval(massdens_elmt))
          maxdensity=maxscal(maxval(massdens_elmt))
        endif
        if(myrank==0)then
          write(logunit,'(a,g0.6,1x,g0.6)')'min, max density (kg/m3): ',mindensity,maxdensity
          flush(logunit)
        endif
        ! Always use positive value for nondimensionalizing
        maxdensity=max(abs(mindensity),abs(maxdensity))
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
          write(logunit,'(a,g0.6,1x,g0.6)')'min, max bulkmod (N/m2): ',minbulkmod,maxbulkmod
          flush(logunit)
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
          write(logunit,'(a,g0.6,1x,g0.6)')'min, max shearmod (N/m2): ',minshearmod,maxshearmod
          flush(logunit)
        endif
      endif
      
      return 

end subroutine set_nondimensional_params 


subroutine calc_nondimensionalisation_vals

  ! USES
  use global
  use math_library_mpi
  use dimensionless
  implicit none 
  ! IO variables
  ! Local variables  



  if(.not.devel_nondim)then
    ! DO NOT nondimensionalize
    if(myrank==0)then
      write(logunit,*)'nondimensionalize: NO'
      flush(logunit)
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
  else
    ! nondimensionalize
    if(myrank==0)then
      write(logunit,*)'nondimensionalize: YES'
      flush(logunit)
    endif
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
  endif

  

  end subroutine calc_nondimensionalisation_vals


  subroutine apply_nondimensionalisation()
    ! USES
    use global !, only: g_coord, NONDIM_L, ISDISP_DOF, massdens_elmt, &
                !      NONDIM_DENSITY, bulkmod_elmt, NONDIM_MOD, shearmod_elmt, &
                !      pole_coord0, pole_coord1, NONDIM_L, axis_range
    use dimensionless
    ! IO variables
    ! Local variables
    ! Code: 

    ! Nondimensionlize
    g_coord=g_coord*NONDIM_L
    if(ISDISP_DOF)then
      massdens_elmt=massdens_elmt*NONDIM_DENSITY
      bulkmod_elmt=bulkmod_elmt*NONDIM_MOD
      shearmod_elmt=shearmod_elmt*NONDIM_MOD

      ym_blk=ym_blk*NONDIM_MOD
      coh_blk=coh_blk*NONDIM_MOD
      rho_blk=rho_blk*NONDIM_DENSITY
      gam_blk=gam_blk*(NONDIM_DENSITY*NONDIM_ACCEL)
    endif

    if(IS_SL)then
      SL0_constant = SL0_constant*NONDIM_L
      rho_water    = rho_water_dim*NONDIM_DENSITY
    endif 
    
    pole_coord0=pole_coord0*NONDIM_L
    pole_coord1=pole_coord1*NONDIM_L
    axis_range=axis_range*NONDIM_L



  end subroutine apply_nondimensionalisation
end module nondimensionalisation 