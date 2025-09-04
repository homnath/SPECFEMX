! relaxation_time.f90
! Last edit: WE Jun 2 2022

module relaxation_time

contains 

!_________________________________________________________
subroutine calc_relaxation_time(tunitfac)
use math_constants !,only:ONE, inftol
#if (USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif
use conversion_constants! , only:  SEC2HOUR, SEC2DAY, SEC2MONTH, SEC2YEAR
                    
                    
use global! , only: nmaxwell, devel_example, imatve_to_imat, &
            !        nmatblk_viscoelas, devel_rtfac, & 
            !        viscosity_blk, shearmod_blk, tunit, logunit, &
            !        myrank, ym_blk


  implicit none 

  ! In/Out variables: 

  ! Local variables: 
  integer          :: iviscoelas,imat,i_elmt,i_gll,i_maxwell,i_mat
  real(kind=kreal) :: tunitfac,ym
  real(kind=kreal) :: min_relaxtime,max_relaxtime

  !allocate(relaxtime(nmaxwell,nmatblk_viscoelas))
  allocate(relaxtime_elmt(nmaxwell,ngll,nelmt))

  ! Relaxation time
  ! Convert relaxation time unit to the time step time unit for consistency
  ! TODO: check for nondimensionalizing the time step and relaxation time
  ! Note: for viscosity it doesn't matter because the only relevant parameter is 
  ! 'tratio' which is dimensionless in itself.
  
  relaxtime_elmt=inftol
  tunitfac=ONE

  ! Can we use a case statement here instead? 
  if(index(tunit,'sec').gt.0 .or. index(tunit,'second').gt.0)then
  ! Do not convert, because original unit of relaxtime is second.
  tunitfac=ONE
  elseif(index(tunit,'hr').gt.0 .or. index(tunit,'hour').gt.0)then
  ! Convert second to hour
  tunitfac=SEC2HOUR
  elseif(index(tunit,'day').gt.0)then
  ! Convert second to day
  tunitfac=SEC2DAY
  elseif(index(tunit,'month').gt.0)then
  ! Convert second to month
  tunitfac=SEC2MONTH
  elseif(index(tunit,'yr').gt.0 .or. index(tunit,'year').gt.0)then
  ! Convert second to year
  tunitfac=SEC2YEAR
  else
  write(logunit,*)'ERROR: unrecognized time unit: ',trim(tunit)
  flush(logunit)
  endif


  ! compute relax time
  iviscoelas=0
  !do i_mat=1,nmatblk_viscoelas
  !!  if(mat_domain(i_mat)==VISCOELASTIC_DOMAIN .or. &
  !!    mat_domain(i_mat)==VISCOELASTIC_TRINFDOMAIN .or. &
  !!    mat_domain(i_mat)==VISCOELASTIC_INFDOMAIN)then
  !!    iviscoelas=iviscoelas+1
  !    imat=imatve_to_imat(i_mat)
  !    ! Variables viscosity_blk and shearmod_blk are NOT nondimensionalized.
  !    ! Therefore the relaxtime will be in seconds.
  !    ! If nondimensionalized it may be better to dimensionalize again to compute
  !    ! the relax time without any confusion.
  !    ! WRONG
  !    !relaxtime(:,iviscoelas)=tunitfac*viscosity_blk(:,iviscoelas)/shearmod_blk(imat)
  !    ! APPROXIMATE
  !    ! CORRECT


  !    if(trim(devel_example).eq.'axial_rod')then
  !    relaxtime(:,i_mat)=tunitfac*viscosity_blk(:,i_mat)/ym_blk(imat)
  !    else  
  !    !relaxtime(:,i_mat)=tunitfac*TWO*viscosity_blk(:,i_mat)/shearmod_blk(imat)
  !    print*,imat,shearmod_blk(imat)
  !    relaxtime(:,i_mat)=devel_rtfac*tunitfac*viscosity_blk(:,i_mat)/shearmod_blk(imat)
  !    endif
  !!  endif
  !enddo
  do i_elmt=1,nelmt
  !  if(mat_domain(i_mat)==VISCOELASTIC_DOMAIN .or. &
  !    mat_domain(i_mat)==VISCOELASTIC_TRINFDOMAIN .or. &
  !    mat_domain(i_mat)==VISCOELASTIC_INFDOMAIN)then
  !    iviscoelas=iviscoelas+1
      !imat=imatve_to_imat(i_mat)
      ! Variables viscosity_blk and shearmod_blk are NOT nondimensionalized.
      ! Therefore the relaxtime will be in seconds.
      ! If nondimensionalized it may be better to dimensionalize again to compute
      ! the relax time without any confusion.
      ! WRONG
      !relaxtime(:,iviscoelas)=tunitfac*viscosity_blk(:,iviscoelas)/shearmod_blk(imat)
      ! APPROXIMATE
      ! CORRECT

      !ym=9.0d0*bulkmod_elmt(i_gll,i_elmt)*shearmod_elmt(i_gll,i_elmt)/(THREE*bulkmod_elmt(i_gll,i_elmt)+shearmod_elmt(i_gll,i_elmt))
      !if(trim(devel_example).eq.'axial_rod')then
      !relaxtime(:,i_mat)=tunitfac*viscosity_blk(:,i_mat)/ym_blk(imat)
      !else  
      !!relaxtime(:,i_mat)=tunitfac*TWO*viscosity_blk(:,i_mat)/shearmod_blk(imat)
      !print*,imat,shearmod_blk(imat)
      !relaxtime(:,i_mat)=devel_rtfac*tunitfac*viscosity_blk(:,i_mat)/shearmod_blk(imat)
      !endif
      do i_gll=1,ngll
        do i_maxwell=1,nmaxwell
          print*,shearmod_elmt(i_gll,i_elmt)
          relaxtime_elmt(i_maxwell,i_gll,i_elmt)=devel_rtfac*tunitfac* &
          viscosity_elmt(i_maxwell,i_gll,i_elmt)/shearmod_elmt(i_gll,i_elmt)
        enddo
      enddo
  !  endif
  enddo

  min_relaxtime=minscal(minval(relaxtime_elmt))
  max_relaxtime=maxscal(maxval(relaxtime_elmt))
  
  if(myrank.eq.0)then
  write(logunit,'(a,g0.6,1x,a,g0.6)')'Relax time => min: ',min_relaxtime,' max: ',max_relaxtime
  write(logunit,'(a)')'Time unit: '//trim(tunit)
  flush(logunit)
  endif

return 
end subroutine calc_relaxation_time
!-------------------------------------------------------------------------------

end module relaxation_time
!===============================================================================
