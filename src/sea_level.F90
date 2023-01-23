module sea_level 
    use global 
    use set_precision
    implicit none 

    contains 

    

! Probably not worth using as it is SLOW
subroutine calc_jacdet2d(i_gll, coord, dshape_quad4, detjac)
    use global 
    use set_precision
    implicit none 
    real(kind=kreal) :: face_normal(NDIM)

    integer :: i_gll

    real(kind=kreal) :: coord(:,:), dshape_quad4(:,:,:), dx_dxi(NDIM), dx_deta(NDIM)
    real(kind=kreal) :: detjac

    dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
    dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))

    ! Normal
    face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
    face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
    face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

    
    detjac=sqrt(dot_product(face_normal,face_normal))
end subroutine  calc_jacdet2d





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




subroutine get_fs_details_nodsq(i_elmt, iface, nfgll, gw)
    use global 
    use set_precision
    use integration 
    use free_surface
    implicit none 

    ! IO variables: 
    integer           :: i_elmt, nfgll ,iface 
    real(kind=kreal)  :: gw(:)


    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_elmt)    
    if(iface==1 .or. iface==3)then
        nfgll             = ngllzx
        gw(1:nfgll)       = gll_weights_zx

      elseif(iface==2 .or. iface==4)then
        nfgll             = ngllyz
        gw(1:nfgll)       = gll_weights_yz

      elseif(iface==5 .or. iface==6)then
        nfgll             = ngllzx
        gw(1:nfgll)       = gll_weights_xy
      else
        !write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        return
    endif
end subroutine get_fs_details_nodsq


subroutine get_fs_details_noweights(i_elmt, iface, nfgll)
    use global 
    use set_precision
    use integration 
    use free_surface
    implicit none 

    ! IO variables: 
    integer           :: i_elmt, nfgll ,iface 

    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_elmt)    
    if(iface==1 .or. iface==3)then
        nfgll             = ngllzx
      elseif(iface==2 .or. iface==4)then
        nfgll             = ngllyz
      elseif(iface==5 .or. iface==6)then
        nfgll             = ngllzx
      else
        !write(errtag,'(a)')'ERROR: wrong face ID for traction!'
        return
    endif
end subroutine get_fs_details_noweights
! #################### END GETTING FS DETAILS #########################






! ################### LOG AND OUTPUT FUNCTIONS  #######################
subroutine start_SL_log(errcode, errtag)

use global
#if(USE_MPI)
use mpi_library
#else 
use serial_library
#endif  
implicit none 

character(len=250) :: errtag
integer :: ios, errcode


! Create file 
if(myrank==0)then
    SL_log_file = trim(file_head)//'SL.log'
    open(unit=SLlogunit,file=trim(SL_log_file),status='replace',action='write',iostat=ios)
    if(ios.ne.0)then
        write(errtag,'(a)')'ERROR: cannot open log file: '//trim(SL_log_file)
        call control_error(errcode,errtag,stdout,myrank)
    endif
endif 

! Write confirmation of file: 
write(*, *)' Created SL log file'
write(SLlogunit, '(A,/,A)')' ****** CREATED SL LOG FILE ****** ', ' '

! Write summary of SL inputs: 

write(SLlogunit,*)
write(SLlogunit,*)'-----------------------------------------------------'
write(SLlogunit,*)'SL IS_SL         :  ', IS_SL
write(SLlogunit,*)'SL DOF           :  ', ISSL_DOF
write(SLlogunit,*)'Save original SL :  ', savedata%sl0
write(SLlogunit,*)'Save SL          :  ', savedata%sl
write(SLlogunit,*)'Save ocean func. :  ', savedata%oceanf
write(SLlogunit,*)'-----------------------------------------------------'
write(SLlogunit,*)


write(SLlogunit, *)'Sea level data read from:  ', trim(slfile)
if(IS_CART_SIM)then
    call summarise_SL_input_cart()
elseif(IS_GLOB_SIM)then 
    write(*,*) 'GLOBAL SIMULATIONS NOT IMPLEMETED YET'
    stop
else
    write(*,*) 'SIMULATION MUST BE GLOBAL OR CARTESIAN'
    stop
endif 
end subroutine start_SL_log




subroutine summarise_SL_input_cart()
    ! Summarises sea level data for cartesian sims
    use global 
    implicit none 

    write(SLlogunit,*)
    write(SLlogunit,*)'Model setup          : Cartesian'

    if (SL0_is_constant)then 
        write(SLlogunit,*)'Type of SL input     : Constant Z value'
        write(SLlogunit,*)'           value     : ', SL0_constant
        write(SLlogunit,*)
    endif 
end subroutine summarise_SL_input_cart




subroutine write_SL0_to_ensight(nodalsl)
    use global 
    use postprocess
    use set_precision
    use free_surface
#if(USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif
implicit none 

    real(kind=kreal) :: nodalsl(:)


    write(SLlogunit,*)'Saving the original SL values'
    write(SLlogunit,*)'  --> Min sea level: ', minscal(minval(nodalsl))
    write(SLlogunit,*)'  --> Max sea level: ', maxscal(maxval(nodalsl))
    

    ! On the free surface
    if(savedata%fsplot)then
      call write_scalar_to_file_freesurf(nnode_fs, nodalsl, &
      ext='sl0',istep=0) 
    endif
    
    if(savedata%fsplot_plane)then
      call write_scalar_to_file_freesurf(nnode_fs, nodalsl, &
      ext='sl0', istep=0,plane=.true.) 
    endif

    write(SLlogunit,*)'  ✓ Saved original sea level '
    write(SLlogunit,*)
    end subroutine write_SL0_to_ensight



subroutine write_OF_to_ensight()
    ! Writes the ocean function to ensight 
    use global 
    use postprocess
    use free_surface
    implicit none 

    write(SLlogunit,*)'Saving the Ocean function: '
    write(SLlogunit,*)' --> total oceanic nodes = ', INT(SUM(nodalOF)), '/', nnode_fs

    

    ! On the free surface
    if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalOF, &
        ext='oceanf',istep=0) 
    endif
    
    if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs, nodalOF, &
        ext='oceanf', istep=0, plane=.true.) 
    endif

    write(SLlogunit,*)'  ✓ Saved ocean function '
    write(SLlogunit,*)
end subroutine write_OF_to_ensight


! ################# END  LOG AND OUTPUT FUNCTIONS  ####################




! #################    INITIAL SETUP FUNCTIONS    #####################
subroutine prepare_sea_level(nodalsl)
    use global 
    use free_surface
    use set_precision
    use math_constants

    implicit none 

    integer :: istattemp, istat
    real(kind=kreal), allocatable :: nodalsl(:)

    write(SLlogunit, *)'Preparing sea level variables...'
    
    istat = 0

    ! If running SL then probably want ocean function: 
    if(IS_SL)then 
        write(SLlogunit,*)'  + number of unique FS nodes: ', nnode_fs
        allocate(oceanf(nelmt_fs, maxngll2d), nodalOF(nnode_fs), stat=istattemp) ! Allocate ocean function and nodal ocean func
        istat=istat+istattemp

        oceanf=ZERO 
        write(SLlogunit, *)'  --> Created ocean function'
    endif 

 
    if(ISSL_DOF)then 
        ! Allocate nodal sea level 
        allocate(nodalsl(nnode_fs), stat=istattemp)
        istat=istat+istattemp
    endif 

    ! Check allocations 
    if(istat/=0)then
        write(logunit,*)'ERROR: cannot allocate memory in prepare_sea_level!'
        flush(logunit)
        stop
    else 
        nodalsl = ZERO
    endif

    write(SLlogunit,*)'  ✓ Prepared sea level. '
    write(SLlogunit,*)

    return 
end subroutine prepare_sea_level




subroutine set_original_sea_level(nodalsl)
    ! Uses
    use set_precision
    use global 
    use integration
    use free_surface
    use math_constants
    implicit none 

    ! IO variables
    real(kind=kreal) :: nodalsl(:)
    ! Local variables 
    integer          :: i_elmt, iface, nfgll, i_gll 
    real(kind=kreal) :: theta, z_coord
   
    real(kind=kreal),dimension(:,:,:), allocatable :: ds_quad4
    real(kind=kreal),dimension(:),     allocatable :: gll_weight
    real(kind=kreal),dimension(:,:),   allocatable :: lag_gll
    real(kind=kreal),dimension(:,:,:), allocatable :: dlag_gll


    allocate(ds_quad4(2,4,maxngll2d))
    allocate(gll_weight(maxngll2d),lag_gll(maxngll2d,maxngll2d),  &
    dlag_gll(2,maxngll2d,maxngll2d))


    ! Code:
    nodalsl = 0.0_kreal

    ! For cartesian: 
    if(IS_CART_SIM)then
        if(SL0_is_constant)then

            ! SL0 is the z coordinate of the sea surface. We therefore
            ! need to calculate the theta value for each point based on the 
            ! the z coordinate of the face 

            ! Loop for each face on the surface: 
            do i_elmt=1, nelmt_fs  

                ! Now with integration weights, ngll etc for face: 
                ! needs nfgll at least
                call get_fs_details(i_elmt, iface, nfgll, gll_weight, ds_quad4)


                ! For each GLL point calculate and store SL0
                do i_gll = 1, nfgll 
                    ! Get Z coordinates for the face and global IDs 
                    z_coord = g_coord(3,  gnum_fs(i_gll,i_elmt))
                    theta   = SL0_constant - z_coord 
                    
                    if(theta.gt.0.0_kreal)then
                        nodalsl(rgnum_fs(i_gll, i_elmt)) = theta
                    endif 
                enddo 
            enddo 
        endif
    endif 

    deallocate(ds_quad4, gll_weight,lag_gll, dlag_gll)
end subroutine set_original_sea_level
! ################### END INITIAL SETUP FUNCTIONS  #####################






subroutine update_ocean_function(nodalice, nodalsl, errcode, errtag, use_orig)
! Routine checks each GLL point on the surface to see if it is part of the ocean set
! see Crawford et al 2018, eqn 31-32.
! Set contains any nodes in which rho_w * theta > rho_i * I 
! If use_orig then will use nodalsl0 and nodalice0 instead of nodalsl and nodalice 

! Bit of an issue here because we need the values of theta, I not their rates (time derivs)
use global 
use free_surface
implicit none 

! IO variables
character(len=250) :: errtag
integer :: ios, errcode
real(kind=kreal), allocatable :: nodalice(:), nodalsl(:)
logical :: use_orig


! Local variables 
integer          :: i_elmt, iface, numf(maxngll2d), i_numf, i_node, nfgll
integer          :: i_gll 
real(kind=kreal) :: theta, I


write(SLlogunit,*)'WARNING: NEED TO ACCURATELY IMPLEMENT OCEAN SET/OCEAN FUNCTION BASED ON ICE CONDITION - see Crawford et al 2018, eqn 31'

! We may want the ocean function but not running simulation 
! - in that case we need to use nodalsl0 rather than u 
if (use_orig)then 
    ! In this case we can calculate ocean function using nodalsl0 and nodalice0: 
    write(SLlogunit,*)'Calculating ocean function with initial values'

    do i_elmt=1, nelmt_fs  

        ! Face number (ie between 1 and 6) and get related properties
        call get_fs_details_noweights(i_elmt, iface, nfgll)

        ! Loop through GLL on the surface: 
        do i_gll = 1, nfgll
            theta =  nodalsl(rgnum_fs(i_gll, i_elmt))
            I     =  nodalice(rgnum_fs(i_gll, i_elmt))

            ! if rho_w theta > rho_ice I then part of ocean set
            if (rho_water * theta .GT. I * rho_ice) then 
                oceanf(i_elmt, i_gll) = 1.0_kreal
                nodalOF(rgnum_fs(i_gll, i_elmt)) = 1.0_kreal
            else
                oceanf(i_elmt, i_gll) = 0.0_kreal
                nodalOF(rgnum_fs(i_gll, i_elmt)) = 0.0_kreal
            endif 

        enddo 
    enddo
else
    write(*,*)'ERROR: UPDATING OCEAN FUNCTION NOT IMPLEMENTED FOR CALCULATIONS YET'
    stop
endif 

write(SLlogunit,*)'  ✓ Updated ocean function. '
write(SLlogunit,*)''

end subroutine update_ocean_function



subroutine calculate_SL_A()
    ! Calculates the area covered by ocean (integral of ocean func
    ! over the solid surface)
    use global
    use element
    use free_surface
    use integration
    use math_constants

    implicit none 

    integer                        :: i_elmtfs, i_gll ! loops
    real(kind=kreal)               :: detjac2d        ! 2d jacobian
    integer                        :: iface           ! face ID for elmt 
    integer                        :: i_elmt          ! face ID for elmt 
    integer                        :: nfgll           ! ngll on 2D face
    real(kind=kreal), allocatable  :: gw(:)           ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4), face_normal(3),& 
                                    dx_dxi(NDIM), dx_deta(NDIM)
    integer :: num4(4)

    ! Code
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))

    write(SLlogunit,*)
    write(SLlogunit,*)'Calculating ocean area'
    SLarea = ZERO 

    do i_elmtfs = 1, nelmt_fs

        call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4)
        
        num4   = gnum4_fs(:, i_elmtfs)
        coord = g_coord(:,num4)

        do i_gll = 1, nfgll 
            
            ! Calculate the magnitude of the 2D jacobian 
            dx_dxi  = matmul(coord,dshape4(1,:,i_gll))
            dx_deta = matmul(coord,dshape4(2,:,i_gll))

            ! Calc normal and therefore jac dec (2D) on the fly
            face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
            face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
            face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
            detjac2d=sqrt(dot_product(face_normal,face_normal))

            SLarea = SLarea + oceanf(i_elmtfs, i_gll)*gw(i_gll)*detjac2d
        enddo ! i_gll
    enddo   ! i_elmtfs

    write(SLlogunit,*)'  --> Area of ocean:    ', SLarea 
    write(SLlogunit,*)'  ✓ Calculated sea level area. '

    ! Escape if no water. 
    if(SLarea.lt.ZERO .or. SLarea.eq.ZERO)then 
        write(*,*)'WARNING: area of ocean = 0 -- NO WATER!!!' 
        write(*,*)'USING SEA LEVEL AREA = 1' 
        write(*,*)'SETTING THETA TF TO  = 0' 

        SLarea = 1 ! cant be zero otherwise divide by zero
        oceanf   = ZERO
        theta_tf = 0.0_kreal
    endif 


    deallocate(gw)
    deallocate(dshape4)
    
end subroutine calculate_SL_A




end module