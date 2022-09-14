module sea_level 
    use global 
    use set_precision
    implicit none 

    contains 

    


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





! #######################  GETTING FS DETAILS #########################
subroutine get_fs_details(i_face, iface, nfgll, gw, dsq4)
    use global 
    use set_precision
    use integration 
    use free_surface
    implicit none 

    ! IO variables: 
    integer           :: i_face, nfgll ,iface 
    real(kind=kreal)  :: gw(:), dsq4(:,:,:)


    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_face)    
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




subroutine get_fs_details_nodsq(i_face, iface, nfgll, gw)
    use global 
    use set_precision
    use integration 
    use free_surface
    implicit none 

    ! IO variables: 
    integer           :: i_face, nfgll ,iface 
    real(kind=kreal)  :: gw(:)


    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_face)    
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


subroutine get_fs_details_noweights(i_face, iface, nfgll)
    use global 
    use set_precision
    use integration 
    use free_surface
    implicit none 

    ! IO variables: 
    integer           :: i_face, nfgll ,iface 

    ! Face number (ie between 1 and 6) and get related properties
    iface = iface_fs(i_face)    
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




subroutine write_SL0_to_ensight()
    use global 
    use postprocess
    use free_surface
    implicit none 

    write(SLlogunit,*)'Saving the original SL values'
    write(SLlogunit,*)'  --> Min sea level: ', minval(nodalsl0)
    write(SLlogunit,*)'  --> Max sea level: ', maxval(nodalsl0)
    

    

    ! On the free surface
    if(savedata%fsplot)then
      call write_scalar_to_file_freesurf(nnode_fs, nodalsl0, &
      ext='sl0',istep=0) 
    endif
    
    if(savedata%fsplot_plane)then
      call write_scalar_to_file_freesurf(nnode_fs, nodalsl0, &
      ext='sl0', istep=0,plane=.true.) 
    endif

    write(SLlogunit,*)'  ✓ Saved original sea level '
    write(SLlogunit,*)
    end subroutine write_SL0_to_ensight
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
        allocate(oceanf(nelmt_fs, maxngll2d)) ! Allocate ocean function 
        oceanf=ZERO 
        write(SLlogunit, *)'  --> Created ocean function'
    endif 


    ! Create store for initial SL 
    if(savedata%sl0)then 
        allocate(nodalsl0(nnode_fs), stat=istattemp)
        nodalsl0 = ZERO
        write(SLlogunit, *)'  --> Created initial SL (nodal)'
        istat=istat+istattemp
    endif

    ! Need to move to ISSL_DOF when implemented properly
    allocate(QSL(maxngll2d, nelmt_fs),           stat=istattemp)! Rphi 
    istat=istat+istattemp
    allocate(storeRphi(maxngll2d, nelmt_fs),     stat=istattemp)! QSL 
    istat=istat+istattemp
    allocate(storeRu(NDIM, maxngll2d, nelmt_fs), stat=istattemp)! Ru 
    istat=istat+istattemp

    QSL        = ZERO
    storeRphi  = ZERO
    storeRu    = ZERO


    if(ISSL_DOF)then 
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
        
        write(SLlogunit, *)'  --> Created Ru matrix'
        write(SLlogunit, *)'  --> Created Rphi matrix'
        write(SLlogunit, *)'  --> Created QSL matrix'
        write(SLlogunit, *)'  --> Created sea level (nodalsl) vector:'
    endif


    write(SLlogunit,*)'  ✓ Prepared sea level. '
    write(SLlogunit,*)


    return 
end subroutine prepare_sea_level




subroutine set_original_sea_level()
    ! Uses
    use set_precision
    use global 
    use integration
    use free_surface
    use math_constants
    implicit none 

    ! IO variables

    ! Local variables 
    integer          :: i_face, iface, nfgll, i_gll 
    real(kind=kreal) :: theta, z_coord
   
    real(kind=kreal),dimension(:,:,:), allocatable :: ds_quad4
    real(kind=kreal),dimension(:),     allocatable :: gll_weight
    real(kind=kreal),dimension(:,:),   allocatable :: lag_gll
    real(kind=kreal),dimension(:,:,:), allocatable :: dlag_gll


    allocate(ds_quad4(2,4,maxngll2d))
    allocate(gll_weight(maxngll2d),lag_gll(maxngll2d,maxngll2d),  &
    dlag_gll(2,maxngll2d,maxngll2d))


    ! Code:
    nodalsl0 = 0.0_kreal

    ! For cartesian: 
    if(IS_CART_SIM)then
        if(SL0_is_constant)then

            ! SL0 is the z coordinate of the sea surface. We therefore
            ! need to calculate the theta value for each point based on the 
            ! the z coordinate of the face 

            ! Loop for each face on the surface: 
            do i_face=1, nelmt_fs  

                ! Now with integration weights, ngll etc for face: 
                call get_fs_details(i_face, iface, nfgll, gll_weight, ds_quad4)


                ! For each GLL point calculate and store SL0
                do i_gll = 1, nfgll 
                    ! Get Z coordinates for the face and global IDs 
                    z_coord = g_coord(3,  gnum_fs(i_gll,i_face))
                    theta   = SL0_constant - z_coord 
                    
                    if(theta.gt.0.0_kreal)then
                        nodalsl0(rgnum_fs(i_gll, i_face)) = theta
                    endif 
                enddo 
            enddo 
        endif
    endif 


    deallocate(ds_quad4, gll_weight,lag_gll, dlag_gll)
end subroutine set_original_sea_level

! ################### END INITIAL SETUP FUNCTIONS  #####################




real(kind=kreal) function gcal(gid, i_gll, i_elmt)
! Calculates G (see overleaf doc.)
use global, only: g0_nodal, grav0_nodal
use math_constants
use set_precision
implicit none

! IO variables
integer, intent(in) :: gid, i_gll, i_elmt 
! Local 
integer :: j ! loop vars
real(kind=kreal) :: theta_tf, phi_tf, u_tf(NDIM), jsum

! TFs are 1 unless running code tests
theta_tf = ONE
phi_tf = ONE
u_tf = ONE

jsum = ZERO ! initialise
do j = 1, NDIM
    jsum = jsum + u_tf(j)*grav0_nodal(j, gid)
enddo 

gcal = g0_nodal(gid)*theta_tf + (oceanf(i_elmt, i_gll)*(phi_tf + jsum))  
end function gcal




subroutine update_ocean_function(u, errcode, errtag, use_s0)
! Routine checks each GLL point on the surface to see if its theta
! Value is 0 or not - if it is 0 then the ocean function becomes 0 
! if SL is larger than 0 then ocean function is 1.
use global 
use free_surface
implicit none 

! IO variables
character(len=250) :: errtag
integer :: ios, errcode
real(kind=kreal), allocatable :: u(:)
logical :: use_s0


! Local variables 
integer          :: i_face, iface, numf(maxngll2d), i_numf, i_node, nfgll
integer          :: i_gll 
real(kind=kreal) :: theta


! We may want the ocean function but not running simulation 
! - in that case we need to use nodalsl0 rather than u 
if (use_s0)then 
    ! In this case we can calculate ocean function using nodalsl0: 
    write(SLlogunit,*)'Calculating ocean function with nodalSL0'

    do i_face=1, nelmt_fs  

        ! Face number (ie between 1 and 6) and get related properties
        call get_fs_details_noweights(i_face, iface, nfgll)


        ! Loop through GLL on the surface: 
        do i_gll = 1, nfgll
            theta = nodalsl0(rgnum_fs(i_gll, i_face))


            if (theta.gt.0.0 ) then 
                ! Sea level is not zero - ocean func is 1 
                oceanf(i_face, i_gll) = 1.0_kreal

            elseif(theta==0.0 ) then 
                ! Sea level is zero - ocean func is 0 
                oceanf(i_face, i_gll) = 0.0_kreal

                else 
                write(errtag,'(a)')'SEA LEVEL VALUE IS NEGATIVE!! '
                return
            endif 

        enddo 
    enddo
else
    write(*, *)'ERROR: UPDATING OCEAN FUNCTION NOT IMPLEMENTED FOR CALCULATIONS'
    stop
endif 

write(SLlogunit,*)'  ✓ Updated ocean function. '

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
    if(SLarea.eq.ZERO)then 
        write(*,*)'ERROR: area of ocean = 0 -- NO WATER!!!' 
        stop 
    endif 


    deallocate(gw)
    deallocate(dshape4)
    
end subroutine calculate_SL_A





subroutine calc_SL_LHS()
    ! Uses 
    use global
    use element
    use free_surface
    use integration
    use math_constants
    implicit none 

    ! The Q matrix is literally just the test function multiplied
    ! by the Jac 2D 

    integer                        :: i_elmtfs, i_gll ! loops
    real(kind=kreal)               :: detjac2d        ! 2d jacobian
    integer                        :: iface           ! face ID for elmt 
    integer                        :: i_elmt          ! face ID for elmt 
    integer                        :: nfgll           ! ngll on 2D face
    real(kind=kreal), allocatable  :: gw(:)           ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4), face_normal(3),& 
                                    dx_dxi(NDIM), dx_deta(NDIM)
    integer :: num4(4), gid, phi_ind, u_ind, j
    real(kind=kreal) :: theta_tf, pi_2d, iRsum, rphival, ruval

    ! Code
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))

    write(SLlogunit,*)
    write(SLlogunit,*)'Calculating LHS matrices'


    theta_tf = ONE


    do i_elmtfs = 1, nelmt_fs
        ! Get details of face
        call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4)
        num4   = gnum4_fs(:, i_elmtfs)
        coord  = g_coord(:,num4)

        ! Get internal sum used for R matrix calculations 
            call calc_R_internal(i_elmtfs, iface, nfgll, gw, dshape4, num4, iRsum)

        do i_gll = 1, nfgll 
            
            gid = gnum_fs(i_gll, i_elmtfs) ! Global ID of node

            ! Calculate the magnitude of the 2D jacobian 
            dx_dxi  = matmul(coord,dshape4(1,:,i_gll))
            dx_deta = matmul(coord,dshape4(2,:,i_gll))

            ! Calc normal and therefore jac dec (2D) on the fly
            face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
            face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
            face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
            detjac2d=sqrt(dot_product(face_normal,face_normal))

            pi_2d = theta_tf * gw(i_gll)*detjac2d ! Weights*jacw

            ! Map QSL diag matrix to vector
            ! QSL(rgnum_fs(i_gll, i_elmtfs)) = QSL(rgnum_fs(i_gll, i_elmtfs)) + theta_tf*pi_2d 
            ! NOTE THAT THE INTEGRAL IS MULTIPLIED BY RHO_W*g so need to multiply also at each GLL pt 
            ! Will assemble in the Petsc loops 
            ! I dont think the addition is needed
            QSL(i_gll, i_elmtfs) = QSL(i_gll, i_elmtfs) + (theta_tf*pi_2d)* g0_nodal(gid)*rho_water
            
            ! Calculate diag R_phi and R_u matrices: 
            rphival = pi_2d*(gcal(gid, i_gll, i_elmtfs) - oceanf(i_elmtfs, i_gll)*iRsum)
            storeRphi(i_gll, i_elmtfs) =    -(rho_water/g0_nodal(gid)) * rphival 



            do j=1,NDIM
                storeRu(j, i_gll, i_elmtfs) = -(rho_water/g0_nodal(gid))* rphival*grav0_nodal(j,gid)
            enddo 


        enddo! i_gll
    enddo   ! i_elmtfs

    write(SLlogunit,*)'  ✓ Calculated LHS'

    deallocate(gw)
    deallocate(dshape4)
end subroutine calc_SL_LHS





subroutine calc_SL_RHS(slload)
    use global 
    use free_surface 
    implicit none 

end subroutine calc_SL_RHS







subroutine calc_R_internal(i_elmtfs, iface, nfgll, gw, dshape4, num4, internalsum)
    ! Calculates internal summation in both R matrices 
    use global
    use element
    use free_surface
    use integration
    use math_constants
    implicit none 
    ! IO variables
    integer                        :: i_elmtfs  ! loops
    integer                        :: iface           ! face ID for elmt 
    integer                        :: nfgll           ! ngll on 2D face
    real(kind=kreal), allocatable  :: gw(:)           ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4)
    integer                        :: num4(4)
    ! RETURNS:
    real(kind=kreal)               :: internalsum

    ! Local 
    real(kind=kreal) :: face_normal(3), dx_dxi(NDIM), dx_deta(NDIM), detjac2d  
    integer ::  i_gll, gid

    ! Code
    
    internalsum = ZERO

    do i_gll = 1, nfgll 
        
        ! Calculate the magnitude of the 2D jacobian 
        dx_dxi  = matmul(coord,dshape4(1,:,i_gll))
        dx_deta = matmul(coord,dshape4(2,:,i_gll))

        ! Calc normal and therefore jac dec (2D) on the fly
        face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
        face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
        face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
        detjac2d=sqrt(dot_product(face_normal,face_normal))

        ! Get global ID of node: 
        gid = gnum_fs(i_gll, i_elmtfs)
        internalsum = internalsum + (gw(i_gll) * detjac2d * gcal(gid, i_gll, i_elmtfs))
    enddo! i_gll


    internalsum = internalsum/SLarea
end subroutine calc_R_internal 

    


end module