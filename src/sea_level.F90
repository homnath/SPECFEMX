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
write(SLlogunit,*)'Save SL          :  ', savedata%sl
write(SLlogunit,*)'Save original ice:  ', savedata%ice0
write(SLlogunit,*)'Save ice         :  ', savedata%ice
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


    subroutine write_ICE0_to_ensight()
        use global 
        use postprocess
        use free_surface
        implicit none 
    
        write(SLlogunit,*)'Saving the original ICE values'
        write(SLlogunit,*)'  --> Min ice level: ', minval(nodalice0)
        write(SLlogunit,*)'  --> Max ice level: ', maxval(nodalice0)
        
    
        ! On the free surface
        if(savedata%fsplot)then
          call write_scalar_to_file_freesurf(nnode_fs, nodalice0, &
          ext='ice0',istep=0) 
        endif
        
        if(savedata%fsplot_plane)then
          call write_scalar_to_file_freesurf(nnode_fs, nodalice0, &
          ext='ice0', istep=0,plane=.true.) 
        endif
    
        write(SLlogunit,*)'  ✓ Saved original ice level '
        write(SLlogunit,*)
        end subroutine write_ICE0_to_ensight


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
        allocate(oceanf(nelmt_fs, maxngll2d), stat=istattemp) ! Allocate ocean function 
        istat=istat+istattemp

        oceanf=ZERO 
        write(SLlogunit, *)'  --> Created ocean function'
    endif 


    ! Create store for initial SL 
    if(savedata%sl0)then 
        allocate(nodalsl0(nnode_fs), nodalice0(nnode_fs), stat=istattemp)
        nodalsl0 = ZERO
        nodalice0 = ZERO
        write(SLlogunit, *)'  --> Created initial SL (nodal)'
        istat=istat+istattemp
    endif

 


    if(ISSL_DOF)then 

        ! Allocate LHS matrices for SL to be assembled into stiffness
        ! matrix: 
        !   QSL     - coupling Q matrix (theta, theta integral)
        !   slc_uu  - U_dot,   U_tilde   coupling  
        !   slc_pu  - Phi_dot, U_tilde   coupling  
        !   slc_ut  - U_dot,   SL_tilde  coupling  

        !   slc_pp  - Phi_dot, Phi_tilde coupling  
        !   slc_pt  - Phi_dot, SL_tilde  coupling  
        !   slc_up  - U_dot,   Phi_tilde coupling  
        !   
        allocate(QSL(maxngll2d, nelmt_fs),                          &
                 slc_uu(NDIM, maxngll2d, NDIM, maxngll2d, nelmt_fs),&
                 slc_pu(NDIM, maxngll2d, maxngll2d, nelmt_fs),      &
                 slc_ut(NDIM, maxngll2d, maxngll2d, nelmt_fs),      &
                 slc_up(NDIM, maxngll2d, maxngll2d, nelmt_fs),      &
                 slc_pp(maxngll2d, maxngll2d, nelmt_fs),            &
                 slc_pt(maxngll2d, maxngll2d, nelmt_fs),            &
                 stat=istattemp) 

        ! Update error alloc status
        istat=istat+istattemp
                 
        ! Initialise LHS arrays
        QSL     = ZERO
        slc_uu  = ZERO
        slc_pu  = ZERO
        slc_ut  = ZERO
        slc_pp  = ZERO
        slc_pt  = ZERO
        slc_up  = ZERO

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
        ! Output confirmation to log. 
        write(SLlogunit, *)'  --> Created QSL matrix'
        write(SLlogunit, *)'  --> Created slc_uu matrix (U_dot, U_tilde)'
        write(SLlogunit, *)'  --> Created slc_pu matrix (Φ_dot, U_tilde)'
        write(SLlogunit, *)'  --> Created slc_ut matrix (U_dot, θ_tilde)'
        write(SLlogunit, *)'  --> Created slc_pp matrix (Φ_dot, Φ_tilde)'
        write(SLlogunit, *)'  --> Created slc_pt matrix (Φ_dot, θ_tilde)'
        write(SLlogunit, *)'  --> Created slc_up matrix (U_dot, Φ_tilde)'
        write(SLlogunit, *)'  --> Created sea level (nodalsl) vector'
    endif



    write(SLlogunit,*)'  ✓ Prepared sea level. '
    write(SLlogunit,*)

    return 
end subroutine prepare_sea_level




subroutine set_original_sea_ice_level()
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
    nodalice0 = 0.0_kreal

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



    ! Set the initial ice levels: 




    deallocate(ds_quad4, gll_weight,lag_gll, dlag_gll)
end subroutine set_original_sea_ice_level

! ################### END INITIAL SETUP FUNCTIONS  #####################






subroutine update_ocean_function(u, errcode, errtag, use_s0)
! Routine checks each GLL point on the surface to see if it is part of the ocean set
! see Crawford et al 2018, eqn 31-32.
! Set contains any nodes in which rho_w * theta > rho_i * I 
! Bit of an issue here because we need the values of theta, I not their rates (time derivs)
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


write(SLlogunit,*)'WARNING: NEED TO ACCURATELY IMPLEMENT OCEAN SET/OCEAN FUNCTION BASED ON ICE CONDITION - see Crawford et al 2018, eqn 31'

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
    write(*, *)'ERROR: UPDATING OCEAN FUNCTION NOT IMPLEMENTED FOR CALCULATIONS YET'
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
    if(SLarea.lt.ZERO .or. SLarea.eq.ZERO)then 
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

    integer                        :: i_elmtfs        ! loops
    real(kind=kreal)               :: detjac2d        ! 2d jacobian
    integer                        :: iface           ! face ID for elmt 
    integer                        :: i_elmt          ! face ID for elmt 
    integer                        :: nfgll           ! ngll on 2D face
    real(kind=kreal), allocatable  :: gw(:)           ! GLL weights 2D
    real(kind=kreal), allocatable  :: dshape4(:,:,:)
    real(kind=kreal)               :: coord(ndim,4), face_normal(3),& 
                                    dx_dxi(NDIM), dx_deta(NDIM)
    integer :: num4(4), gid, phi_ind, u_ind, j,k, abg, xyg, gid_abg, gid_xyg, i_dim
    real(kind=kreal) :: theta_tf, pi_2d_abg, pi_2d_xyg, area_inv, &
                        ival, iival, rho_over_g, phi_tf, u_tf(NDIM)

    real(kind=kreal) :: g0abg, grav_abgj, Cabg, utfj, g0xyg, Cxyg, v1,rho_Ag


    ! Code
    allocate(gw(maxngll2d))
    allocate(dshape4(2,4,maxngll2d))

    write(SLlogunit,*)
    write(SLlogunit,*)'Calculating LHS matrices'

    ! Test functions. 
    theta_tf = ONE
    u_tf = ONE
    phi_tf = ONE
    ! Inverse area
    area_inv = ONE/SLarea

    ! Loop through each element on the free surface 
    do i_elmtfs = 1, nelmt_fs
        ! Get details of face
        call get_fs_details(i_elmtfs, iface, nfgll, gw, dshape4)
        num4   = gnum4_fs(:, i_elmtfs)
        coord  = g_coord(:,num4)

       
        do abg = 1, nfgll ! ABG
            gid_abg = gnum_fs(abg, i_elmtfs) ! Global ID of node ABG

            ! Get Jacobian_2d x weights for ABG 
            dx_dxi  = matmul(coord,dshape4(1,:,abg))
            dx_deta = matmul(coord,dshape4(2,:,abg))
            face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
            face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
            face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)
            pi_2d_abg   =  gw(abg) * sqrt(dot_product(face_normal,face_normal)) ! Weights*jacw

            
            ! Get the values here because they are repeated lots 
            g0abg      = g0_nodal(gid_abg)      ! g0 abg 
            Cabg       = oceanf(i_elmtfs, abg)  ! Ocean func abg


            ! Non-coupling Kmat terms (ie SL integral)
            QSL(abg, i_elmtfs) = QSL(abg, i_elmtfs) - (theta_tf * pi_2d_abg * g0abg * rho_water)


            ! Factor of rho/g outside of integral 
            rho_over_g =  (-rho_water/g0abg)       ! -rho/g
            rho_Ag     =  (rho_over_g / SLarea)    ! -rho/(g*Area)


            ! UPDATE THE ABG-ABG INDICES: 
            ! Phi_dot theta_tilde
            slc_pt(abg, abg, i_elmtfs) = slc_pt(abg, abg, i_elmtfs) + (g0abg * pi_2d_abg * theta_tf  * rho_over_g)
            ! Phi_dot phi_tilde
            slc_pp(abg, abg, i_elmtfs) = slc_pp(abg, abg, i_elmtfs) + (phi_tf * Cabg * pi_2d_abg  * rho_over_g)
            

            do j=1,NDIM
                grav_abgj = grav0_nodal(j, gid_abg)

                ! Phi_dot u_tilde
                slc_pu(j, abg, abg, i_elmtfs) = slc_pu(j, abg, abg, i_elmtfs) + (Cabg * pi_2d_abg *  u_tf(j) * grav_abgj * rho_over_g) 

                ! u_dot theta_tilde
                slc_ut(j, abg, abg, i_elmtfs) = slc_ut(j, abg, abg, i_elmtfs) + (pi_2d_abg * g0abg * theta_tf * grav_abgj * rho_over_g)

                ! u_dot phi_tilde
                slc_up(j, abg, abg, i_elmtfs) = slc_up(j, abg, abg, i_elmtfs) + (Cabg * pi_2d_abg * phi_tf * grav_abgj * rho_over_g)

                do k=1,NDIM
                    ! u_dot u_phi 
                    slc_uu(j, abg, k, abg, i_elmtfs) = slc_uu(j, abg, k, abg, i_elmtfs) + (Cabg * pi_2d_abg * grav_abgj *  u_tf(k) * grav0_nodal(k, gid_abg) * rho_over_g)
                enddo !k
            enddo  ! j 





            ! Diagonal+non-diagonal components of Kmat coupling SL 
            do xyg = 1, nfgll !XYG
                gid_xyg = gnum_fs(xyg, i_elmtfs) ! Global ID of node XYG

                ! Get Jacobian_2d x weights for XYG 
                dx_dxi  = matmul(coord,dshape4(1,:,xyg))
                dx_deta = matmul(coord,dshape4(2,:,xyg))
                face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3) 
                face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
                face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

                
                pi_2d_xyg  =  gw(xyg)*sqrt(dot_product(face_normal,face_normal)) ! Weights*jacw
                g0xyg      =  g0_nodal(gid_xyg)      ! g0 abg 
                Cxyg       =  oceanf(i_elmtfs, xyg)  ! Ocean func abg



                ! Note here that ABG is the index of the variable
                ! and XYG is the test function so when we assemble the 
                ! matrix, ABG should be the column index 

                slc_pt(abg, xyg, i_elmtfs) = slc_pt(abg, xyg, i_elmtfs) - (g0xyg * theta_tf * pi_2d_abg * Cabg * pi_2d_xyg *rho_Ag )

                slc_pp(abg, xyg, i_elmtfs) = slc_pp(abg, xyg, i_elmtfs) - (Cxyg * phi_tf * pi_2d_abg * Cabg * pi_2d_xyg * rho_Ag)


                do j=1,NDIM
                    v1 = pi_2d_abg * Cabg * grav0_nodal(j, gid_abg) * pi_2d_xyg 

                    slc_pu(j, abg, xyg, i_elmtfs) = slc_pu(j, abg, xyg, i_elmtfs) - (pi_2d_abg * Cabg * pi_2d_xyg * Cxyg * u_tf(j) *  grav0_nodal(j, gid_xyg) * rho_Ag )

                    slc_ut(j, abg, xyg, i_elmtfs) = slc_ut(j, abg, xyg, i_elmtfs) - (v1 * g0xyg * theta_tf * rho_Ag) 

                    slc_up(j, abg, xyg, i_elmtfs) = slc_up(j, abg, xyg, i_elmtfs) - (v1 * Cxyg * phi_tf * rho_Ag)

                    do k = 1, NDIM 
                        slc_uu(j, abg, k, xyg, i_elmtfs) = slc_uu(j, abg, k, xyg, i_elmtfs)  - (v1 * Cxyg * u_tf(k) *  grav0_nodal(k, gid_xyg) * rho_Ag)
                    enddo ! k 
                enddo ! j

            enddo! xyg
            
        enddo! abg
    enddo! i_elmtfs

    write(SLlogunit,*)'  ✓ Calculated LHS'

    deallocate(gw)
    deallocate(dshape4)
end subroutine calc_SL_LHS




end module